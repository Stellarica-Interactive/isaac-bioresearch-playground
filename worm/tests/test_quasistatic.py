"""Quasi-static body mechanics.

The test that matters is
:func:`test_travel_does_not_depend_on_the_timestep`. The solver this replaces
fails it: the same scripted gait covers 7.216 BL at 240 Hz and 1.442 BL at 60 Hz,
because a segment's velocity relaxes in 2.35 ms against a 4.17 ms step, so drag
erases momentum faster than the simulator advances. Everything else here guards
properties that would make a locomotion number meaningless if they broke.
"""

from __future__ import annotations

import numpy as np
import pytest

from worm.body.drag import DragParameters
from worm.body.geometry import BodyPlan
from worm.body.quasistatic import QuasiStaticBody

PLAN = BodyPlan()


def _wave_torques(body: QuasiStaticBody, t: float, *, amplitude: float) -> np.ndarray:
    """A travelling wave of joint torque: one wavelength, half a hertz."""
    joints = np.arange(body.plan.n_joints)
    phase = 2.0 * np.pi * (0.5 * t - joints / body.plan.n_joints)
    return amplitude * np.sin(phase)


def _run(
    *,
    hz: float,
    seconds: float,
    amplitude: float = 2.0e-5,
    ratio: float | None = None,
) -> QuasiStaticBody:
    kwargs = {} if ratio is None else {"ratio": ratio}
    body = QuasiStaticBody(PLAN, drag=DragParameters(**kwargs))
    dt = 1.0 / hz
    for step in range(int(seconds * hz)):
        body.step(_wave_torques(body, step * dt, amplitude=amplitude), dt_s=dt)
    return body


def _distance_bl(body: QuasiStaticBody, start: np.ndarray) -> float:
    moved = body.segment_centres().mean(axis=0) - start
    return float(np.linalg.norm(moved)) / PLAN.total_length_m


def test_the_force_balance_is_satisfied() -> None:
    """Nothing pushes the animal but the ground it pushes against.

    A non-zero residual means momentum is entering from somewhere, which is
    exactly the failure that made a frozen body appear to rotate for an hour
    (model_assumptions 5AE).
    """
    body = QuasiStaticBody(PLAN)
    for step in range(50):
        torques = _wave_torques(body, step / 60.0, amplitude=2.0e-5)
        force, torque = body.residual_force_and_torque(torques)
        assert force < 1e-12, f"net force {force:.3e} N"
        assert abs(torque) < 1e-12, f"net torque {torque:.3e} N m"
        body.step(torques, dt_s=1.0 / 60.0)


def test_travel_does_not_depend_on_the_timestep() -> None:
    """The property the inertial solver does not have.

    Dropping inertia removes the stability constraint entirely, so the step size
    only has to resolve the change in shape. Three rates spanning 8x should agree
    to a few per cent; the solver this replaces differs fivefold between two of
    them.
    """
    start = QuasiStaticBody(PLAN).segment_centres().mean(axis=0)
    distances = {hz: _distance_bl(_run(hz=hz, seconds=8.0), start) for hz in (60, 240, 480)}
    assert min(distances.values()) > 0.01, f"barely moved at all: {distances}"
    spread = max(distances.values()) - min(distances.values())
    assert spread < 0.05 * max(distances.values()), f"timestep-dependent: {distances}"


def test_a_travelling_wave_moves_the_body() -> None:
    start = QuasiStaticBody(PLAN).segment_centres().mean(axis=0)
    assert _distance_bl(_run(hz=240, seconds=8.0), start) > 0.05


def test_isotropic_drag_produces_no_net_travel() -> None:
    """The mechanism, isolated.

    Anisotropy is the whole reason undulation becomes thrust: with equal
    coefficients along and across the body, every push forward is matched by an
    equal slip sideways. A solver that still travels under isotropic drag is
    getting its motion from somewhere other than the medium.
    """
    start = QuasiStaticBody(PLAN).segment_centres().mean(axis=0)
    isotropic = _distance_bl(_run(hz=240, seconds=8.0, ratio=1.0), start)
    anisotropic = _distance_bl(_run(hz=240, seconds=8.0), start)
    assert isotropic < 0.1 * anisotropic, f"{isotropic=:.4f} {anisotropic=:.4f}"


def test_no_drive_means_no_motion() -> None:
    """With no muscle torque the body must sit exactly still -- the control the
    inertial solver passes too, and the one that proves nothing drifts."""
    body = QuasiStaticBody(PLAN)
    start = body.segment_centres().mean(axis=0)
    for _ in range(600):
        body.step(np.zeros(PLAN.n_joints), dt_s=1.0 / 60.0)
    assert _distance_bl(body, start) < 1e-9
    assert abs(body.heading) < 1e-12


def test_joints_respect_the_shared_limit() -> None:
    """The same 60 degree ceiling the Isaac articulation enforces. A body free to
    fold further would not be comparable with the committed runs."""
    body = QuasiStaticBody(PLAN)
    for _ in range(2000):
        body.step(np.full(PLAN.n_joints, 1e-3), dt_s=1.0 / 60.0)
    assert np.all(np.abs(body.joint_angles) <= np.radians(60.0) + 1e-12)


def test_the_body_keeps_its_length() -> None:
    """Segment spacing is kinematic here, so length cannot drift -- unlike an
    articulation, where a solver can stretch joints under load."""
    body = _run(hz=240, seconds=5.0)
    spacing = np.linalg.norm(np.diff(body.nodes(), axis=0), axis=1)
    assert np.allclose(spacing, body.segment_length_m, rtol=1e-12)


@pytest.mark.parametrize("ratio", (3.0, 10.0, 20.0))
def test_more_anisotropy_moves_further(ratio: float) -> None:
    """Monotonic in the one drag parameter with published measurements behind it
    (Rabets et al. 2014 measured 3-10; the committed default is 20)."""
    start = QuasiStaticBody(PLAN).segment_centres().mean(axis=0)
    distance = _distance_bl(_run(hz=240, seconds=8.0, ratio=ratio), start)
    baseline = _distance_bl(_run(hz=240, seconds=8.0, ratio=1.0), start)
    assert distance > baseline


# -- how the body turns while it crawls ------------------------------------
#
# Watched in the viewport, the scripted gait looked like it was rocking the whole
# animal in rotation rather than holding a course. It is not: `heading` is
# segment 0's direction and segment 0 is the head, so the large number is a head
# swing. These guard both halves of that -- that the course really does hold, and
# that the head/body distinction stays visible in the API.

SETTLE_CYCLES = 6
GAIT_HZ = 0.5


def _wave_n(body: QuasiStaticBody, t: float, *, amplitude: float, waves: float) -> np.ndarray:
    """``_wave_torques`` with the number of wavelengths on the body exposed."""
    joints = np.arange(body.plan.n_joints)
    phase = 2.0 * np.pi * (GAIT_HZ * t - waves * joints / body.plan.n_joints)
    return amplitude * np.sin(phase)


def _turning(*, waves: float, hz: float = 60.0, seconds: float = 24.0):
    """Per-cycle head swing, body swing and course drift, after settling.

    The first cycles are discarded: the body starts straight and winding the wave
    onto it turns the frame once, permanently. That is a reorientation, not a
    swing, and counting it is how this measurement was first reported wrongly.
    """
    body = QuasiStaticBody(PLAN, drag=DragParameters())
    dt = 1.0 / hz
    heads, bodies = [], []
    for step in range(int(seconds * hz)):
        body.step(_wave_n(body, step * dt, amplitude=2.0e-5, waves=waves), dt_s=dt)
        heads.append(body.heading)
        mean = body._directions().mean(axis=0)
        bodies.append(np.arctan2(mean[1], mean[0]))
    heads = np.unwrap(np.array(heads))
    bodies = np.unwrap(np.array(bodies))

    per_cycle = int(hz / GAIT_HZ)
    first = SETTLE_CYCLES * per_cycle
    head_swings, body_swings, levels = [], [], []
    for start in range(first, len(heads) - per_cycle + 1, per_cycle):
        window = heads[start : start + per_cycle]
        head_swings.append(window.max() - window.min())
        levels.append(window.mean())
        other = bodies[start : start + per_cycle]
        body_swings.append(other.max() - other.min())
    drift = float(np.mean(np.diff(levels))) if len(levels) > 1 else 0.0
    return (
        float(np.degrees(np.mean(head_swings))),
        float(np.degrees(np.mean(body_swings))),
        float(np.degrees(drift)),
    )


def test_the_course_does_not_drift() -> None:
    """The swing is an oscillation about a fixed heading, not a slow turn.

    This is what separates "waves its head while crawling straight" from "crawls
    in a circle", which is what it looked like it was doing. Measured drift is
    0.00 deg per cycle; the bound is loose on purpose, since a real defect here
    would be degrees, not hundredths.
    """
    _, _, drift = _turning(waves=1.0)
    assert abs(drift) < 0.5, f"course drifts {drift:+.3f} deg per cycle"


def test_the_head_swings_several_times_more_than_the_body() -> None:
    """``heading`` is the head's direction, and must not be read as the body's.

    Measured: 61.8 deg of head against 13.2 deg of body axis per cycle, which is
    the right shape for a crawling worm. A solver whose two numbers converged
    would be rotating the whole animal.
    """
    head, axis, _ = _turning(waves=1.0)
    assert head > 3.0 * axis, f"head {head:.1f} deg vs body {axis:.1f} deg"


def test_more_waves_on_the_body_means_less_yaw() -> None:
    """The signature of zero-net-torque counter-rotation.

    Nothing external torques the animal, so the frame must counter-rotate as the
    bending wave carries mass around. One wavelength does that asymmetrically --
    front and back of the wave are not mirror images -- and more wavelengths
    cancel more completely. If the yaw did not fall, it would be the solver
    rather than the mechanics.
    """
    one, _, _ = _turning(waves=1.0)
    two, _, _ = _turning(waves=2.0)
    assert two < 0.8 * one, f"1 wave {one:.1f} deg, 2 waves {two:.1f} deg"


# -- joints at their limit --------------------------------------------------
#
# A joint limit is an internal constraint, equal and opposite across the joint,
# so it cannot torque the animal as a whole. Solving the joints freely and then
# clipping breaks that: the rigid rates are solved jointly with qdot, so
# discarding part of qdot leaves a velocity field that no longer balances, and
# the leftover rotation accumulates every step. `test_joints_respect_the_shared
# _limit` above drives at exactly the torque that triggers it and passes, because
# it checks the angles and never looks at the heading.

SATURATING_TORQUE = 1.0e-3


def _settle(torque: float, *, seconds: float = 30.0, hz: float = 240.0):
    """Hold a uniform drive and report the heading rate over the last third."""
    body = QuasiStaticBody(PLAN)
    drive = np.full(PLAN.n_joints, torque)
    dt = 1.0 / hz
    headings = []
    for _ in range(int(seconds * hz)):
        body.step(drive, dt_s=dt)
        headings.append(body.heading)
    tail = np.array(headings[-int(seconds * hz / 3) :])
    rate = float(np.degrees(tail[-1] - tail[0])) / (seconds / 3.0)
    at_limit = int(np.sum(np.abs(body.joint_angles) >= np.radians(60.0) - 1e-9))
    return body, rate, at_limit


def test_a_saturated_body_does_not_spin() -> None:
    """A body held in a static shape must stop, and did not.

    A uniform drive bends the body into an arc and holds it there. That shape
    generates no thrust -- anisotropic drag converts *undulation* into force and
    there is no undulation -- so once the arc stops changing the body must stop.
    Before the joints were locked out of the solve it spun at -26814 deg/s,
    sustained, for as long as the run lasted.
    """
    body, rate, at_limit = _settle(SATURATING_TORQUE)
    assert at_limit == PLAN.n_joints, f"only {at_limit} joints saturated; test is not testing it"
    assert abs(rate) < 0.1, f"spins at {rate:+.1f} deg/s with {at_limit} joints at the limit"


def test_the_force_balance_holds_with_joints_at_their_limit() -> None:
    """The residual, re-asserted in the regime where clipping used to break it.

    `test_the_force_balance_is_satisfied` runs an unsaturated gait, so it never
    reached this path -- which is why a solver that spun at 74 revolutions per
    second passed every test in this file.
    """
    body, _, at_limit = _settle(SATURATING_TORQUE, seconds=10.0)
    assert at_limit == PLAN.n_joints
    drive = np.full(PLAN.n_joints, SATURATING_TORQUE)
    force, torque = body.residual_force_and_torque(drive)
    assert force < 1e-12, f"net force {force:.3e} N"
    assert abs(torque) < 1e-12, f"net torque {torque:.3e} N m"


def test_a_joint_at_its_limit_still_moves_back_off_it() -> None:
    """Locking must be one-sided, or the body would stick permanently.

    A joint is held only against the drive that pins it. Reverse the drive and it
    has to come straight off the limit, otherwise the first strong contraction of
    a run would freeze that joint for good.
    """
    body, _, at_limit = _settle(SATURATING_TORQUE, seconds=10.0)
    assert at_limit == PLAN.n_joints
    pinned = body.joint_angles.copy()
    for _ in range(240):
        body.step(np.full(PLAN.n_joints, -SATURATING_TORQUE), dt_s=1.0 / 240.0)
    assert np.all(body.joint_angles < pinned - 1e-6), "joints stayed pinned under reversed drive"


def test_a_straight_body_is_clear_of_itself() -> None:
    """The self-distance metric, on the one configuration with an obvious answer.

    Non-neighbouring segments on a straight body are at least a full segment
    apart, so the gap is positive. Immediate neighbours share a node and always
    overlap, which is why the metric skips them -- without that it would report a
    straight worm as interpenetrating.
    """
    assert QuasiStaticBody(PLAN).min_self_distance_m() > 0.0


def test_a_closed_coil_interpenetrates() -> None:
    """What the solver permits, stated as a measurement rather than assumed.

    The geometric argument that a coil clears itself is about *opposite sides* of
    the loop, 32 mm apart on a 6.5 mm thick body, and it answers the wrong
    question. What binds is the nose meeting the tail, where the body is thinnest
    -- 1.09 mm of radius each, so 2.18 mm of clearance and no more. Curl the body
    into a closed loop and it passes through itself. There is no contact model to
    stop it; this records that, so adding one has a test to turn.
    """
    body = QuasiStaticBody(PLAN)
    body.joint_angles = np.full(PLAN.n_joints, 2.0 * np.pi / PLAN.n_joints)
    nodes = body.nodes()
    extent = float(np.linalg.norm(nodes[-1] - nodes[0])) / PLAN.total_length_m
    assert extent < 0.1, f"not actually coiled: extent {extent:.2f}"
    assert body.min_self_distance_m() < 0.0, "a closed coil should interpenetrate"


# -- a drag law with a threshold --------------------------------------------
#
# Linear drag inverts to v = F/c, so no force is too small to move the body and
# nothing is ever stuck. Observed in the viewport as a body holding a stable
# shape while sliding steadily across the ground. Two departures from linearity
# are available, both off by default so every committed measurement is
# reproducible: a sublinear exponent and a yield force. These check that the
# nonlinear solve is a solve rather than an approximation that happens to look
# plausible. See model_assumptions 5AH.

SUBLINEAR = DragParameters(exponent=0.6)
YIELDING = DragParameters(yield_force=1.0e-4)


def test_the_default_drag_law_is_linear() -> None:
    """Everything measured in this project assumed it, so it has to stay."""
    assert DragParameters().is_linear
    assert not SUBLINEAR.is_linear
    assert not YIELDING.is_linear
    assert not DragParameters(exponent=0.6, yield_force=1.0e-4).is_linear


def test_the_sublinear_law_stiffens_below_its_reference_speed() -> None:
    """That is the whole mechanism: a groove resists more as you slow down."""
    reference = SUBLINEAR.reference_speed
    slow, fast = np.array([reference / 100.0]), np.array([reference * 100.0])
    assert SUBLINEAR.speed_factor(slow)[0] > 1.0
    assert SUBLINEAR.speed_factor(fast)[0] < 1.0
    assert np.isclose(SUBLINEAR.speed_factor(np.array([reference]))[0], 1.0)
    # The linear law must be exactly inert, not merely close to it.
    assert np.all(DragParameters().speed_factor(slow) == 1.0)


def test_the_yield_force_saturates() -> None:
    """A yield term is a force, not a coefficient: past the regularisation speed
    it must stop growing, or it is just more viscosity under another name."""
    speed = np.array([YIELDING.yield_speed * 1e3, YIELDING.yield_speed * 1e4])
    force = YIELDING.yield_coefficient(speed) * speed
    assert np.allclose(force, YIELDING.yield_force, rtol=1e-3)
    assert np.all(DragParameters().yield_coefficient(speed) == 0.0)


@pytest.mark.parametrize("drag", (SUBLINEAR, YIELDING), ids=("sublinear", "yielding"))
def test_the_force_balance_holds_under_a_nonlinear_law(drag: DragParameters) -> None:
    """The fixed point has to actually converge.

    The balance is solved with coefficients that depend on the velocities it
    produces, so it is closed by iteration. An iteration that stopped early
    would leave a velocity field that does not balance -- the same class of
    defect as clipping the joint limits, which spun the body at 26814 deg/s
    while every other test passed.
    """
    body = QuasiStaticBody(PLAN, drag=drag)
    worst = 0.0
    for step in range(40):
        torques = _wave_torques(body, step / 60.0, amplitude=5.0e-5)
        worst = max(worst, body.residual_fraction(torques))
        body.step(torques, dt_s=1.0 / 60.0)
    # Relative, not absolute: an absolute bound on a nonlinear solve just encodes
    # the drag magnitude. Written first as the linear path's 1e-12 N, which
    # failed at 1.6e-8 N -- correctly, and the bound was wrong rather than the
    # solver. At the committed tolerance the achieved fraction is 5e-9, so this
    # has two orders of headroom and would still catch a stalled iteration.
    assert worst < 1.0e-7, f"{worst:.3e} of the drag force fails to cancel"


@pytest.mark.parametrize("drag", (SUBLINEAR, YIELDING), ids=("sublinear", "yielding"))
def test_a_nonlinear_law_still_needs_anisotropy_to_travel(drag: DragParameters) -> None:
    """The mechanism check, repeated for each law.

    Thrust must still come from the medium resisting sideways motion more than
    lengthwise motion. A nonlinear law that travelled under isotropic drag would
    be generating motion out of its own nonlinearity.
    """
    start = QuasiStaticBody(PLAN).segment_centres().mean(axis=0)
    isotropic = QuasiStaticBody(
        PLAN, drag=DragParameters(ratio=1.0, exponent=drag.exponent, yield_force=drag.yield_force)
    )
    anisotropic = QuasiStaticBody(PLAN, drag=drag)
    # Four seconds is two gait cycles, which separates the two cases by an order
    # of magnitude. Kept short on purpose: the nonlinear solve runs up to 80
    # passes per step, so simulated seconds here are expensive in a way they are
    # not for the linear law.
    for step in range(int(4.0 * 120)):
        for body in (isotropic, anisotropic):
            body.step(_wave_torques(body, step / 120.0, amplitude=5.0e-5), dt_s=1.0 / 120.0)
    flat = _distance_bl(isotropic, start)
    curved = _distance_bl(anisotropic, start)
    assert flat < 0.2 * curved, f"{flat=:.4f} {curved=:.4f}"


@pytest.mark.parametrize("drag", (SUBLINEAR, YIELDING), ids=("sublinear", "yielding"))
def test_a_nonlinear_law_still_sits_still_with_no_drive(drag: DragParameters) -> None:
    """A threshold must not become a source. With no muscle torque the body has
    to stay exactly where it is, as it does under the linear law."""
    body = QuasiStaticBody(PLAN, drag=drag)
    start = body.segment_centres().mean(axis=0)
    for _ in range(120):
        body.step(np.zeros(PLAN.n_joints), dt_s=1.0 / 60.0)
    assert _distance_bl(body, start) < 1e-9
    assert abs(body.heading) < 1e-12


# -- per-segment orientation ------------------------------------------------
#
# `heading` is segment 0's angle and nothing else's. Driving an articulation
# whose root link is the MIDDLE of the chain with the head's angle displaced the
# rendered body 47.9 mm and rotated it by a figure that oscillated at the gait
# frequency, which is what a user spent a session reporting as the animal
# rotating. See model_assumptions 5AJ.
#
# Note what these tests cannot do: the defect was in the articulation, and this
# suite cannot launch Isaac. `tools/check_kinematic_drive.py` is the real guard
# and has to be run by hand. These cover only the solver side.


def test_segment_angles_start_at_the_heading() -> None:
    """Segment 0 is the head, so its angle is the heading by definition."""
    body = QuasiStaticBody(PLAN)
    body.heading = 0.7
    body.joint_angles = np.linspace(-0.2, 0.2, PLAN.n_joints)
    angles = body.segment_angles()
    assert angles.shape == (PLAN.n_segments,)
    assert angles[0] == pytest.approx(0.7)


def test_segment_angles_diverge_from_the_heading_along_the_body() -> None:
    """The whole point: no other segment's angle is the heading.

    A bent body's segments point in different directions, so anything that needs
    one segment's orientation has to ask for that segment. Writing `heading` for
    all of them is the 5AJ defect.
    """
    body = QuasiStaticBody(PLAN)
    body.heading = 0.0
    body.joint_angles = np.full(PLAN.n_joints, np.radians(5.0))
    angles = np.degrees(body.segment_angles())
    assert angles[0] == pytest.approx(0.0)
    assert angles[-1] == pytest.approx(5.0 * PLAN.n_joints)
    # Monotonic under a uniform bend, and spanning a large angle: the middle
    # segment is already 55 degrees off the head on a gently curled body.
    assert np.all(np.diff(angles) > 0)
    assert abs(angles[PLAN.n_segments // 2]) > 50.0


def test_segment_angles_agree_with_the_direction_vectors() -> None:
    """One source of truth. ``_directions`` is built from this, and the force
    balance uses ``_directions``, so a disagreement would mean the solver and
    anything driving a display were working from different bodies."""
    body = QuasiStaticBody(PLAN)
    body.heading = -1.1
    rng = np.random.default_rng(0)
    body.joint_angles = rng.uniform(-0.3, 0.3, PLAN.n_joints)
    angles = body.segment_angles()
    expected = np.stack([np.cos(angles), np.sin(angles)], axis=1)
    assert np.allclose(body._directions(), expected)


# -- self-contact -----------------------------------------------------------
#
# The connectome curls the body until `clear` reaches -1.39 mm: it passes through
# itself by a little over a tail diameter. `self_contact` adds stiff drag along
# the line between two non-neighbouring segments that have closed to touching,
# acting only on approach. Off by default, because every number in
# model_assumptions was measured without it.
#
# The band width is bounded above by a measurement rather than chosen freely: the
# scripted gait's closest self-approach is 1.57 mm, so a wider band would engage
# during ordinary crawling and make every locomotion figure contingent on this
# model. See model_assumptions 5AF.7.

CURLING_TORQUE = 6.0e-5


def _worst_clearance(*, contact: bool, seconds: float = 12.0) -> float:
    """Closest the body comes to itself while curling, in metres."""
    body = QuasiStaticBody(PLAN, self_contact=contact)
    drive = np.full(PLAN.n_joints, CURLING_TORQUE)
    worst = 1.0
    for _ in range(int(seconds * 240)):
        body.step(drive, dt_s=1.0 / 240.0)
        worst = min(worst, body.min_self_distance_m())
    return worst


def test_without_contact_the_body_folds_through_itself() -> None:
    """The permission being removed, asserted so the fix has something to fix.

    A uniform drive curls the body until non-neighbouring segments overlap. The
    solver knows only drag, so nothing stops them.
    """
    assert _worst_clearance(contact=False) < -1.0e-3, "expected the body to interpenetrate"


def test_self_contact_keeps_the_body_clear_of_itself() -> None:
    """The same curl, with contact on, must not interpenetrate at all.

    Measured over a longer curl: worst clearance goes from -5.80 mm to +0.79 mm.
    Asserting merely "less overlap" would pass on the first version of this,
    which improved -6.29 mm to -5.18 mm and did not work.
    """
    assert _worst_clearance(contact=True) > 0.0


def test_self_contact_does_not_disturb_the_gait() -> None:
    """A contact model that changed the gait would make every locomotion number
    in model_assumptions contingent on two arbitrary constants.

    It does not, and that is a property of the band width rather than luck: the
    scripted gait never comes within 1.57 mm of itself, which is outside the
    band, so no contact ever engages. Measured, the speed is identical to four
    decimal places with contact on and off.
    """
    start = QuasiStaticBody(PLAN).segment_centres().mean(axis=0)
    distances = []
    for contact in (False, True):
        body = QuasiStaticBody(PLAN, drag=DragParameters(), self_contact=contact)
        for step in range(int(8.0 * 120)):
            body.step(_wave_torques(body, step / 120.0, amplitude=5.0e-5), dt_s=1.0 / 120.0)
        distances.append(_distance_bl(body, start))
    assert distances[0] == pytest.approx(distances[1], rel=1e-9), (
        f"contact changed the gait: {distances[0]:.6f} vs {distances[1]:.6f} BL"
    )
