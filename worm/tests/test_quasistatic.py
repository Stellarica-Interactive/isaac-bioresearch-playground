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
