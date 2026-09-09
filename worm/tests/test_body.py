"""Body plan, muscles and drag — all testable without launching a simulator.

That separation is the point: the biomechanics are where the modelling decisions
live, so they should be checkable in milliseconds rather than only visible once
something is wobbling in a viewport.

The drag tests matter most. An undulating chain on an isotropic surface goes
nowhere, so if the anisotropy is wrong the worm will wriggle in place and the
symptom will look like a neural problem.
"""

from __future__ import annotations

import numpy as np
import pytest

from worm.body.drag import DragParameters, GroundDrag
from worm.body.geometry import (
    DEFAULT_SEGMENTS,
    QUADRANTS,
    BodyPlan,
    muscle_map,
    muscle_name,
    parse_muscle,
    segment_of_muscle,
)
from worm.body.muscles import MuscleModel, MuscleParameters, sine_wave_drive
from worm.importers.naming import body_wall_muscle_ids


@pytest.fixture(scope="module")
def plan() -> BodyPlan:
    return BodyPlan()


class TestBodyPlan:
    def test_segment_count_matches_the_muscle_rows(self, plan: BodyPlan) -> None:
        """24 is derived from anatomy, not chosen for convenience."""
        assert plan.n_segments == DEFAULT_SEGMENTS == 24
        assert plan.n_joints == 23

    def test_segments_tile_the_body_without_gaps(self, plan: BodyPlan) -> None:
        segments = plan.segments()
        assert len(segments) == plan.n_segments
        total = sum(s.length_m for s in segments)
        assert total == pytest.approx(plan.total_length_m)
        # Not strict: the tail is one shorter by construction.
        for a, b in zip(segments, segments[1:], strict=False):
            gap = b.centre_x_m - a.centre_x_m
            assert gap == pytest.approx(plan.segment_length_m)

    def test_body_tapers_at_both_ends(self, plan: BodyPlan) -> None:
        radii = [s.radius_m for s in plan.segments()]
        middle = radii[len(radii) // 2]
        assert radii[0] < middle
        assert radii[-1] < middle
        assert max(radii) == pytest.approx(middle, rel=0.05)

    def test_all_radii_are_positive(self, plan: BodyPlan) -> None:
        assert all(s.radius_m > 0 for s in plan.segments())

    def test_mass_is_physically_sensible(self, plan: BodyPlan) -> None:
        """A real adult is about 1 microgram; scaled 100x in length it is 1e6 times
        the volume, so about a gram."""
        assert 1e-4 < plan.total_mass_kg < 1e-2

    def test_scaling_changes_size_but_not_shape(self) -> None:
        small, big = BodyPlan(length_scale=1.0), BodyPlan(length_scale=100.0)
        ratio = big.total_length_m / small.total_length_m
        assert ratio == pytest.approx(100.0)
        for a, b in zip(small.segments(), big.segments(), strict=True):
            assert b.radius_m / a.radius_m == pytest.approx(100.0)


class TestMuscleMapping:
    @pytest.mark.parametrize(
        ("name", "expected"),
        [("MDL01", ("DL", 1)), ("MVR23", ("VR", 23)), ("MDR24", ("DR", 24))],
    )
    def test_parses_body_wall_muscles(self, name: str, expected: tuple[str, int]) -> None:
        assert parse_muscle(name) == expected

    @pytest.mark.parametrize("name", ["AVAL", "pm4VL", "vm1AL", "MANAL", "mc1V"])
    def test_rejects_everything_that_is_not_a_body_wall_muscle(self, name: str) -> None:
        """Pharyngeal, vulval and anal muscles must not be swept into the body."""
        assert parse_muscle(name) is None

    def test_round_trip(self) -> None:
        for quadrant in QUADRANTS:
            for row in (1, 7, 24):
                assert parse_muscle(muscle_name(quadrant, row)) == (quadrant, row)

    def test_at_24_segments_the_mapping_is_one_to_one(self, plan: BodyPlan) -> None:
        """The reason for choosing 24: every muscle row gets its own segment."""
        for row in range(1, 25):
            assert segment_of_muscle(muscle_name("DL", row), plan) == row - 1

    def test_real_muscle_names_all_map(self, plan: BodyPlan) -> None:
        mapped = muscle_map(sorted(body_wall_muscle_ids()), plan)
        assert sum(len(v) for v in mapped.values()) == 95
        # 24 + 24 + 23 + 24: the animal really is missing MVL24.
        assert sorted(len(v) for v in mapped.values()) == [23, 24, 24, 24]

    def test_every_segment_is_driven(self, plan: BodyPlan) -> None:
        mapped = muscle_map(sorted(body_wall_muscle_ids()), plan)
        for quadrant, segments in mapped.items():
            missing = set(range(plan.n_segments)) - set(segments)
            # Only the ventral-left quadrant has a gap, and it is a real one.
            assert missing in ({23}, set()), (quadrant, missing)


class TestMuscleModel:
    def test_starts_relaxed(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        assert np.all(m.activation == 0.0)
        assert np.allclose(m.joint_torques(), 0.0)

    def test_activation_lags_its_drive(self, plan: BodyPlan) -> None:
        """Muscle is slow. One time constant gets it 63% of the way, not all."""
        m = MuscleModel(plan)
        drive = np.ones_like(m.activation)
        m.step(drive, dt_ms=m.params.activation_tau_ms)
        assert np.allclose(m.activation, 1.0 - np.exp(-1.0), atol=1e-9)

    def test_activation_converges_and_stays_bounded(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        drive = np.ones_like(m.activation)
        for _ in range(500):
            m.step(drive, dt_ms=10.0)
        assert np.allclose(m.activation, 1.0, atol=1e-6)
        assert np.all(m.activation <= 1.0)

    def test_activation_is_stable_at_any_timestep(self, plan: BodyPlan) -> None:
        """Exact exponential step, so a huge dt saturates rather than diverging."""
        m = MuscleModel(plan)
        m.step(np.ones_like(m.activation), dt_ms=1e6)
        assert np.all(np.isfinite(m.activation))
        assert np.all(m.activation <= 1.0)

    def test_drive_is_clamped_to_the_unit_interval(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        m.step(np.full_like(m.activation, 5.0), dt_ms=1e6)
        assert np.all(m.activation <= 1.0 + 1e-9)
        m.reset()
        m.step(np.full_like(m.activation, -5.0), dt_ms=1e6)
        assert np.all(m.activation >= -1e-9)

    def test_wrong_shaped_drive_is_rejected(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        with pytest.raises(ValueError, match="quadrants x segments"):
            m.step(np.zeros(plan.n_segments), dt_ms=1.0)

    def test_dorsal_contraction_bends_dorsally(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        drive = np.zeros_like(m.activation)
        for q in ("DL", "DR"):
            drive[QUADRANTS.index(q)] = 1.0
        m.step(drive, dt_ms=1e6)
        assert np.all(m.joint_torques() > 0)

    def test_ventral_contraction_bends_the_other_way(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        drive = np.zeros_like(m.activation)
        for q in ("VL", "VR"):
            drive[QUADRANTS.index(q)] = 1.0
        m.step(drive, dt_ms=1e6)
        assert np.all(m.joint_torques() < 0)

    def test_co_contraction_produces_no_net_torque_but_real_stiffness(
        self, plan: BodyPlan
    ) -> None:
        """Pulling both sides equally stiffens the joint without bending it.

        Real behaviour, and it falls out of the antagonist model rather than being
        added to it.
        """
        m = MuscleModel(plan)
        m.step(np.ones_like(m.activation), dt_ms=1e6)
        assert np.allclose(m.joint_torques(), 0.0, atol=1e-12)
        assert np.all(m.stiffness() > 0)

    def test_passive_stiffness_opposes_bending(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        angles = np.full(plan.n_joints, 0.5)
        assert np.all(m.joint_torques(joint_angles=angles) < 0)

    def test_passive_damping_opposes_motion(self, plan: BodyPlan) -> None:
        m = MuscleModel(plan)
        velocities = np.full(plan.n_joints, 1.0)
        assert np.all(m.joint_torques(joint_velocities=velocities) < 0)

    def test_torque_scales_with_the_parameter(self, plan: BodyPlan) -> None:
        drive = np.zeros((len(QUADRANTS), plan.n_segments))
        drive[QUADRANTS.index("DL")] = 1.0
        weak = MuscleModel(plan, MuscleParameters(peak_torque_scale=1e-5))
        strong = MuscleModel(plan, MuscleParameters(peak_torque_scale=1e-4))
        for m in (weak, strong):
            m.step(drive, dt_ms=1e6)
        assert np.allclose(strong.joint_torques() / weak.joint_torques(), 10.0)

    def test_one_torque_per_joint(self, plan: BodyPlan) -> None:
        assert MuscleModel(plan).joint_torques().shape == (plan.n_joints,)


class TestScriptedWave:
    """The control. Not biology, and deleted from any biological experiment."""

    def test_shape_and_range(self, plan: BodyPlan) -> None:
        drive = sine_wave_drive(plan, t_ms=0.0)
        assert drive.shape == (len(QUADRANTS), plan.n_segments)
        assert np.all(drive >= 0.0) and np.all(drive <= 1.0)

    def test_a_muscle_can_only_pull(self, plan: BodyPlan) -> None:
        """Half-wave rectified: dorsal and ventral are never both driven."""
        for t in np.linspace(0, 4000, 40):
            drive = sine_wave_drive(plan, t_ms=float(t))
            dorsal = drive[[QUADRANTS.index(q) for q in ("DL", "DR")]].max(axis=0)
            ventral = drive[[QUADRANTS.index(q) for q in ("VL", "VR")]].max(axis=0)
            assert np.all(np.minimum(dorsal, ventral) < 1e-12)

    def test_it_is_a_wave_not_a_synchronous_twitch(self, plan: BodyPlan) -> None:
        """Every segment must not do the same thing at the same time, or the body
        curls and uncurls instead of propagating a bend."""
        m = MuscleModel(plan, MuscleParameters(activation_tau_ms=1.0))
        m.step(sine_wave_drive(plan, t_ms=250.0), dt_ms=50.0)
        torque = m.joint_torques()
        assert np.any(torque > 0) and np.any(torque < 0)

    def test_the_wave_travels_head_to_tail(self, plan: BodyPlan) -> None:
        """Direction matters: a tail-to-head wave drives backward locomotion."""

        def peak_segment(t_ms: float) -> int:
            drive = sine_wave_drive(plan, t_ms=t_ms)
            return int(np.argmax(drive[QUADRANTS.index("DL")]))

        # Sample within a single cycle so the peak does not wrap around.
        early, late = peak_segment(0.0), peak_segment(200.0)
        assert late > early

    def test_it_is_periodic(self, plan: BodyPlan) -> None:
        period_ms = 1000.0 / 0.5
        a = sine_wave_drive(plan, t_ms=137.0)
        b = sine_wave_drive(plan, t_ms=137.0 + period_ms)
        assert np.allclose(a, b)


class TestGroundDrag:
    """The model that turns undulation into forward motion.

    Without anisotropy the worm oscillates in place, and the failure looks like a
    neural problem rather than a mechanical one.
    """

    def straight_body(self, plan: BodyPlan) -> np.ndarray:
        return np.stack(
            [np.array([s.centre_x_m for s in plan.segments()]), np.zeros(plan.n_segments)],
            axis=1,
        )

    def test_tangent_of_a_straight_body_points_along_it(self, plan: BodyPlan) -> None:
        drag = GroundDrag(plan)
        tangent = drag.tangents(self.straight_body(plan))
        assert np.allclose(tangent, np.array([1.0, 0.0]))

    def test_tangents_are_unit_vectors(self, plan: BodyPlan) -> None:
        drag = GroundDrag(plan)
        pos = self.straight_body(plan)
        pos[:, 1] = 0.01 * np.sin(np.linspace(0, 4 * np.pi, plan.n_segments))
        assert np.allclose(np.linalg.norm(drag.tangents(pos), axis=1), 1.0)

    def test_sideways_motion_is_resisted_far_more_than_forward(
        self, plan: BodyPlan
    ) -> None:
        """The whole point. Ratio of resistance must equal the drag ratio."""
        drag = GroundDrag(plan, DragParameters(tangential=1.0, ratio=20.0))
        pos = self.straight_body(plan)
        along = np.tile(np.array([1.0, 0.0]), (plan.n_segments, 1))
        across = np.tile(np.array([0.0, 1.0]), (plan.n_segments, 1))

        f_along = np.linalg.norm(drag.forces(pos, along), axis=1)
        f_across = np.linalg.norm(drag.forces(pos, across), axis=1)
        assert np.allclose(f_across / f_along, 20.0)

    def test_drag_always_opposes_velocity(self, plan: BodyPlan) -> None:
        """It can only remove energy. A drag model that can add energy makes any
        locomotion result meaningless."""
        rng = np.random.default_rng(0)
        drag = GroundDrag(plan)
        pos = self.straight_body(plan)
        pos[:, 1] = 0.02 * rng.standard_normal(plan.n_segments)
        for _ in range(20):
            vel = rng.standard_normal((plan.n_segments, 2))
            power = np.einsum("ij,ij->i", drag.forces(pos, vel), vel)
            assert np.all(power <= 1e-12)

    def test_no_motion_no_force(self, plan: BodyPlan) -> None:
        drag = GroundDrag(plan)
        pos = self.straight_body(plan)
        assert np.allclose(drag.forces(pos, np.zeros_like(pos)), 0.0)

    def test_force_is_linear_in_velocity(self, plan: BodyPlan) -> None:
        drag = GroundDrag(plan)
        pos = self.straight_body(plan)
        vel = np.tile(np.array([0.3, 0.7]), (plan.n_segments, 1))
        assert np.allclose(drag.forces(pos, 2 * vel), 2 * drag.forces(pos, vel))

    def test_isotropic_drag_is_expressible_and_shows_why_it_fails(
        self, plan: BodyPlan
    ) -> None:
        """With ratio 1 the surface cannot distinguish directions at all, which is
        the configuration in which an undulating body goes nowhere."""
        drag = GroundDrag(plan, DragParameters(tangential=1.0, ratio=1.0))
        pos = self.straight_body(plan)
        along = np.tile(np.array([1.0, 0.0]), (plan.n_segments, 1))
        across = np.tile(np.array([0.0, 1.0]), (plan.n_segments, 1))
        assert np.allclose(
            np.linalg.norm(drag.forces(pos, along), axis=1),
            np.linalg.norm(drag.forces(pos, across), axis=1),
        )

    def test_degenerate_geometry_does_not_produce_nan(self, plan: BodyPlan) -> None:
        """Coincident segment centres must not divide by zero mid-simulation."""
        drag = GroundDrag(plan)
        pos = np.zeros((plan.n_segments, 2))
        forces = drag.forces(pos, np.ones_like(pos))
        assert np.all(np.isfinite(forces))

    def test_mismatched_shapes_are_rejected(self, plan: BodyPlan) -> None:
        drag = GroundDrag(plan)
        with pytest.raises(ValueError, match="must match"):
            drag.forces(self.straight_body(plan), np.zeros((3, 2)))
