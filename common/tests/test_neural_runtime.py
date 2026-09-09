"""Neural runtime: correctness against analytic solutions, stability, determinism.

Most of these run on hand-built two- and three-cell networks where the right
answer can be worked out on paper. A simulator that agrees with a closed-form
solution on a trivial network is not proven correct on a real one, but a
simulator that *disagrees* is definitely wrong, and these catch that.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import cast

import numpy as np
import pytest

from common.data.schemas import (
    Cell,
    CellCategory,
    Confidence,
    Connection,
    Connectome,
    Provenance,
    Scope,
    Sign,
    SynapseType,
)
from common.neural.neuron_models import (
    GradedLeakyIntegrator,
    GradedLeakyIntegratorParameters,
    prepare,
    sigmoid,
)
from common.neural.runtime import Integrator, NeuralRuntime
from common.neural.synapses import (
    NetworkMatrices,
    UnknownSignPolicy,
    WeightScaling,
    build_matrices,
)

PROV = Provenance(
    source_id="test",
    kind="connectome",
    citation="Test fixture",
    url="https://example.invalid/",
    license="CC0",
)


def net(cells: list[str], conns: list[Connection]) -> Connectome:
    return Connectome(
        id="test",
        organism="Test",
        sex="hermaphrodite",
        stage="adult",
        scope=Scope.WHOLE_ANIMAL,
        cells=tuple(Cell(id=c, category=CellCategory.NEURON) for c in cells),
        connections=tuple(conns),
        provenance=PROV,
    )


def signed(pre: str, post: str, weight: int, sign: Sign) -> Connection:
    return Connection(
        pre, post, SynapseType.CHEMICAL, weight, sign=sign, sign_confidence=Confidence.PREDICTED
    )


P = GradedLeakyIntegratorParameters()


# ---------------------------------------------------------------------------
# Analytic checks
# ---------------------------------------------------------------------------


class TestIsolatedNeuron:
    """One cell, no connections: an exponential relaxation with a known solution."""

    def build(self, integrator: Integrator, dt: float = 0.5) -> NeuralRuntime:
        return NeuralRuntime.build(
            net(["A", "B"], []),
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            dt_ms=dt,
            integrator=integrator,
        )

    def test_rests_at_the_leak_potential(self) -> None:
        rt = self.build(Integrator.RK4)
        assert rt.voltage("A") == pytest.approx(P.e_leak_mv, abs=1e-9)

    def test_stays_at_rest_without_input(self) -> None:
        rt = self.build(Integrator.RK4)
        rt.run(1000)
        assert rt.voltage("A") == pytest.approx(P.e_leak_mv, abs=1e-6)

    @pytest.mark.parametrize(
        "integrator", [Integrator.EULER, Integrator.RK4, Integrator.EXPONENTIAL]
    )
    def test_steady_state_under_injected_current(self, integrator: Integrator) -> None:
        """V_inf = E_leak + I / G_leak, exactly."""
        rt = self.build(integrator, dt=0.5)
        rt.inject("A", 0.05)  # pA
        rt.run(3000)  # many time constants (tau = 100 ms)
        expected = P.e_leak_mv + 0.05 / P.g_leak_ns
        assert rt.voltage("A") == pytest.approx(expected, rel=1e-4)

    def test_relaxation_follows_the_exponential_solution(self) -> None:
        """V(t) = V_inf + (V0 - V_inf) exp(-t/tau), tau = C / G_leak = 100 ms."""
        v0 = P.e_leak_mv
        v_inf = P.e_leak_mv + 0.05 / P.g_leak_ns
        tau = P.c_m_pf / P.g_leak_ns
        rt = self.build(Integrator.RK4, dt=0.25)
        rt.inject("A", 0.05)
        for t in (10.0, 50.0, 100.0, 250.0):
            rt.run(t - rt.t_ms)
            expected = v_inf + (v0 - v_inf) * math.exp(-t / tau)
            assert rt.voltage("A") == pytest.approx(expected, rel=1e-5), t

    def test_membrane_time_constant_is_100_ms(self) -> None:
        rt = self.build(Integrator.RK4)
        assert rt.fastest_time_constant_ms() == pytest.approx(100.0, rel=1e-9)


class TestGapJunctions:
    def test_coupling_pulls_two_cells_together(self) -> None:
        """A gap junction drags an undriven cell toward its driven neighbour.

        Note they do *not* fully equalize: the junction has finite conductance, so a
        steady current from A through the junction and out through B's leak requires
        a standing voltage difference. Stronger coupling shrinks that difference,
        which is what this asserts.
        """
        base = NeuralRuntime.build(
            net(["A", "B"], []), unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=0.5
        )
        base.inject("A", 0.2)
        base.run(2000)
        undriven_rest = base.voltage("B")

        gaps = {}
        for weight in (1, 10, 100):
            c = net(["A", "B"], [Connection("A", "B", SynapseType.ELECTRICAL, weight)])
            rt = NeuralRuntime.build(
                c,
                unknown_sign=UnknownSignPolicy.EXCLUDE,
                dt_ms=0.5,
                integrator=Integrator.EXPONENTIAL,
            )
            rt.inject("A", 0.2)
            rt.run(2000)
            gaps[weight] = abs(rt.voltage("A") - rt.voltage("B"))
            assert rt.voltage("B") > undriven_rest, weight

        assert gaps[1] > gaps[10] > gaps[100]
        assert gaps[100] < 0.05

    def test_gap_current_is_antisymmetric(self) -> None:
        """Whatever leaves one cell enters the other; the pair conserves current."""
        c = net(["A", "B"], [Connection("A", "B", SynapseType.ELECTRICAL, 3)])
        rt = NeuralRuntime.build(c, unknown_sign=UnknownSignPolicy.EXCLUDE)
        rt.inject("A", 0.5)
        rt.run(50)
        model = cast(GradedLeakyIntegrator, rt.model)
        terms = model.currents(rt.state[0], rt.state[1], rt.i_ext_pa, rt.network)
        assert terms["gap"].sum() == pytest.approx(0.0, abs=1e-12)

    def test_self_gap_junction_is_dropped(self) -> None:
        """Electrically inert in a single-compartment model: driving force is zero."""
        c = net(["A"], [Connection("A", "A", SynapseType.ELECTRICAL, 5)])
        m = build_matrices(
            c,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            g_syn_ps=100,
            g_gap_ps=100,
            e_exc_mv=0,
            e_inh_mv=-45,
            e_leak_mv=-35,
        )
        assert m.g_gap[0, 0] == 0.0

    def test_gap_matrix_is_symmetric(self) -> None:
        c = net(["A", "B", "C"], [Connection("A", "C", SynapseType.ELECTRICAL, 2)])
        m = build_matrices(
            c,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            g_syn_ps=100,
            g_gap_ps=100,
            e_exc_mv=0,
            e_inh_mv=-45,
            e_leak_mv=-35,
        )
        assert np.array_equal(m.g_gap, m.g_gap.T)


class TestChemicalSynapses:
    def two_cell(self, sign: Sign) -> NeuralRuntime:
        c = net(["A", "B"], [signed("A", "B", 10, sign)])
        return NeuralRuntime.build(
            c, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=0.2, integrator=Integrator.RK4
        )

    def test_excitatory_input_depolarizes_the_target(self) -> None:
        rt = self.two_cell(Sign.EXCITATORY)
        before = rt.voltage("B")
        rt.inject("A", 2.0)
        rt.run(2000)
        assert rt.voltage("B") > before

    def test_inhibitory_input_hyperpolarizes_the_target(self) -> None:
        rt = self.two_cell(Sign.INHIBITORY)
        before = rt.voltage("B")
        rt.inject("A", 2.0)
        rt.run(2000)
        assert rt.voltage("B") < before

    def test_signal_flows_only_in_the_synapse_direction(self) -> None:
        """Chemical synapses are directed; driving the target must not move the source."""
        rt = self.two_cell(Sign.EXCITATORY)
        before_a = rt.voltage("A")
        rt.inject("B", 2.0)
        rt.run(500)
        assert rt.voltage("A") == pytest.approx(before_a, abs=1e-9)

    def test_reversal_potential_is_per_connection_not_per_neuron(self) -> None:
        """One neuron can excite one target and inhibit another.

        The receptor deciding a synapse's sign is postsynaptic, so this must be
        representable. Most published implementations collapse it to one value per
        presynaptic cell and cannot express it.
        """
        c = net(
            ["A", "B", "C"],
            [signed("A", "B", 5, Sign.EXCITATORY), signed("A", "C", 5, Sign.INHIBITORY)],
        )
        rt = NeuralRuntime.build(
            c, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=0.2, integrator=Integrator.RK4
        )
        b0, c0 = rt.voltage("B"), rt.voltage("C")
        rt.inject("A", 2.0)
        rt.run(2000)
        assert rt.voltage("B") > b0
        assert rt.voltage("C") < c0


class TestEquilibrium:
    def test_initial_state_is_a_true_equilibrium(self) -> None:
        """The solved resting state must have (almost) zero voltage derivative."""
        c = net(
            ["A", "B", "C"],
            [
                signed("A", "B", 4, Sign.EXCITATORY),
                signed("B", "C", 7, Sign.INHIBITORY),
                Connection("A", "C", SynapseType.ELECTRICAL, 2),
            ],
        )
        rt = NeuralRuntime.build(c, unknown_sign=UnknownSignPolicy.EXCLUDE)
        d = rt.model.derivatives(rt.state, rt.i_ext_pa, rt.network)
        assert np.max(np.abs(d[0])) < 1e-9

    def test_sigmoids_start_at_their_midpoint(self) -> None:
        c = net(["A", "B"], [signed("A", "B", 3, Sign.EXCITATORY)])
        rt = NeuralRuntime.build(c, unknown_sign=UnknownSignPolicy.EXCLUDE)
        model = cast(GradedLeakyIntegrator, rt.model)
        assert model.v_threshold_mv is not None
        phi = sigmoid(rt.state[0], model.v_threshold_mv, model.params.beta_per_mv)
        assert np.allclose(phi, 0.5)

    def test_resting_activation_matches_the_analytic_value(self) -> None:
        s_star = P.s_at_half_activation
        assert s_star == pytest.approx(
            0.5 * P.a_r_per_ms / (0.5 * P.a_r_per_ms + P.a_d_per_ms)
        )
        assert s_star == pytest.approx(1.0 / 11.0)


class TestSigmoid:
    def test_midpoint_is_a_half(self) -> None:
        assert sigmoid(np.array([5.0]), np.array([5.0]), 0.125)[0] == pytest.approx(0.5)

    def test_saturates_without_overflow(self) -> None:
        """The naive form overflows for large negative arguments; this must not."""
        v = np.array([-1e6, 1e6])
        out = sigmoid(v, np.zeros(2), 0.125)
        assert np.all(np.isfinite(out))
        assert out[0] == pytest.approx(0.0)
        assert out[1] == pytest.approx(1.0)

    def test_monotonic(self) -> None:
        v = np.linspace(-100, 100, 500)
        out = sigmoid(v, np.zeros(500), 0.125)
        assert np.all(np.diff(out) >= 0)


# ---------------------------------------------------------------------------
# Numerical behaviour
# ---------------------------------------------------------------------------


class TestIntegrators:
    def network(self) -> Connectome:
        return net(
            ["A", "B", "C", "D"],
            [
                signed("A", "B", 8, Sign.EXCITATORY),
                signed("B", "C", 5, Sign.INHIBITORY),
                signed("C", "D", 6, Sign.EXCITATORY),
                Connection("A", "D", SynapseType.ELECTRICAL, 3),
            ],
        )

    def run(self, integrator: Integrator, dt: float, duration: float = 40.0) -> np.ndarray:
        rt = NeuralRuntime.build(
            self.network(),
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            dt_ms=dt,
            integrator=integrator,
        )
        rt.inject("A", 1.5)
        rt.run(duration)
        return rt.state.copy()

    def test_all_three_agree_at_a_small_timestep(self) -> None:
        reference = self.run(Integrator.RK4, 0.005)
        for integrator in (Integrator.EULER, Integrator.EXPONENTIAL):
            got = self.run(integrator, 0.005)
            assert np.max(np.abs(got[0] - reference[0])) < 0.01, integrator

    def test_rk4_converges_faster_than_euler(self) -> None:
        reference = self.run(Integrator.RK4, 0.005)
        err_euler = np.max(np.abs(self.run(Integrator.EULER, 0.5)[0] - reference[0]))
        err_rk4 = np.max(np.abs(self.run(Integrator.RK4, 0.5)[0] - reference[0]))
        assert err_rk4 < err_euler

    def test_exponential_survives_a_timestep_that_breaks_euler(self) -> None:
        """The point of having it: stability where explicit methods diverge.

        The stiff network here has a membrane time constant of about 0.025 ms, so
        dt = 5 ms is two hundred times over the explicit stability limit.
        """
        stiff = net(["A", "B"], [Connection("A", "B", SynapseType.ELECTRICAL, 400)])

        euler = NeuralRuntime.build(
            stiff,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            dt_ms=5.0,
            integrator=Integrator.EULER,
        )
        euler.inject("A", 1.0)
        with np.errstate(over="ignore", invalid="ignore"):
            try:
                euler.run(500)
                diverged = float(np.max(np.abs(euler.state[0]))) > 1e6
            except FloatingPointError:
                diverged = True
        assert diverged, "forward Euler should be unusable at 200x its stability limit"

        exp = NeuralRuntime.build(
            stiff,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            dt_ms=5.0,
            integrator=Integrator.EXPONENTIAL,
        )
        exp.inject("A", 1.0)
        exp.run(500)
        assert np.all(np.isfinite(exp.state))
        assert np.max(np.abs(exp.state[0])) < 200.0


class TestStability:
    def test_report_flags_an_oversized_timestep(self) -> None:
        stiff = net(["A", "B"], [Connection("A", "B", SynapseType.ELECTRICAL, 400)])
        rt = NeuralRuntime.build(
            stiff, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=5.0, integrator=Integrator.EULER
        )
        report = rt.stability_report()
        assert report.stable is False
        assert report.recommended_dt_ms < report.dt_ms
        assert report.dt_over_tau > 100

    def test_report_accepts_a_reasonable_timestep(self) -> None:
        rt = NeuralRuntime.build(net(["A"], []), unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=1.0)
        assert rt.stability_report().stable is True

    def test_exponential_is_reported_as_unconditionally_stable(self) -> None:
        stiff = net(["A", "B"], [Connection("A", "B", SynapseType.ELECTRICAL, 400)])
        rt = NeuralRuntime.build(
            stiff,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            dt_ms=1000.0,
            integrator=Integrator.EXPONENTIAL,
        )
        assert rt.stability_report().stable is True

    def test_divergence_raises_with_an_actionable_message(self) -> None:
        """When it blows up, the error must say why and what to do about it."""
        stiff = net(["A", "B"], [Connection("A", "B", SynapseType.ELECTRICAL, 400)])
        rt = NeuralRuntime.build(
            stiff, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=50.0, integrator=Integrator.EULER
        )
        rt.inject("A", 5.0)
        # The overflow is the point of the test; do not let numpy warn about it.
        with (
            np.errstate(over="ignore", invalid="ignore"),
            pytest.raises(FloatingPointError, match="EXPONENTIAL"),
        ):
            rt.run(5000)


class TestDeterminism:
    def build(self) -> NeuralRuntime:
        c = net(
            ["A", "B", "C"],
            [signed("A", "B", 4, Sign.EXCITATORY), signed("B", "C", 3, Sign.INHIBITORY)],
        )
        rt = NeuralRuntime.build(c, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=0.1)
        rt.inject("A", 1.0)
        return rt

    def test_two_runs_are_bit_identical(self) -> None:
        a, b = self.build(), self.build()
        a.run(100)
        b.run(100)
        assert np.array_equal(a.state, b.state)

    def test_step_size_composition_is_exact(self) -> None:
        """1000 steps of dt equals 500 then 500, to the bit."""
        a, b = self.build(), self.build()
        a.run(100)
        b.run(50)
        b.run(50)
        assert np.array_equal(a.state, b.state)
        assert a.step_count == b.step_count


class TestCheckpointing:
    def build(self) -> NeuralRuntime:
        c = net(
            ["A", "B", "C"],
            [
                signed("A", "B", 4, Sign.EXCITATORY),
                signed("B", "C", 3, Sign.INHIBITORY),
                Connection("A", "C", SynapseType.ELECTRICAL, 2),
            ],
        )
        return NeuralRuntime.build(c, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=0.1)

    def test_restore_continues_identically(self, tmp_path: Path) -> None:
        """Save, stop, restore, continue: numerically identical to never stopping."""
        a = self.build()
        a.inject("A", 1.0)
        a.run(50)
        a.save_checkpoint(tmp_path / "ck.npz")
        a.run(50)

        b = self.build()
        b.load_checkpoint(tmp_path / "ck.npz")
        b.run(50)

        assert np.array_equal(a.state, b.state)
        assert a.t_ms == pytest.approx(b.t_ms)
        assert a.step_count == b.step_count

    def test_external_input_is_part_of_the_state(self, tmp_path: Path) -> None:
        a = self.build()
        a.inject("A", 3.25)
        a.save_checkpoint(tmp_path / "ck.npz")
        b = self.build()
        b.load_checkpoint(tmp_path / "ck.npz")
        assert b.i_ext_pa[b.network.index("A")] == pytest.approx(3.25)

    def test_restoring_into_a_different_network_is_refused(self, tmp_path: Path) -> None:
        a = self.build()
        a.run(10)
        a.save_checkpoint(tmp_path / "ck.npz")

        lesioned = net(
            ["A", "B", "C"],
            [signed("A", "B", 4, Sign.EXCITATORY)],  # the B->C synapse removed
        )
        b = NeuralRuntime.build(lesioned, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=0.1)
        with pytest.raises(ValueError, match="different network"):
            b.load_checkpoint(tmp_path / "ck.npz")

    def test_a_lesion_experiment_can_opt_in(self, tmp_path: Path) -> None:
        """Restoring a neural state into a deliberately altered network is allowed,
        but only by asking for it."""
        a = self.build()
        a.run(10)
        a.save_checkpoint(tmp_path / "ck.npz")
        lesioned = net(["A", "B", "C"], [signed("A", "B", 4, Sign.EXCITATORY)])
        b = NeuralRuntime.build(lesioned, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=0.1)
        b.load_checkpoint(tmp_path / "ck.npz", allow_different_network=True)
        assert b.step_count == a.step_count

    def test_cell_set_mismatch_is_refused(self, tmp_path: Path) -> None:
        a = self.build()
        a.save_checkpoint(tmp_path / "ck.npz")
        other = NeuralRuntime.build(
            net(["A", "B"], []), unknown_sign=UnknownSignPolicy.EXCLUDE
        )
        with pytest.raises(ValueError, match="different set of cells"):
            other.load_checkpoint(tmp_path / "ck.npz")


# ---------------------------------------------------------------------------
# Sign handling
# ---------------------------------------------------------------------------


class TestUnknownSignPolicy:
    def unsigned(self) -> Connectome:
        return net(["A", "B"], [Connection("A", "B", SynapseType.CHEMICAL, 5)])

    def matrices(self, policy: UnknownSignPolicy) -> NetworkMatrices:
        return build_matrices(
            self.unsigned(),
            unknown_sign=policy,
            g_syn_ps=100,
            g_gap_ps=100,
            e_exc_mv=0.0,
            e_inh_mv=-45.0,
            e_leak_mv=-35.0,
        )

    def test_exclude_drops_the_connection(self) -> None:
        m = self.matrices(UnknownSignPolicy.EXCLUDE)
        assert m.g_syn.sum() == 0.0
        assert m.census.excluded == 1

    def test_excitatory_policy_assigns_the_excitatory_reversal(self) -> None:
        m = self.matrices(UnknownSignPolicy.EXCITATORY)
        assert m.e_rev[m.index("B"), m.index("A")] == 0.0
        assert m.census.from_policy == 1

    def test_neutral_policy_shunts_at_the_leak_potential(self) -> None:
        m = self.matrices(UnknownSignPolicy.NEUTRAL)
        assert m.e_rev[m.index("B"), m.index("A")] == -35.0

    def test_policy_is_required(self) -> None:
        """No default: half of real connections are unsigned, so this must be chosen."""
        with pytest.raises(TypeError):
            NeuralRuntime.build(self.unsigned())  # type: ignore[call-arg]

    def test_census_counts_add_up(self) -> None:
        c = net(
            ["A", "B", "C"],
            [
                signed("A", "B", 1, Sign.EXCITATORY),
                signed("A", "C", 1, Sign.INHIBITORY),
                signed("B", "C", 1, Sign.MIXED),
                Connection("C", "A", SynapseType.CHEMICAL, 1),
            ],
        )
        m = build_matrices(
            c,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            g_syn_ps=100,
            g_gap_ps=100,
            e_exc_mv=0,
            e_inh_mv=-45,
            e_leak_mv=-35,
        )
        assert (m.census.excitatory, m.census.inhibitory, m.census.mixed) == (1, 1, 1)
        assert m.census.excluded == 1
        assert m.census.total == 4
        assert m.census.measured_fraction == pytest.approx(0.75)

    def test_mixed_sign_shunts_rather_than_picking_a_side(self) -> None:
        """The source says the net effect is undetermined; we do not resolve it."""
        c = net(["A", "B"], [signed("A", "B", 5, Sign.MIXED)])
        m = build_matrices(
            c,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            g_syn_ps=100,
            g_gap_ps=100,
            e_exc_mv=0,
            e_inh_mv=-45,
            e_leak_mv=-35,
        )
        assert m.e_rev[m.index("B"), m.index("A")] == -35.0


class TestWeightScaling:
    def c(self) -> Connectome:
        return net(["A", "B"], [signed("A", "B", 9, Sign.EXCITATORY)])

    def scaled(self, scaling: WeightScaling) -> float:
        m = build_matrices(
            self.c(),
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            g_syn_ps=100,
            g_gap_ps=100,
            e_exc_mv=0,
            e_inh_mv=-45,
            e_leak_mv=-35,
            weight_scaling=scaling,
        )
        return float(m.g_syn.max())

    def test_linear(self) -> None:
        assert self.scaled(WeightScaling.LINEAR) == pytest.approx(9 * 0.1)

    def test_sqrt(self) -> None:
        assert self.scaled(WeightScaling.SQRT) == pytest.approx(3 * 0.1)

    def test_binary_discards_the_count(self) -> None:
        assert self.scaled(WeightScaling.BINARY) == pytest.approx(0.1)


class TestRecording:
    def test_records_at_the_requested_interval(self) -> None:
        rt = NeuralRuntime.build(
            net(["A", "B"], []), unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=1.0
        )
        rec = rt.record(["A"], every=10)
        rt.run(100)
        times, states = rec.as_arrays()
        assert len(times) == 10
        assert states.shape == (10, 2, 1)


class TestModelProtocol:
    def test_graded_integrator_satisfies_the_protocol(self) -> None:
        from common.neural.neuron_models import NeuronModel

        m = GradedLeakyIntegrator()
        assert isinstance(m, NeuronModel)

    def test_prepare_caches_thresholds(self) -> None:
        c = net(["A", "B"], [signed("A", "B", 2, Sign.EXCITATORY)])
        m = build_matrices(
            c,
            unknown_sign=UnknownSignPolicy.EXCLUDE,
            g_syn_ps=100,
            g_gap_ps=100,
            e_exc_mv=0,
            e_inh_mv=-45,
            e_leak_mv=-35,
        )
        model = GradedLeakyIntegrator()
        assert model.v_threshold_mv is None
        assert prepare(model, m).v_threshold_mv is not None
