"""``fixed_point`` must find an equilibrium, or refuse.

Model assumptions 5AN.8: the loop's eigenvalues were once taken wherever a
trajectory happened to be after twenty seconds, which was the top of a burst on
a 22-second cycle, and every conclusion drawn from them had to be withdrawn.
``fixed_point`` exists so that cannot happen again, so the two things it must do
are tested on maps whose answer is known: find the equilibrium of a map that has
one, including when the trajectory never reaches it, and raise for a map that
has none.
"""

from __future__ import annotations

import numpy as np
import pytest

from tools.analyse_loop import fixed_point, jacobian


class _AffineLoop:
    """A stand-in with the ``Loop`` interface: ``x -> A x + b``.

    ``sizes`` follows ``Loop``: voltages, activations, muscles, joints. The
    activation block is clipped to [0, 1] by ``fixed_point``, so the equilibrium
    is placed inside it.
    """

    def __init__(self, a: np.ndarray, b: np.ndarray) -> None:
        self.a, self.b = a, b
        self.dt = 0.01
        self.sizes = (2, 2, 1, 1)
        self.x = np.zeros(b.size)

    def joints(self, x: np.ndarray) -> np.ndarray:
        n, _, m, j = self.sizes
        return x[2 * n + m : 2 * n + m + j]

    def step(self) -> None:
        self.x = self.a @ self.x + self.b

    def get(self) -> np.ndarray:
        return self.x.copy()

    def map(self, x: np.ndarray) -> np.ndarray:
        return self.a @ x + self.b


def _rotation(angle: float, radius: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return radius * np.array([[c, -s], [s, c]])


def _system(radius: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A 6-state map whose equilibrium is known, with an oscillating block."""
    a = np.zeros((6, 6))
    a[0:2, 0:2] = _rotation(0.3, radius)  # the voltages oscillate
    a[2:4, 2:4] = 0.5 * np.eye(2)
    a[4, 4] = 0.8
    a[5, 5] = 0.9
    x_star = np.array([-30.0, -20.0, 0.3, 0.6, 0.2, 0.1])
    b = x_star - a @ x_star
    return a, b, x_star


def test_finds_a_stable_equilibrium() -> None:
    a, b, x_star = _system(radius=0.95)
    x, jac = fixed_point(_AffineLoop(a, b), settle_s=0.5, average_s=1.0)
    assert np.allclose(x, x_star, atol=1e-8)
    assert np.allclose(jac, a, atol=1e-6)


def test_finds_an_unstable_equilibrium_the_trajectory_never_reaches() -> None:
    """The committed loop's case: the trajectory runs away from its equilibrium."""
    a, b, x_star = _system(radius=1.02)
    loop = _AffineLoop(a, b)
    x, _ = fixed_point(loop, settle_s=0.5, average_s=1.0)
    assert np.allclose(x, x_star, atol=1e-8)
    # And the trajectory itself is nowhere near it: settling would have lied.
    assert np.abs(loop.get() - x_star).max() > 1.0


def test_start_continues_a_given_equilibrium() -> None:
    a, b, x_star = _system(radius=0.95)
    x, _ = fixed_point(_AffineLoop(a, b), start=x_star + 0.01)
    assert np.allclose(x, x_star, atol=1e-8)


def test_refuses_a_map_with_no_equilibrium() -> None:
    """``x -> x + c`` has no fixed point; returning anything would be a lie."""
    a = np.eye(6)
    b = np.full(6, 1.0e-3)
    with pytest.raises((RuntimeError, np.linalg.LinAlgError)):
        fixed_point(_AffineLoop(a, b), settle_s=0.1, average_s=0.1, max_iter=3)


def test_jacobian_of_an_affine_map_is_its_matrix() -> None:
    a, b, x_star = _system(radius=0.9)
    assert np.allclose(jacobian(_AffineLoop(a, b), x_star), a, atol=1e-6)


def test_a_lesion_removes_the_cells_and_nothing_else() -> None:
    """An ablated cell transmits nothing, chemically or electrically; every other
    connection, and every other cell's tuning, is left as built.

    The runner's --lesion holds a cell's chemical output at zero and leaves its gap
    junctions, through which a killed cell would still pass current. The loop's
    ablation cuts both (model_assumptions 5AO.10).
    """
    import argparse

    from tools.analyse_loop import Loop

    def build(lesion: str) -> Loop:
        return Loop(
            argparse.Namespace(
                unknown_sign="exclude",
                physics_hz=240.0,
                neural_dt_ms=1.0,
                torque_scale=1e-3,
                proprioceptive_mv=400.0,
                proprio_rate=0.7,
                command_mv=20.0,
                settle_s=0.0,
                open_loop=False,
                param=[],
                lesion=lesion,
            )
        )

    intact, lesioned = build("").runtime, build("SMD")
    net = lesioned.runtime.network
    assert lesioned.lesioned == ["SMDDL", "SMDDR", "SMDVL", "SMDVR"]
    idx = [net.index(c) for c in lesioned.lesioned]
    keep = [i for i in range(net.n) if i not in idx]
    assert not net.g_syn[:, idx].any() and not net.g_gap[:, idx].any()
    assert not net.g_gap[idx, :].any() and not net.gap_row_sum[idx].any()
    assert intact.network.g_syn[:, idx].any() and intact.network.g_gap[:, idx].any()
    np.testing.assert_array_equal(
        net.g_syn[np.ix_(keep, keep)], intact.network.g_syn[np.ix_(keep, keep)]
    )
    np.testing.assert_array_equal(
        net.g_gap[np.ix_(keep, keep)], intact.network.g_gap[np.ix_(keep, keep)]
    )
    # The survivors keep the sigmoid thresholds they were built with.
    np.testing.assert_array_equal(
        lesioned.runtime.model.v_threshold_mv, intact.model.v_threshold_mv
    )


def test_a_rebuilt_lesion_builds_the_network_without_the_cells() -> None:
    """The rebuilt ablation leaves the cells out of the network altogether, and
    no proprioceptive population still feeds them (model_assumptions 5AO.10)."""
    import argparse

    from tools.analyse_loop import Loop

    loop = Loop(
        argparse.Namespace(
            unknown_sign="exclude",
            physics_hz=240.0,
            neural_dt_ms=1.0,
            torque_scale=1e-3,
            proprioceptive_mv=400.0,
            proprio_rate=0.7,
            command_mv=20.0,
            settle_s=0.0,
            open_loop=False,
            param=[],
            head_proprio=True,
            head_classes=("SMD", "RMD"),
            lesion="SMD",
            lesion_mode="rebuilt",
        )
    )
    cells = set(loop.runtime.network.cell_ids)
    assert loop.lesioned == ["SMDDL", "SMDDR", "SMDVL", "SMDVR"]
    assert not cells & set(loop.lesioned)
    for pop in (loop.proprio, *loop.extra_proprio):
        assert set(pop.targets) <= cells
        assert len(pop.sensed_segment) == len(pop.targets)
    assert any(c.startswith("RMD") for c in loop.extra_proprio[0].targets)
    loop.step()  # feedback reaches only cells that exist


def test_the_loop_can_run_cells_on_their_measured_biophysics() -> None:
    """``intrinsic`` hands named cells to their conductance model (5AO.12): they
    start at the measured resting potential, their sigmoid follows the network's
    offset convention, and the loop steps with them."""
    import argparse

    from tools.analyse_loop import Loop

    loop = Loop(
        argparse.Namespace(
            unknown_sign="exclude",
            physics_hz=240.0,
            neural_dt_ms=1.0,
            torque_scale=1e-3,
            proprioceptive_mv=400.0,
            proprio_rate=0.7,
            command_mv=20.0,
            settle_s=0.0,
            open_loop=False,
            param=[],
            intrinsic="rmd_nicoletti2019",
        )
    )
    assert loop.hybrid is not None
    assert set(loop.hybrid.cells) == {"RMDDL", "RMDDR", "RMDVL", "RMDVR", "RMDL", "RMDR"}
    rest = float(loop.hybrid.model.resting_state()[0, 0])
    offset = loop.runtime.model.params.v_threshold_offset_mv
    np.testing.assert_allclose(
        loop.runtime.model.v_threshold_mv[loop.hybrid.indices], rest + offset
    )
    np.testing.assert_allclose(loop.runtime.state[0, loop.hybrid.indices], rest)
    for _ in range(24):  # 0.1 s
        loop.step()
    v = loop.runtime.state[0, loop.hybrid.indices]
    assert np.all(np.isfinite(v))
    np.testing.assert_array_equal(v, loop.hybrid.state[0])  # the conductance model owns them
