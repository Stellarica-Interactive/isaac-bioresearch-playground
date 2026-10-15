"""Extrasynaptic links: couplings the connectome cannot show (model_assumptions 5AO.11).

What has to hold: a link adds exactly the conductance asked for, with its sign's
reversal potential; a link onto an existing connection combines with it without
changing the current that connection carries; links naming absent cells are
skipped; and the network is re-solved, as a build with the links would be.
"""

from __future__ import annotations

import numpy as np
import pytest

from worm.neural.config import (
    ExtrasynapticLink,
    add_extrasynaptic,
    build_runtime,
    smd_rme_loop,
)


@pytest.fixture(scope="module")
def runtime_factory():
    def make():
        runtime, _ = build_runtime("cook_2019_herm", unknown_sign="exclude")
        return runtime

    return make


def test_a_link_adds_its_conductance_with_its_sign(runtime_factory) -> None:
    runtime = runtime_factory()
    net = runtime.network
    i, j = net.index("RMEV"), net.index("SMDVL")
    assert net.g_syn[i, j] == 0.0  # not connected in Cook 2019
    added = add_extrasynaptic(runtime, [ExtrasynapticLink("SMDVL", "RMEV", 1.5, "excitatory")])
    assert added == ["SMDVL->RMEV"]
    assert runtime.network.g_syn[i, j] == pytest.approx(1.5)
    assert runtime.network.e_rev[i, j] == pytest.approx(runtime.model.params.e_exc_mv)


def test_a_link_onto_a_connection_combines_with_it_exactly(runtime_factory) -> None:
    """SMDVR already synapses onto RMED. Adding an inhibitory link there must give
    the same current, at every voltage, as the two connections side by side."""
    runtime = runtime_factory()
    i, j = runtime.network.index("RMED"), runtime.network.index("SMDVR")
    g_old, e_old = runtime.network.g_syn[i, j], runtime.network.e_rev[i, j]
    assert g_old > 0.0
    e_inh = runtime.model.params.e_inh_mv
    add_extrasynaptic(runtime, [ExtrasynapticLink("SMDVR", "RMED", 0.5, "inhibitory")])
    g_new, e_new = runtime.network.g_syn[i, j], runtime.network.e_rev[i, j]
    for v in (-70.0, -40.0, 0.0):
        assert g_new * (e_new - v) == pytest.approx(g_old * (e_old - v) + 0.5 * (e_inh - v))


def test_links_to_absent_cells_are_skipped_and_the_network_is_resolved(runtime_factory) -> None:
    runtime = runtime_factory()
    before = runtime.model.v_threshold_mv.copy()
    added = add_extrasynaptic(
        runtime,
        [
            ExtrasynapticLink("NOT_A_CELL", "RMED", 1.0, "excitatory"),
            ExtrasynapticLink("SMDDL", "RMED", 2.0, "excitatory"),
        ],
    )
    assert added == ["SMDDL->RMED"]
    rmed = runtime.network.index("RMED")
    assert runtime.model.v_threshold_mv[rmed] != pytest.approx(before[rmed])
    np.testing.assert_array_equal(runtime.state, runtime.model.initial_state(runtime.network))


def test_link_validation() -> None:
    runtime, _ = build_runtime("cook_2019_herm", unknown_sign="exclude")
    with pytest.raises(ValueError):
        add_extrasynaptic(runtime, [ExtrasynapticLink("SMDDL", "RMED", 1.0, "modulatory")])
    with pytest.raises(ValueError):
        add_extrasynaptic(runtime, [ExtrasynapticLink("SMDDL", "RMED", -1.0, "excitatory")])


def test_the_smd_rme_loop_pairs_the_cells_recorded_together() -> None:
    """SMDD drives RMED and SMDV drives RMEV (Shen et al. 2016); "same" sends each
    RME's inhibition back to its own SMD cells, "all" to all four."""
    same = smd_rme_loop(1.0, 2.0, pairing="same")
    excite = {(link.pre, link.post) for link in same if link.sign == "excitatory"}
    assert excite == {
        ("SMDDL", "RMED"),
        ("SMDDR", "RMED"),
        ("SMDVL", "RMEV"),
        ("SMDVR", "RMEV"),
    }
    inhibit = {(link.pre, link.post) for link in same if link.sign == "inhibitory"}
    assert inhibit == {(post, pre) for pre, post in excite}
    assert all(link.g_ns == 1.0 for link in same if link.sign == "excitatory")
    assert all(link.g_ns == 2.0 for link in same if link.sign == "inhibitory")

    every = smd_rme_loop(1.0, 2.0, pairing="all")
    assert len([link for link in every if link.sign == "inhibitory"]) == 8
    with pytest.raises(ValueError):
        smd_rme_loop(1.0, 1.0, pairing="crossed")
