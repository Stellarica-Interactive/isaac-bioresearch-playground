"""Scaling gap junctions by class must be the same model as scaling ``g_gap_ps``.

model_assumptions 5AN.20 found the forward crawl only with every gap junction
twenty times weaker. :func:`rescale_gap_junctions` exists to ask which kind of
junction that actually needed, and the question is only fair if scaling all
three kinds by one factor is *exactly* the model that changing ``g_gap_ps`` by
that factor builds -- conductances, re-solved sigmoid midpoints and starting
state alike. Otherwise a per-class result would differ from the global one for
a reason that is not biology.
"""

from __future__ import annotations

import numpy as np
import pytest

from common.data.schemas import CellCategory
from common.neural.synapses import UnknownSignPolicy
from worm.importers.naming import body_wall_muscle_ids
from worm.loader import load
from worm.neural.config import RUNTIME_OVERLAYS, build_runtime, rescale_gap_junctions


@pytest.fixture(scope="module")
def connectome():
    c, _ = load("cook_2019_herm", annotations=RUNTIME_OVERLAYS)
    return c


def _build(connectome, overrides=None):
    muscles = body_wall_muscle_ids()
    cells = tuple(
        c.id for c in connectome.cells if c.category is CellCategory.NEURON or c.id in muscles
    )
    runtime, _ = build_runtime(
        "cook_2019_herm",
        unknown_sign=UnknownSignPolicy("exclude"),
        cells=cells,
        connectome=connectome,
        parameter_overrides=overrides,
    )
    return runtime


def test_scaling_every_class_is_scaling_g_gap_ps(connectome) -> None:
    rebuilt = _build(connectome, {"g_gap_ps": 5.0})
    scaled = _build(connectome)
    rescale_gap_junctions(
        scaled, connectome, {"neuron": 0.05, "muscle": 0.05, "neuron-muscle": 0.05}
    )
    assert np.allclose(scaled.network.g_gap, rebuilt.network.g_gap, rtol=1e-12, atol=0.0)
    assert np.allclose(scaled.network.gap_row_sum, rebuilt.network.gap_row_sum, rtol=1e-12)
    assert np.allclose(scaled.model.v_threshold_mv, rebuilt.model.v_threshold_mv, atol=1e-9)
    assert np.allclose(scaled.state, rebuilt.state, atol=1e-9)


def test_one_class_touches_only_its_own_junctions(connectome) -> None:
    base = _build(connectome)
    scaled = _build(connectome)
    rescale_gap_junctions(scaled, connectome, {"muscle": 0.1})
    category = {c.id: c.category for c in connectome.cells}
    muscle = np.array([category[c] is CellCategory.MUSCLE for c in base.network.cell_ids])
    both = muscle[:, None] & muscle[None, :]
    present = base.network.g_gap > 0
    assert np.any(both & present), "the body-wall muscles have gap junctions between them"
    assert np.allclose(scaled.network.g_gap[both], 0.1 * base.network.g_gap[both])
    assert np.array_equal(scaled.network.g_gap[~both], base.network.g_gap[~both])
    # Symmetric before, symmetric after: a gap junction is one pore, not two.
    assert np.array_equal(scaled.network.g_gap, scaled.network.g_gap.T)


def test_factor_one_changes_nothing(connectome) -> None:
    base = _build(connectome)
    scaled = _build(connectome)
    rescale_gap_junctions(scaled, connectome, {"neuron": 1.0, "muscle": 1.0, "neuron-muscle": 1.0})
    assert np.array_equal(scaled.network.g_gap, base.network.g_gap)
    assert np.allclose(scaled.model.v_threshold_mv, base.model.v_threshold_mv, atol=1e-12)


def test_an_unknown_class_is_refused(connectome) -> None:
    with pytest.raises(ValueError, match="unknown gap-junction classes"):
        rescale_gap_junctions(_build(connectome), connectome, {"muscles": 0.5})
