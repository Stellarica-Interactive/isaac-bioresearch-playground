"""Display geometry: the skin, the ground texture, the connectome drawing.

None of it touches the simulation, so these check only that the shapes are right
and the outputs well-formed and deterministic.
"""

from __future__ import annotations

import struct
import zlib

import numpy as np
import pytest

from worm.body.geometry import BodyPlan
from worm.body.visuals import (
    ACTIVITY_STOPS,
    Beside,
    ConnectomeDrawing,
    ConnectomePanel,
    activity_colours,
    agar_texture,
    body_radius,
    connectome_layout,
    encode_png,
    linear,
    place,
    ring_colours,
    skin_geometry,
    skin_rings,
    skin_topology,
    sphere,
    spine,
)

PLAN = BodyPlan()


def straight_centres(plan: BodyPlan = PLAN) -> np.ndarray:
    return np.array([[s.centre_x_m, 0.0, plan.max_radius_m] for s in plan.segments()])


def bent_centres(plan: BodyPlan = PLAN) -> np.ndarray:
    x = np.array([s.centre_x_m for s in plan.segments()])
    y = 0.1 * plan.total_length_m * np.sin(2 * np.pi * x / (0.7 * plan.total_length_m))
    return np.column_stack([x, y, np.zeros_like(x)])


def test_the_spine_passes_through_every_segment_centre() -> None:
    centres = bent_centres()
    pts, tangent, _, _ = spine(centres, PLAN, per_segment=4)
    # Four samples a span, and a span of half a segment before the first centre.
    np.testing.assert_allclose(pts[4:-4:4], centres, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(tangent, axis=1), 1.0)


def test_the_skin_reaches_both_tips_and_follows_the_profile() -> None:
    centres = straight_centres()
    pts, _, _, fraction = spine(centres, PLAN)
    assert fraction[0] == 0.0 and fraction[-1] == 1.0
    # Half a segment beyond the first and last centres: the whole body length.
    assert np.isclose(pts[-1, 0] - pts[0, 0], PLAN.total_length_m, rtol=1e-6)
    n_around = 12
    points, normals = skin_geometry(centres, PLAN, n_around=n_around)
    n_rings = skin_rings(PLAN.n_segments)
    assert points.shape == normals.shape == (n_rings * n_around + 2, 3)
    rings = points[: n_rings * n_around].reshape(n_rings, n_around, 3)
    radii = np.linalg.norm(rings - pts[:, None, :], axis=2)
    np.testing.assert_allclose(radii, body_radius(fraction, PLAN)[:, None] * np.ones(n_around))


def test_skin_normals_are_unit_and_point_outwards() -> None:
    centres = bent_centres()
    points, normals = skin_geometry(centres, PLAN, n_around=16)
    pts = spine(centres, PLAN)[0]
    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1.0)
    rings = points[:-2].reshape(len(pts), 16, 3) - pts[:, None, :]
    assert np.all(np.sum(rings * normals[:-2].reshape(len(pts), 16, 3), axis=2) > 0.0)


def test_the_skin_topology_indexes_every_vertex_and_nothing_else() -> None:
    n_rings, n_around = skin_rings(PLAN.n_segments), 16
    counts, indices = skin_topology(n_rings, n_around)
    assert counts.sum() == len(indices)
    assert indices.min() == 0 and indices.max() == n_rings * n_around + 1
    assert set(np.unique(indices)) == set(range(n_rings * n_around + 2))


def test_the_neuron_sphere_is_closed_unit_and_faces_outwards() -> None:
    vertices, counts, indices = sphere(n_lat=4, n_lon=8)
    np.testing.assert_allclose(np.linalg.norm(vertices, axis=1), 1.0)
    assert counts.sum() == len(indices)
    assert set(np.unique(indices)) == set(range(len(vertices)))
    start = np.r_[0, np.cumsum(counts)[:-1]]
    for s, n in zip(start, counts, strict=True):
        a, b, c = vertices[indices[s : s + 3]]
        centre = vertices[indices[s : s + n]].mean(axis=0)
        assert np.dot(np.cross(b - a, c - a), centre) > 0.0  # wound outwards
    # Closed: every edge is shared by exactly two faces.
    edges: dict[tuple[int, int], int] = {}
    for s, n in zip(start, counts, strict=True):
        face = indices[s : s + n]
        for i in range(n):
            key = tuple(sorted((int(face[i]), int(face[(i + 1) % n]))))
            edges[key] = edges.get(key, 0) + 1
    assert set(edges.values()) == {2}


def test_ring_colours_cover_every_vertex() -> None:
    colours = ring_colours(np.ones((PLAN.n_segments, 3)) * 0.5, 9, 6)
    assert colours.shape == (9 * 6 + 2, 3)


def test_the_agar_texture_is_deterministic_and_mid_toned() -> None:
    a, b = agar_texture(128, seed=3), agar_texture(128, seed=3)
    assert a.shape == (128, 128, 3) and a.dtype == np.uint8
    np.testing.assert_array_equal(a, b)
    assert 110 < a.mean() < 200


def test_png_round_trips() -> None:
    rgb = agar_texture(32, seed=1)
    png = encode_png(rgb)
    assert png[:8] == bytes([0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A])  # PNG signature
    width, height = struct.unpack(">II", png[16:24])
    assert (width, height) == (32, 32)
    start = png.index(b"IDAT") + 4
    length = struct.unpack(">I", png[start - 8 : start - 4])[0]
    raw = zlib.decompress(png[start : start + length])
    rows = np.frombuffer(raw, dtype=np.uint8).reshape(32, 1 + 32 * 3)
    np.testing.assert_array_equal(rows[:, 1:].reshape(32, 32, 3), rgb)


def test_the_layout_puts_motor_neurons_where_their_synapses_are() -> None:
    layout = connectome_layout(
        ["VB5", "AVAL", "PLML", "ALML"],
        {"VB5": "VB", "AVAL": "AVA", "PLML": "PLM", "ALML": "ALM"},
        {"VB5": 11.0},
        PLAN.n_segments,
    )
    assert np.isclose(layout["VB5"][0], 11.5 / PLAN.n_segments)
    assert layout["AVAL"][0] < 0.2  # head
    assert layout["PLML"][0] > 0.85  # tail
    assert 0.35 < layout["ALML"][0] < 0.55  # mid-body
    for s, across, up in layout.values():
        assert 0.0 <= s <= 1.0 and -1.0 <= across <= 1.0 and -1.0 <= up <= 1.0


def test_cells_are_placed_in_the_body_over_it_or_lifted_above_it() -> None:
    centres = straight_centres()
    s, across, up = np.array([0.5, 0.1]), np.array([1.0, -0.5]), np.array([0.0, 0.5])
    inside = place(s, across, up, centres, PLAN, scale=0.8)
    assert np.all(np.abs(inside[:, 1]) <= 0.8 * body_radius(s, PLAN) + 1e-12)
    lifted = place(s, across, up, centres, PLAN, scale=0.8, lift_m=0.02)
    np.testing.assert_allclose(lifted[:, 2] - inside[:, 2], 0.02)
    # Over the top of the skin, so a camera above sees them within the outline.
    surface = place(s, across, up, centres, PLAN, on_surface=True)
    assert np.all(surface[:, 2] > centres[0, 2] + body_radius(s, PLAN))


def test_synapse_paths_run_along_the_body_from_cell_to_cell() -> None:
    names = ("AVAL", "VB5")
    layout = connectome_layout(names, {"VB5": "VB", "AVAL": "AVA"}, {"VB5": 11.0}, PLAN.n_segments)
    drawing = ConnectomeDrawing(
        names=names, layout=layout, edges=(("AVAL", "VB5", 1.0),), indices=np.arange(2)
    )
    centres = straight_centres()
    cells, paths = drawing.positions(centres, PLAN)
    assert cells.shape == (2, 3) and paths.shape == (drawing.points_per_edge, 3)
    s, across, up = (np.array([layout[n][k] for n in names]) for k in range(3))
    np.testing.assert_allclose(cells, place(s, across, up, centres, PLAN))
    np.testing.assert_allclose(paths[[0, -1]], cells, atol=1e-12)
    # Head to mid-body without doubling back.
    steps = np.diff(paths[:, 0]) * np.sign(cells[1, 0] - cells[0, 0])
    assert np.all(steps >= 0.0)


def test_linear_undoes_the_srgb_encoding() -> None:
    np.testing.assert_allclose(linear(np.array([0.0, 1.0])), [0.0, 1.0])
    assert linear(np.array([0.5]))[0] == pytest.approx(0.2140, abs=1e-4)
    assert np.all(np.diff(linear(np.linspace(0.0, 1.0, 50))) > 0.0)


def test_activity_colours_run_from_rest_to_active_without_greying() -> None:
    c = activity_colours(np.array([0.0, 1.0, -1.0, 2.0]))
    np.testing.assert_allclose(c[0], ACTIVITY_STOPS[0])
    np.testing.assert_allclose(c[1], ACTIVITY_STOPS[-1])
    np.testing.assert_allclose(c[[2, 3]], c[[0, 1]])  # clipped
    # Saturated all the way: a straight blue-to-gold blend greys out half-way.
    ramp = activity_colours(np.linspace(0.0, 1.0, 101))
    assert np.all(ramp.max(axis=1) - ramp.min(axis=1) > 0.25)


def test_activity_level_runs_from_rest_to_the_ceiling() -> None:
    drawing = ConnectomeDrawing(
        names=("A", "B", "C"),
        layout={},
        edges=(),
        indices=np.array([2, 0, 1]),
        rest=0.1,
        ceiling=0.5,
    )
    # Network order: index 0 is B, 1 is C, 2 is A.
    np.testing.assert_allclose(drawing.level(np.array([0.5, 0.3, 0.0])), [0.0, 1.0, 0.5])


def test_the_drawing_beside_the_worm_starts_above_it_and_turns_with_it() -> None:
    beside = Beside(distance_m=0.03, tau_s=1.0)
    centres = straight_centres()
    first = beside.offset(centres, 0.0)
    assert first[1] == pytest.approx(0.03) and abs(first[0]) < 1e-12  # the top of the view

    def turned(angle: float) -> np.ndarray:
        c, s = np.cos(angle), np.sin(angle)
        xy = centres[:, :2] @ np.array([[c, s], [-s, c]])
        return np.column_stack([xy, centres[:, 2]])

    previous = first
    for k in range(1, 200):  # half a turn over 20 s
        offset = beside.offset(turned(np.radians(0.9 * k)), 0.1 * k)
        assert np.linalg.norm(offset[:2]) == pytest.approx(0.03)
        assert np.linalg.norm(offset - previous) < 0.003  # never jumps across the body
        previous = offset
    # The same side of the body, which is now the bottom of the view.
    assert previous[1] < -0.025


def test_the_panel_lights_the_active_cell_and_its_synapse() -> None:
    names = ["VB5", "AVAL"]
    layout = connectome_layout(names, {"VB5": "VB", "AVAL": "AVA"}, {"VB5": 11.0}, PLAN.n_segments)
    panel = ConnectomePanel.build(layout, names, [("AVAL", "VB5", 1.0)])
    quiet, busy = panel.draw(np.array([0.0, 0.0])), panel.draw(np.array([0.0, 1.0]))
    assert busy.shape == (200, 640, 4) and busy.dtype == np.uint8
    np.testing.assert_array_equal(busy, panel.draw(np.array([0.0, 1.0])))
    # AVAL's dot turns gold; VB5's stays the resting blue.
    gold = np.round(np.array(ACTIVITY_STOPS[-1]) * 255)
    rest = np.round(np.array(ACTIVITY_STOPS[0]) * 255)
    aval = panel.dot_pixel[panel.dot_cell == 1]
    vb5 = panel.dot_pixel[panel.dot_cell == 0]
    np.testing.assert_allclose(busy.reshape(-1, 4)[aval, :3], np.tile(gold, (len(aval), 1)), atol=1)
    np.testing.assert_allclose(busy.reshape(-1, 4)[vb5, :3], np.tile(rest, (len(vb5), 1)), atol=1)
    # The synapse out of it brightens; the background stays as it was.
    line = np.setdiff1d(panel.line_pixel, panel.dot_pixel)
    assert (
        busy.reshape(-1, 4)[line, :3].astype(int).sum()
        > quiet.reshape(-1, 4)[line, :3].astype(int).sum()
    )
    np.testing.assert_array_equal(busy[0, 0], quiet[0, 0])
