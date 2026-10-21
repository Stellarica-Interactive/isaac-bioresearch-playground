"""Geometry for presenting the worm: a skin, a ground, and the connectome drawn on it.

Display only. Nothing here reaches the physics or the nervous system, so a run
looks different with it and measures identically. It is pure numpy so it can be
tested without Isaac Sim; ``worm/isaac/visuals.py`` turns it into USD.

* :func:`skin_geometry` and :func:`skin_topology` -- a smooth tube swept along
  the spine through the segment centres, so the animal reads as one body rather
  than a chain of capsules.
* :func:`agar_texture` and :func:`encode_png` -- a procedural agar-plate texture,
  written as a PNG with nothing but ``zlib``, so no binary enters the repository.
* :func:`connectome_layout`, :func:`connectome_edges`, :func:`edge_samples` and
  :func:`place`, gathered in :class:`ConnectomeDrawing` -- the connectome placed
  in the body: each neuron at a position along it, each synapse a path along it,
  bending with the animal.
* :class:`ConnectomePanel` -- the same, flat, for an on-screen panel.
* :func:`activity_colours` and :func:`linear` -- colours for activity, and their
  conversion to the linear values a renderer's ``displayColor`` expects.

**The layout is a display choice, not anatomy data.** The repository holds no
soma coordinates. Motor neurons sit where their own muscle synapses are
(:func:`~worm.body.neural_bridge.neuron_body_positions`, measured); every other
neuron is placed by the region its class lives in -- the nerve ring and head
ganglia, a handful of mid-body cells, the tail ganglia -- with a deterministic
scatter so the cells do not sit on top of each other.

Everything that runs per frame is vectorised and anything fixed for a run is
worked out once: at 240 Hz physics, a frame has milliseconds to spare.
"""

from __future__ import annotations

import hashlib
import struct
import zlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from worm.body.geometry import BodyPlan

if TYPE_CHECKING:
    from common.data.schemas import Connectome

# -- the body ------------------------------------------------------------------


def _catmull_rom(points: np.ndarray, per_span: int) -> np.ndarray:
    """A smooth curve through ``points``, ``per_span`` samples between each pair."""
    p = np.vstack([2 * points[0] - points[1], points, 2 * points[-1] - points[-2]])
    p0, p1, p2, p3 = p[:-3, None], p[1:-2, None], p[2:-1, None], p[3:, None]
    t = (np.arange(per_span) / per_span)[None, :, None]
    curve = 0.5 * (
        2 * p1
        + (-p0 + p2) * t
        + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t**2
        + (-p0 + 3 * p1 - 3 * p2 + p3) * t**3
    )
    return np.vstack([curve.reshape(-1, p.shape[1]), p[-2]])


def spine(
    centres: np.ndarray, plan: BodyPlan, per_segment: int = 4
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Points along the body from nose to tail tip, with their frames.

    ``centres`` are the segment centres, head first, shape ``(n_segments, 3)``.
    The curve is extended half a segment beyond the first and last centres so the
    skin reaches the actual tips. Returns ``(points, tangent, normal, fraction)``:
    ``normal`` lies in the plane of the body (the side ``_probe_contact`` calls
    dorsal), and ``fraction`` is the position along the body, 0 at the nose.
    """
    c = np.asarray(centres, dtype=np.float64)
    head = c[0] + 0.5 * (c[0] - c[1])
    tail = c[-1] + 0.5 * (c[-1] - c[-2])
    pts = _catmull_rom(np.vstack([head, c, tail]), per_segment)
    tangent = np.gradient(pts, axis=0)
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1, keepdims=True), 1e-12)
    normal = np.stack([-tangent[:, 1], tangent[:, 0], np.zeros(len(pts))], axis=1)
    normal /= np.maximum(np.linalg.norm(normal, axis=1, keepdims=True), 1e-12)
    seg = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    fraction = seg / max(seg[-1], 1e-12)
    return pts, tangent, normal, fraction


def body_radius(fraction: np.ndarray, plan: BodyPlan, tip: float = 0.18) -> np.ndarray:
    """The body's radius along its length: the segment profile, never quite zero."""
    f = np.clip(np.asarray(fraction, dtype=np.float64), 0.0, 1.0)
    return plan.max_radius_m * np.maximum(np.sin(np.pi * f) ** 0.4, tip)


def skin_topology(n_rings: int, n_around: int) -> tuple[np.ndarray, np.ndarray]:
    """Face vertex counts and indices for a tube of rings, capped at both ends.

    Vertices are the rings in order, ``n_around`` each, then one tip vertex per
    end. Fixed for a given body, so it is built once and only the points move.
    """
    counts, indices = [], []
    for r in range(n_rings - 1):
        for k in range(n_around):
            a = r * n_around + k
            b = r * n_around + (k + 1) % n_around
            counts.append(4)
            indices += [a, b, b + n_around, a + n_around]
    nose, tail = n_rings * n_around, n_rings * n_around + 1
    last = (n_rings - 1) * n_around
    for k in range(n_around):
        counts.append(3)
        indices += [nose, (k + 1) % n_around, k]
        counts.append(3)
        indices += [tail, last + k, last + (k + 1) % n_around]
    return np.array(counts, dtype=np.int32), np.array(indices, dtype=np.int32)


def skin_geometry(
    centres: np.ndarray, plan: BodyPlan, *, per_segment: int = 4, n_around: int = 16
) -> tuple[np.ndarray, np.ndarray]:
    """The skin's vertices for the body as it is now, and their outward normals.

    Vertices are as :func:`skin_topology` indexes them: the rings, then a tip at
    each end. One spine serves both, since this runs every frame.
    """
    pts, tangent, normal, fraction = spine(centres, plan, per_segment)
    radius = body_radius(fraction, plan)
    up = np.array([0.0, 0.0, 1.0])
    theta = np.linspace(0.0, 2 * np.pi, n_around, endpoint=False)
    around = np.cos(theta)[None, :, None] * normal[:, None, :] + np.sin(theta)[None, :, None] * up
    ring = (pts[:, None, :] + radius[:, None, None] * around).reshape(-1, 3)
    nose = pts[0] + (pts[0] - pts[1]) / max(np.linalg.norm(pts[0] - pts[1]), 1e-12) * radius[0]
    tail = (
        pts[-1] + (pts[-1] - pts[-2]) / max(np.linalg.norm(pts[-1] - pts[-2]), 1e-12) * radius[-1]
    )
    points = np.vstack([ring, nose, tail])
    normals = np.vstack([around.reshape(-1, 3), -tangent[0], tangent[-1]])
    return points, normals


def sphere(n_lat: int = 4, n_lon: int = 8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A unit sphere as a small mesh: ``(vertices, face counts, face indices)``.

    A pole at each end and ``n_lat - 1`` rings of ``n_lon`` between, faces wound
    to face outwards. Enough to read as round at the size a neuron is drawn, and
    small enough that three hundred of them update as one cheap mesh.
    """
    theta = 2 * np.pi * np.arange(n_lon) / n_lon
    rings = [
        np.column_stack(
            [
                np.sin(phi) * np.cos(theta),
                np.sin(phi) * np.sin(theta),
                np.full(n_lon, np.cos(phi)),
            ]
        )
        for phi in np.pi * np.arange(1, n_lat) / n_lat
    ]
    vertices = np.vstack([[0.0, 0.0, 1.0], *rings, [0.0, 0.0, -1.0]])
    bottom, last = len(vertices) - 1, 1 + (n_lat - 2) * n_lon
    counts: list[int] = []
    indices: list[int] = []
    for k in range(n_lon):
        counts.append(3)
        indices += [0, 1 + k, 1 + (k + 1) % n_lon]
    for r in range(n_lat - 2):
        upper, lower = 1 + r * n_lon, 1 + (r + 1) * n_lon
        for k in range(n_lon):
            k1 = (k + 1) % n_lon
            counts.append(4)
            indices += [upper + k, lower + k, lower + k1, upper + k1]
    for k in range(n_lon):
        counts.append(3)
        indices += [bottom, last + (k + 1) % n_lon, last + k]
    return vertices, np.array(counts, dtype=np.int32), np.array(indices, dtype=np.int32)


def skin_rings(n_segments: int, per_segment: int = 4) -> int:
    """How many rings :func:`skin_points` produces for a body of ``n_segments``."""
    return (n_segments + 1) * per_segment + 1


def ring_colours(segment_colours: np.ndarray, n_rings: int, n_around: int) -> np.ndarray:
    """Per-vertex colours for the skin, interpolated from one colour per segment."""
    seg = np.asarray(segment_colours, dtype=np.float64)
    at = np.linspace(0.0, len(seg) - 1, n_rings)
    rings = np.stack([np.interp(at, np.arange(len(seg)), seg[:, k]) for k in range(3)], axis=1)
    per_vertex = np.repeat(rings, n_around, axis=0)
    return np.vstack([per_vertex, rings[0], rings[-1]])


# -- the ground ----------------------------------------------------------------


def agar_texture(size: int = 512, seed: int = 0) -> np.ndarray:
    """An agar plate seen from above: a warm mid-tone, faintly mottled, finely grained.

    Smoothed random fields at a few scales over a warm translucent base, so the
    eye has texture to read motion against without a grid that competes with the
    animal. Deterministic in ``seed``. ``(size, size, 3)`` uint8, tileable.
    """
    rng = np.random.default_rng(seed)
    field = np.zeros((size, size))
    for octave, weight in ((8, 0.55), (32, 0.3), (128, 0.15)):
        octave = min(octave, size)
        coarse = rng.standard_normal((octave, octave))
        spectrum = np.fft.fft2(coarse)
        padded = np.zeros((size, size), dtype=complex)
        h = octave // 2
        padded[:h, :h] = spectrum[:h, :h]
        padded[:h, -h:] = spectrum[:h, -h:]
        padded[-h:, :h] = spectrum[-h:, :h]
        padded[-h:, -h:] = spectrum[-h:, -h:]
        layer = np.real(np.fft.ifft2(padded))
        field += weight * layer / max(np.abs(layer).max(), 1e-12)
    grain = rng.standard_normal((size, size)) * 0.04
    base = np.array([0.66, 0.60, 0.47])
    tone = 0.10 * field + grain
    rgb = np.clip(base[None, None, :] + tone[:, :, None] * np.array([1.0, 0.95, 0.8]), 0.0, 1.0)
    return (rgb * 255.0 + 0.5).astype(np.uint8)


def encode_png(rgb: np.ndarray) -> bytes:
    """An RGB or RGBA uint8 array as PNG bytes, with nothing but ``zlib``."""
    a = np.ascontiguousarray(rgb, dtype=np.uint8)
    h, w, channels = a.shape
    colour_type = {3: 2, 4: 6}[channels]
    raw = b"".join(b"\x00" + a[y].tobytes() for y in range(h))

    def chunk(kind: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + kind
            + data
            + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)
        )

    header = struct.pack(">IIBBBBB", w, h, 8, colour_type, 0, 0, 0)
    return (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", header)
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )


# -- the connectome --------------------------------------------------------------

#: Neurons outside the head, by class, and roughly where along the body their cell
#: bodies sit (0 nose, 1 tail). Display placement only; see the module docstring.
_MID_BODY = {
    "ALM": 0.45,
    "AVM": 0.40,
    "BDU": 0.30,
    "CAN": 0.50,
    "HSN": 0.55,
    "PDE": 0.72,
    "PVD": 0.75,
    "PVM": 0.62,
    "SDQ": 0.45,
}
_TAIL = {
    "ALN", "DVA", "DVB", "DVC", "LUA", "PDA", "PDB", "PHA", "PHB", "PHC", "PLM", "PLN",
    "PQR", "PVC", "PVN", "PVP", "PVQ", "PVR", "PVT", "PVW",
}  # fmt: skip
_PHARYNX = {"I1", "I2", "I3", "I4", "I5", "I6", "M1", "M2", "M3", "M4", "M5", "MC", "MI", "NSM"}


def _jitter(name: str, k: int) -> float:
    """A deterministic number in [-1, 1) from a cell's name, ``k`` selecting which."""
    digest = hashlib.sha256(f"{name}:{k}".encode()).digest()
    return int.from_bytes(digest[:4], "big") / 2**31 - 1.0


def connectome_layout(
    neurons: Sequence[str],
    classes: Mapping[str, str],
    motor_positions: Mapping[str, float],
    n_segments: int,
) -> dict[str, tuple[float, float, float]]:
    """Where each neuron is drawn: ``(fraction along the body, across, up)``.

    ``across`` and ``up`` run from -1 to 1 within the body's cross-section; positive
    ``across`` is the dorsal side. ``motor_positions`` are segment indices from
    :func:`~worm.body.neural_bridge.neuron_body_positions`, which places a motor
    neuron where its own muscle synapses are -- the one measured coordinate here.
    """
    out: dict[str, tuple[float, float, float]] = {}
    for name in neurons:
        cls = classes.get(name, "")
        if name in motor_positions and name[:2] in ("DA", "DB", "DD", "VA", "VB", "VD", "AS"):
            s = (motor_positions[name] + 0.5) / n_segments
            # The ventral cord runs along the ventral side.
            across, up = -0.55 + 0.15 * _jitter(name, 1), 0.3 * _jitter(name, 2)
        elif cls in _MID_BODY:
            s = _MID_BODY[cls] + 0.03 * _jitter(name, 0)
            across, up = 0.6 * _jitter(name, 1), 0.6 * _jitter(name, 2)
        elif cls in _TAIL or name in _TAIL:
            s = 0.92 + 0.04 * _jitter(name, 0)
            across, up = 0.55 * _jitter(name, 1), 0.55 * _jitter(name, 2)
        elif cls in _PHARYNX or name in _PHARYNX:
            s = 0.06 + 0.04 * _jitter(name, 0)
            across, up = 0.25 * _jitter(name, 1), 0.25 * _jitter(name, 2)
        else:
            # The head: nerve ring and ganglia, the bulk of the nervous system.
            s = 0.11 + 0.05 * _jitter(name, 0)
            across, up = 0.7 * _jitter(name, 1), 0.7 * _jitter(name, 2)
        out[name] = (float(np.clip(s, 0.0, 1.0)), float(across), float(up))
    return out


def connectome_edges(
    g_syn: np.ndarray, cell_ids: Sequence[str], keep: Sequence[str], max_edges: int = 600
) -> list[tuple[str, str, float]]:
    """The strongest chemical connections among ``keep``, as ``(pre, post, nS)``.

    ``g_syn`` is the network's ``[post, pre]`` conductance matrix. A connectome has
    thousands of edges; drawn all at once they are a smear, so only the strongest
    are kept.
    """
    index = {c: i for i, c in enumerate(cell_ids)}
    idx = np.array([index[c] for c in keep if c in index])
    sub = g_syn[np.ix_(idx, idx)]
    posts, pres = np.nonzero(sub)
    weights = sub[posts, pres]
    order = np.argsort(-weights)[:max_edges]
    names = [cell_ids[i] for i in idx]
    return [(names[pres[k]], names[posts[k]], float(weights[k])) for k in order]


#: Activity colours at rest, half-way and fully active, sRGB. Through a warm rose
#: rather than straight from blue to gold, a blend that passes through grey.
ACTIVITY_STOPS = ((0.16, 0.26, 0.72), (0.86, 0.30, 0.42), (1.0, 0.84, 0.18))


def activity_colours(
    level: np.ndarray, stops: Sequence[Sequence[float]] = ACTIVITY_STOPS
) -> np.ndarray:
    """sRGB colours for activity ``level`` from 0 (rest) to 1 (fully active), ``(n, 3)``.

    Low levels are lifted (a 0.6 power) so a cell part-way active is visibly not
    at rest.
    """
    x = np.clip(np.asarray(level, dtype=np.float64), 0.0, 1.0) ** 0.6
    table = np.asarray(stops, dtype=np.float64)
    position = x * (len(table) - 1)
    k = np.minimum(position.astype(int), len(table) - 2)
    w = (position - k)[:, None]
    return table[k] * (1.0 - w) + table[k + 1] * w


def linear(rgb: np.ndarray) -> np.ndarray:
    """sRGB colours as linear values, the space RTX reads ``displayColor`` in.

    Set unconverted, a colour picked by eye renders paler and greyer than chosen
    -- the skin's activity tint washed out to near-white that way (found by
    screenshot). Textures are unaffected: the renderer decodes those itself.
    """
    c = np.clip(np.asarray(rgb, dtype=np.float64), 0.0, 1.0)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


#: The panel's background, and a synapse's colour before activity tints it, sRGB.
_PANEL_BACKGROUND = (0.05, 0.06, 0.09)
_PANEL_SYNAPSE = (0.45, 0.55, 0.75)


@dataclass(frozen=True)
class ConnectomePanel:
    """The connectome flat, for an on-screen panel.

    Horizontal is position along the body, nose at the left; vertical spreads the
    cells within their region. Cells are dots coloured by their activity; synapses
    are faint curves that take on their presynaptic cell's colour as it
    activates; the background is translucent and dark. Which pixels each cell and
    synapse covers is fixed, so :meth:`build` works it out once and :meth:`draw`
    only colours them.
    """

    width: int
    height: int
    #: The empty panel, ``(height * width, 4)`` uint8.
    background: np.ndarray
    #: Every pixel any synapse crosses, as ``y * width + x``.
    line_pixel: np.ndarray
    #: Per synapse sample: which synapse, and which of ``line_pixel`` it lands on.
    sample_edge: np.ndarray
    sample_slot: np.ndarray
    #: Per synapse: its presynaptic cell, and its opacity with that cell at rest.
    edge_pre: np.ndarray
    edge_alpha: np.ndarray
    #: Per dot pixel: which cell, and the pixel.
    dot_cell: np.ndarray
    dot_pixel: np.ndarray

    @classmethod
    def build(
        cls,
        layout: Mapping[str, tuple[float, float, float]],
        names: Sequence[str],
        edges: Sequence[tuple[str, str, float]],
        *,
        width: int = 640,
        height: int = 200,
    ) -> ConnectomePanel:
        margin = 12
        where = np.array([layout[n] for n in names], dtype=np.float64).reshape(-1, 3)
        x = margin + where[:, 0] * (width - 2 * margin)
        y = height / 2 + (0.45 * where[:, 1] + 0.25 * where[:, 2]) * (height / 2 - margin)
        order = {n: i for i, n in enumerate(names)}
        pre = np.array([order[a] for a, _, _ in edges], dtype=np.int64)
        post = np.array([order[b] for _, b, _ in edges], dtype=np.int64)
        g = np.array([w for _, _, w in edges], dtype=np.float64)
        bow = np.array([_jitter(f"{a}>{b}", 3) for a, b, _ in edges]) * 0.25 * (x[post] - x[pre])
        t = np.linspace(0.0, 1.0, 24)[None, :]
        ex = (1 - t) * x[pre, None] + t * x[post, None]
        ey = (1 - t) * y[pre, None] + t * y[post, None] - 4 * t * (1 - t) * bow[:, None] * 0.2
        ex, ey = np.rint(ex).astype(np.int64), np.rint(ey).astype(np.int64)
        on = (ex >= 0) & (ex < width) & (ey >= 0) & (ey < height)
        yy, xx = np.mgrid[-3:4, -3:4]
        disc = (xx**2 + yy**2) <= 7
        dx = np.rint(x).astype(np.int64)[:, None] + xx[disc][None, :]
        dy = np.rint(y).astype(np.int64)[:, None] + yy[disc][None, :]
        dot_on = (dx >= 0) & (dx < width) & (dy >= 0) & (dy < height)
        line_pixel, sample_slot = np.unique((ey * width + ex)[on], return_inverse=True)
        background = np.empty((width * height, 4))
        background[:, :3] = _PANEL_BACKGROUND
        background[:, 3] = 0.82
        return cls(
            width=width,
            height=height,
            background=_bytes(background),
            line_pixel=line_pixel,
            sample_edge=np.broadcast_to(np.arange(len(pre))[:, None], ex.shape)[on],
            sample_slot=sample_slot.ravel(),
            edge_pre=pre,
            edge_alpha=np.minimum(0.18 + 0.1 * g, 0.5),
            dot_cell=np.broadcast_to(np.arange(len(names))[:, None], dx.shape)[dot_on],
            dot_pixel=(dy * width + dx)[dot_on],
        )

    def draw(self, level: np.ndarray) -> np.ndarray:
        """The panel, ``(height, width, 4)`` uint8, for each cell's activity
        ``level`` -- 0 rest to 1 fully active, in the ``names`` order it was built
        with."""
        level = np.clip(np.asarray(level, dtype=np.float64), 0.0, 1.0)
        colours = activity_colours(level)
        glow = np.minimum(1.0, 1.5 * level[self.edge_pre])
        plain = np.asarray(_PANEL_SYNAPSE)
        edge_rgb = plain + (colours[self.edge_pre] - plain) * glow[:, None]
        alpha = (self.edge_alpha + 0.3 * glow)[self.sample_edge]
        # Every sample laid over the background at its own opacity, combined in
        # log space so the result does not depend on the order synapses are drawn.
        m = len(self.line_pixel)
        cover = 1.0 - np.exp(np.bincount(self.sample_slot, np.log1p(-alpha), minlength=m))
        weight = np.maximum(np.bincount(self.sample_slot, alpha, minlength=m), 1e-12)
        tint = np.stack(
            [
                np.bincount(self.sample_slot, alpha * edge_rgb[self.sample_edge, k], minlength=m)
                for k in range(3)
            ],
            axis=1,
        )
        line = np.asarray(_PANEL_BACKGROUND) * (1.0 - cover[:, None])
        line += tint / weight[:, None] * cover[:, None]
        img = self.background.copy()
        img[self.line_pixel, :3] = _bytes(line)
        img[self.dot_pixel, :3] = _bytes(colours[self.dot_cell])
        img[self.dot_pixel, 3] = 255
        return img.reshape(self.height, self.width, 4)


def _bytes(values: np.ndarray) -> np.ndarray:
    """Colour values 0..1 as uint8."""
    return (np.clip(values, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)


@dataclass(frozen=True)
class ConnectomeDrawing:
    """What to draw of the connectome: which cells, where, and which synapses.

    ``indices`` are the cells' positions in the network, for reading their
    activation each frame. ``rest`` and ``ceiling`` are the synaptic activation of
    a cell at its own rest and the most any cell can reach; the colour scale runs
    between them, so a colour means the same thing in every frame. Each synapse is
    drawn as ``points_per_edge`` samples along its path (:func:`edge_samples`).
    """

    names: tuple[str, ...]
    layout: dict[str, tuple[float, float, float]]
    edges: tuple[tuple[str, str, float], ...]
    indices: np.ndarray
    rest: float = 0.0
    ceiling: float = 1.0
    points_per_edge: int = 10

    @cached_property
    def body_coordinates(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(s, across, up)`` for every cell, then for every sample of every
        synapse's path. Fixed for the drawing, so worked out once."""
        cells = np.array([self.layout[n] for n in self.names], dtype=np.float64).reshape(-1, 3)
        s, across, up = edge_samples(self.layout, self.edges, points_per_edge=self.points_per_edge)
        return (
            np.concatenate([cells[:, 0], s]),
            np.concatenate([cells[:, 1], across]),
            np.concatenate([cells[:, 2], up]),
        )

    def positions(
        self,
        centres: np.ndarray,
        plan: BodyPlan,
        *,
        scale: float = 0.75,
        lift_m: float = 0.0,
        on_surface: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """World positions for the body as it is now: the cells, ``(n_cells, 3)``,
        and the synapse paths, ``(n_edges * points_per_edge, 3)``. One call to
        :func:`place` for both."""
        s, across, up = self.body_coordinates
        xyz = place(s, across, up, centres, plan, scale=scale, lift_m=lift_m, on_surface=on_surface)
        return xyz[: len(self.names)], xyz[len(self.names) :]

    def level(self, activation: np.ndarray) -> np.ndarray:
        """Each drawn neuron's activity, 0 at rest to 1 at the ceiling, from the
        network's activation array (every cell, network order)."""
        act = np.asarray(activation, dtype=np.float64)[self.indices]
        return np.clip((act - self.rest) / max(self.ceiling - self.rest, 1e-12), 0.0, 1.0)

    @classmethod
    def build(
        cls,
        connectome: Connectome,
        cell_ids: Sequence[str],
        g_syn: np.ndarray,
        plan: BodyPlan,
        *,
        max_edges: int = 600,
        activation_range: tuple[float, float] = (0.0, 1.0),
    ) -> ConnectomeDrawing:
        from common.data.schemas import CellCategory
        from worm.body.neural_bridge import neuron_body_positions

        present = set(cell_ids)
        neurons = tuple(
            c.id for c in connectome.cells if c.category is CellCategory.NEURON and c.id in present
        )
        classes = {c.id: (c.class_name or "") for c in connectome.cells}
        layout = connectome_layout(
            neurons, classes, neuron_body_positions(connectome, plan), plan.n_segments
        )
        edges = tuple(connectome_edges(g_syn, cell_ids, neurons, max_edges))
        index = {c: i for i, c in enumerate(cell_ids)}
        return cls(
            names=neurons,
            layout=layout,
            edges=edges,
            indices=np.array([index[n] for n in neurons], dtype=np.int64),
            rest=activation_range[0],
            ceiling=activation_range[1],
        )


def place(
    s: np.ndarray,
    across: np.ndarray,
    up: np.ndarray,
    centres: np.ndarray,
    plan: BodyPlan,
    *,
    scale: float = 0.75,
    lift_m: float = 0.0,
    on_surface: bool = False,
) -> np.ndarray:
    """World positions for body coordinates ``(s, across, up)``, ``(n, 3)``.

    ``on_surface`` puts everything just above the top of the skin instead of
    within the body, so that from a camera looking down it is drawn inside the
    body's outline and over it, which reads as inside without needing the skin to
    be transparent.
    """
    pts, _tangent, normal, fraction = spine(centres, plan)
    s = np.asarray(s, dtype=np.float64)
    k = np.clip(np.searchsorted(fraction, s), 1, len(pts) - 1)
    w = (s - fraction[k - 1]) / np.maximum(fraction[k] - fraction[k - 1], 1e-12)
    centre = pts[k - 1] + (pts[k] - pts[k - 1]) * w[:, None]
    nrm = normal[k - 1] + (normal[k] - normal[k - 1]) * w[:, None]
    nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-12)
    radius = body_radius(s, plan)
    world = centre + (radius * scale * np.asarray(across))[:, None] * nrm
    if on_surface:
        world[:, 2] = centre[:, 2] + 1.08 * radius
    else:
        world[:, 2] += radius * scale * np.asarray(up)
    world[:, 2] += lift_m
    return world


def edge_samples(
    layout: Mapping[str, tuple[float, float, float]],
    edges: Sequence[tuple[str, str, float]],
    *,
    points_per_edge: int = 10,
    bow: float = 0.35,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Each synapse as a path along the body: ``(s, across, up)`` at each of its
    ``points_per_edge`` samples, one synapse after another.

    A straight chord between two cells cuts across a bent body; an axon follows
    it. So the path runs along the body from one cell to the other, its
    cross-section coordinates blended between the two and bowed sideways by a
    deterministic amount. Placed with :func:`place` it bends with the animal;
    draw it as Catmull-Rom curves.
    """
    if not edges:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    t = np.linspace(0.0, 1.0, points_per_edge)[None, :]
    a = np.array([layout[pre] for pre, _, _ in edges], dtype=np.float64)
    b = np.array([layout[post] for _, post, _ in edges], dtype=np.float64)
    jitter = np.array([_jitter(f"{pre}>{post}", 5) for pre, post, _ in edges])
    sway = bow * jitter[:, None] * np.sin(np.pi * t)
    s = a[:, :1] + (b[:, :1] - a[:, :1]) * t
    across = np.clip(a[:, 1:2] + (b[:, 1:2] - a[:, 1:2]) * t + sway, -1.0, 1.0)
    up = a[:, 2:3] + (b[:, 2:3] - a[:, 2:3]) * t
    return s.ravel(), across.ravel(), up.ravel()


@dataclass
class Beside:
    """Where to draw a copy of the body so the camera overhead sees it beside the
    animal rather than on top of it.

    A fixed distance off the body's heading -- head minus tail, smoothed over
    ``tau_s`` so the offset does not sway with every undulation -- on whichever
    side faced +y, the top of the view, when it was first asked. It keeps that
    side as the animal turns, so it never jumps across the body.
    """

    distance_m: float
    tau_s: float = 3.0
    heading: np.ndarray | None = None
    side: float = 0.0
    _t: float | None = None

    def offset(self, centres: np.ndarray, t: float) -> np.ndarray:
        """The world offset ``(3,)`` for the body as it is at time ``t``, seconds."""
        chord = np.asarray(centres, dtype=np.float64)[0, :2] - np.asarray(centres)[-1, :2]
        norm = float(np.linalg.norm(chord))
        if norm > 1e-12 and np.isfinite(norm):
            chord = chord / norm
            if self.heading is None:
                self.heading = chord
            else:
                dt = max(t - self._t, 0.0) if self._t is not None else 0.0
                h = self.heading + (1.0 - np.exp(-dt / self.tau_s)) * (chord - self.heading)
                self.heading = h / max(float(np.linalg.norm(h)), 1e-12)
        self._t = t
        if self.heading is None:
            return np.zeros(3)
        left = np.array([-self.heading[1], self.heading[0]])
        if self.side == 0.0:
            self.side = 1.0 if left[1] >= 0.0 else -1.0
        xy = self.side * self.distance_m * left
        return np.array([xy[0], xy[1], 0.0])


def skin_tint(painted: np.ndarray | None, n_segments: int) -> np.ndarray:
    """Skin colours: a flesh tone, shifted towards the tint's orange or blue.

    ``painted`` is what :class:`~worm.isaac.stage.ActivityTint` gave each segment,
    which blends from its neutral grey; the shift away from that grey is added to
    the flesh tone, so the skin carries the same signal in its own colour.
    """
    flesh = np.array([0.74, 0.56, 0.47])
    if painted is None:
        return np.tile(flesh, (n_segments, 1))
    neutral = np.array([0.55, 0.55, 0.58])
    return np.clip(flesh[None, :] + 0.9 * (np.asarray(painted) - neutral[None, :]), 0.0, 1.0)
