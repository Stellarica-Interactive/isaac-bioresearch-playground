"""The shape of the worm, and how muscles map onto it.

Deliberately free of any Isaac Sim or USD import. The body plan is a
specification that can be reasoned about, tested and sanity-checked without
launching a simulator; ``worm/isaac/`` turns it into a USD articulation.

Why 24 segments
---------------

Not because it is a round number. The adult hermaphrodite has 95 body wall
muscles in four quadrants running the length of the animal -- 24 dorsal-left,
24 dorsal-right, 23 ventral-left and 24 ventral-right. At 24 segments each one
receives exactly one muscle per quadrant, so ``MDL07`` drives segment 7's
dorsal-left and nothing has to be invented to make the mapping work. The one
gap, ``MVL24``, is a real gap in the animal rather than in the model.

Why the model is planar
-----------------------

*C. elegans* crawls **on its side**. Its dorsoventral bending therefore appears
as side-to-side waves in the plane of the agar, and a planar chain with one
rotational degree of freedom per joint captures the gait. This is the same
simplification Boyle, Berri and Cohen make in the standard neuromechanical
model. It discards roll, and it discards the small lateral component of real
crawling.

Why the body is scaled up
-------------------------

A real adult is about 1 mm long and 65 um wide. PhysX default tolerances are
tuned for objects around a metre, and sub-millimetre rigid bodies behave badly
without retuning. So the simulated body is scaled by
:data:`DEFAULT_LENGTH_SCALE` and the drag coefficients are scaled to match.
**This makes the dynamics dimensionally consistent but not quantitatively the
animal's** -- see ``docs/model_assumptions.md``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

#: Segments along the body, head to tail. One per body wall muscle row.
DEFAULT_SEGMENTS = 24

#: Real adult hermaphrodite dimensions, in metres.
REAL_LENGTH_M = 1.0e-3
REAL_MAX_RADIUS_M = 32.5e-6

#: Scale factor applied to every length. See the module docstring: PhysX is not
#: well behaved at sub-millimetre scale with default tolerances.
DEFAULT_LENGTH_SCALE = 100.0

#: Body wall muscle names are ``M`` + quadrant + two-digit row, e.g. ``MDL07``.
_MUSCLE_RE = re.compile(r"^M(?P<quadrant>[DV][LR])(?P<row>\d{2})$")

#: Quadrant order used for every muscle-indexed array in this package.
QUADRANTS = ("DL", "DR", "VL", "VR")

#: Which quadrants pull the body dorsally, and which ventrally. In the planar
#: model these are the two antagonists of each joint.
DORSAL = ("DL", "DR")
VENTRAL = ("VL", "VR")


@dataclass(frozen=True, slots=True)
class Segment:
    """One rigid link of the chain."""

    index: int
    """0 at the head, increasing toward the tail."""

    length_m: float
    radius_m: float
    centre_x_m: float
    """Position along the body axis of this segment's centre, head at x = 0."""

    @property
    def volume_m3(self) -> float:
        return float(np.pi * self.radius_m**2 * self.length_m)


@dataclass(frozen=True)
class BodyPlan:
    """A complete, simulator-independent description of the body.

    Segment radii follow a smooth taper that is fattest near the middle and
    narrows toward head and tail. That shape is qualitative -- it is chosen to
    look and behave like a worm, not fitted to measurements of one -- so it is an
    engineering simplification, recorded as such.
    """

    n_segments: int = DEFAULT_SEGMENTS
    length_scale: float = DEFAULT_LENGTH_SCALE
    density_kg_m3: float = 1000.0
    """Roughly water. The animal is mostly water, so this is a fair first guess."""

    @property
    def total_length_m(self) -> float:
        return REAL_LENGTH_M * self.length_scale

    @property
    def max_radius_m(self) -> float:
        return REAL_MAX_RADIUS_M * self.length_scale

    @property
    def segment_length_m(self) -> float:
        return self.total_length_m / self.n_segments

    @property
    def n_joints(self) -> int:
        """One fewer than the segments: joints sit between neighbours."""
        return self.n_segments - 1

    def radius_at(self, index: int) -> float:
        """Taper: full width through the middle, narrowing at both ends.

        ``sin`` raised to a small power gives a body that is blunt through the
        middle third and tapers over the last few segments, which is roughly the
        real profile. Qualitative, not fitted.
        """
        t = (index + 0.5) / self.n_segments
        return float(self.max_radius_m * np.sin(np.pi * t) ** 0.4)

    def segments(self) -> tuple[Segment, ...]:
        """Build the segment list. Not cached: 24 elements is not worth the
        lifetime complexity of memoising a method onto a frozen instance."""
        length = self.segment_length_m
        return tuple(
            Segment(
                index=i,
                length_m=length,
                radius_m=self.radius_at(i),
                centre_x_m=(i + 0.5) * length,
            )
            for i in range(self.n_segments)
        )

    def masses_kg(self) -> np.ndarray:
        return np.array([s.volume_m3 * self.density_kg_m3 for s in self.segments()])

    @property
    def total_mass_kg(self) -> float:
        return float(self.masses_kg().sum())

    def describe(self) -> str:
        return (
            f"{self.n_segments} segments, {self.n_joints} joints\n"
            f"  length {self.total_length_m * 1e3:.1f} mm, "
            f"max radius {self.max_radius_m * 1e3:.2f} mm "
            f"({self.length_scale:g}x real size)\n"
            f"  mass {self.total_mass_kg * 1e3:.3f} g at "
            f"{self.density_kg_m3:g} kg/m^3"
        )


# ---------------------------------------------------------------------------
# Muscles
# ---------------------------------------------------------------------------


def parse_muscle(name: str) -> tuple[str, int] | None:
    """``"MDL07"`` -> ``("DL", 7)``. ``None`` if it is not a body wall muscle."""
    m = _MUSCLE_RE.match(name)
    if not m:
        return None
    return m.group("quadrant"), int(m.group("row"))


def muscle_name(quadrant: str, row: int) -> str:
    return f"M{quadrant}{row:02d}"


def segment_of_muscle(name: str, plan: BodyPlan) -> int | None:
    """Which segment a muscle drives.

    Muscle rows are numbered from 1 at the head. With 24 segments the mapping is
    one-to-one; with any other segment count the rows are distributed evenly and
    the result is an interpolation, which is a modelling choice rather than
    anatomy. Prefer 24.
    """
    parsed = parse_muscle(name)
    if parsed is None:
        return None
    _, row = parsed
    index = int((row - 1) * plan.n_segments / 24)
    return min(index, plan.n_segments - 1)


def muscle_map(muscle_names: list[str], plan: BodyPlan) -> dict[str, list[int]]:
    """Group muscle names by quadrant, ordered by the segment each drives.

    Returns one list of segment indices per quadrant, aligned with
    ``muscle_names`` filtered to that quadrant. Muscles that are not body wall
    muscles are skipped, which is how pharyngeal and vulval muscles get excluded
    without a special case anywhere else.
    """
    out: dict[str, list[int]] = {q: [] for q in QUADRANTS}
    for name in muscle_names:
        parsed = parse_muscle(name)
        if parsed is None:
            continue
        quadrant, _ = parsed
        segment = segment_of_muscle(name, plan)
        if segment is not None:
            out[quadrant].append(segment)
    return out
