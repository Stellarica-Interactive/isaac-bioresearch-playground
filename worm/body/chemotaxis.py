"""Chemosensation: a food gradient becomes current in the amphid neurons.

    food = FoodSource(x_m=0.15, y_m=0.0)
    sense = Chemosensation.build(cells=runtime.network.cell_ids)
    sense.step(food.concentration(nose_xy), dt_ms=4.17)
    runtime.inject_many(sense.currents())

The part that is easy to get wrong
----------------------------------

**The worm does not steer toward food.** It has no mechanism for sensing that food
is to its left. Its head is a couple of hundred microns across and the gradient
across that width is almost nothing.

What it does instead is a **biased random walk**. It runs roughly straight, then
performs a *pirouette* -- a reversal and a sharp turn -- and sets off in a new,
roughly random direction. The bias is entirely in *when* it pirouettes: with
concentration rising it suppresses them and keeps going; with concentration
falling it pirouettes sooner and tries another direction. That alone climbs a
gradient, with no steering anywhere in it (Pierce-Shimomura, Morse & Lockery 1999).

So the quantity these neurons report is **the time derivative of concentration**,
sampled as the animal moves, not the concentration itself. A model that injects
current proportional to absolute concentration is describing a different animal.
This module therefore adapts: each sensor tracks a slow running estimate of recent
concentration and responds to the departure from it.

Who senses what, and in which direction
---------------------------------------

**Measured:** AWA and AWC respond to volatile odorants and AWC^on is one of the two
cells whose biophysics we have (Bargmann, Hartwieg & Horvitz 1993). ASE carries
water-soluble attractants, and ASEL and ASER are functionally asymmetric --
unusual in any animal -- with ASEL responding to increases in salt and ASER to
decreases (Suzuki et al. 2008). AWC is famously an **OFF** cell: it fires on
odour *removal*, not arrival.

**Ours:** the shape of the gradient in the world, the adaptation time constant, and
the gain. All three are ASSUMED, and the gain is expressed as a depolarisation
rather than a current for the reasons in :mod:`common.neural.stimulus`.

What this cannot do yet
-----------------------

Produce chemotaxis. The behaviour needs runs and pirouettes, which need working
forward locomotion and reversals, and `docs/model_assumptions.md` §5C establishes
that we have neither. This is the sensory half: a gradient, a derivative, and
current arriving in the right cells with the right signs. Expect the same outcome
as touch (§5D) -- a correctly directed response too small to see.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

#: Depolarisation a full-scale concentration change should produce, mV. ASSUMED.
DEFAULT_RESPONSE_MV = 20.0

#: Adaptation time constant, ms. ASSUMED. It sets what "recently" means when the
#: cell compares now against a moment ago, and so the timescale of gradient the
#: animal can follow. Real adaptation in AWC is on the order of seconds.
DEFAULT_ADAPTATION_MS = 2000.0

#: How far the odour reaches, as a fraction of body length. ASSUMED: a real
#: bacterial lawn is a diffusion profile in a specific arena, not an exponential.
DEFAULT_DECAY_BODY_LENGTHS = 3.0


@dataclass(frozen=True)
class Sensor:
    """One chemosensory cell, and which direction of change excites it."""

    cell: str
    sign: float
    """+1 for a cell excited by *rising* concentration, -1 by falling."""

    note: str = ""


#: Sensors grouped by *modality*, because they must not be driven together.
#:
#: AWA/AWC answer volatile odorants and ASE water-soluble compounds -- different
#: molecules, sensed by different cells. Driving all six from one "food" scalar
#: makes the ON and OFF cells fight, and measurably inverts the answer: the AWC
#: pathway alone gives AIY-AIB of +0.148 mV on a rising gradient (run, the correct
#: response) while all six together give -0.057 mV (turn). See
#: docs/model_assumptions.md 5J.
MODALITIES: dict[str, tuple[str, ...]] = {
    "volatile": ("AWAL", "AWAR", "AWCL", "AWCR"),
    "soluble": ("ASEL", "ASER"),
    "awc": ("AWCL", "AWCR"),
}

#: The amphid chemosensory neurons this model uses, with response direction.
#:
#: Cell identities and their ON/OFF character are **measured**; which molecules
#: each prefers is measured too, though this model lumps "food" into one quantity
#: rather than distinguishing diacetyl from salt.
CHEMOSENSORS: tuple[Sensor, ...] = (
    Sensor("AWAL", +1.0, "volatile attractants; ON to odour onset"),
    Sensor("AWAR", +1.0, "volatile attractants; ON to odour onset"),
    Sensor("AWCL", -1.0, "volatile attractants; OFF cell -- fires on removal"),
    Sensor("AWCR", -1.0, "volatile attractants; OFF cell -- fires on removal"),
    Sensor("ASEL", +1.0, "water-soluble; responds to increases (Suzuki et al. 2008)"),
    Sensor("ASER", -1.0, "water-soluble; responds to decreases -- the ASE asymmetry"),
)


@dataclass(frozen=True)
class FoodSource:
    """A patch of bacteria, as a concentration field over the plane.

    Exponential falloff from a point. That is **not** what a bacterial lawn looks
    like -- a real one is a finite patch with a diffusion profile shaped by the
    arena -- but the behaviour depends on the gradient's sign and rough scale
    rather than its exact form, and an exponential makes the scale explicit.
    """

    x_m: float
    y_m: float
    decay_m: float = 0.1 * DEFAULT_DECAY_BODY_LENGTHS
    strength: float = 1.0

    def concentration(self, position_xy: np.ndarray) -> float:
        """Concentration at a point, peaking at ``strength`` on the source."""
        delta = np.asarray(position_xy, dtype=np.float64).reshape(-1)[:2] - np.array(
            [self.x_m, self.y_m]
        )
        return float(self.strength * np.exp(-np.linalg.norm(delta) / self.decay_m))


@dataclass
class Chemosensation:
    """Turns a concentration *history* into current in the amphid neurons.

    Stateful, and necessarily so: the response is to change, so the cell has to
    remember what it was smelling a moment ago.
    """

    sensors: tuple[Sensor, ...] = CHEMOSENSORS
    response_mv: float = DEFAULT_RESPONSE_MV
    adaptation_ms: float = DEFAULT_ADAPTATION_MS
    per_cell_pa: Mapping[str, float] | None = None
    """Current for a full-scale response, per cell. Built with
    :func:`common.neural.stimulus.solve_for_depolarisation`, because one current
    means different things to cells whose conductances differ."""

    baseline: float = field(default=0.0, init=False)
    deviation: float = field(default=0.0, init=False)
    _started: bool = field(default=False, init=False)

    @classmethod
    def build(
        cls,
        cells: tuple[str, ...] | None = None,
        *,
        modality: str | None = "awc",
        response_mv: float = DEFAULT_RESPONSE_MV,
        adaptation_ms: float = DEFAULT_ADAPTATION_MS,
        per_cell_pa: Mapping[str, float] | None = None,
    ) -> Chemosensation:
        """``cells`` restricts to sensors present in the network.

        ``modality`` selects which pathway the gradient drives, and defaults to
        ``"awc"`` rather than to everything. One concentration scalar is one
        molecule, and the AWC and ASE pathways answer different ones; driving both
        from it inverts the result (see :data:`MODALITIES`). Pass ``None`` for all
        six only when the stimulus really is meant to be every modality at once.
        """
        sensors = CHEMOSENSORS
        if modality is not None:
            if modality not in MODALITIES:
                raise KeyError(f"modality must be one of {sorted(MODALITIES)}, not {modality!r}")
            chosen = set(MODALITIES[modality])
            sensors = tuple(s for s in sensors if s.cell in chosen)
        if cells is not None:
            present = set(cells)
            sensors = tuple(s for s in sensors if s.cell in present)
        return cls(
            sensors=sensors,
            response_mv=response_mv,
            adaptation_ms=adaptation_ms,
            per_cell_pa=per_cell_pa,
        )

    @property
    def cells(self) -> tuple[str, ...]:
        return tuple(s.cell for s in self.sensors)

    def step(self, concentration: float, dt_ms: float) -> None:
        """Advance the adaptation state with a new sample.

        The baseline relaxes toward the current concentration with time constant
        ``adaptation_ms``, and the response is the gap between them. A cell sitting
        in a uniform field goes quiet however strong the field is, which is what
        adaptation is for and is why a real worm can follow a gradient across four
        orders of magnitude of absolute concentration.
        """
        if not self._started:
            # Start adapted, so arriving in a field is not itself a stimulus.
            self.baseline = float(concentration)
            self._started = True
        alpha = 1.0 - np.exp(-dt_ms / self.adaptation_ms)
        self.deviation = float(concentration) - self.baseline
        self.baseline += alpha * self.deviation

    def currents(self) -> dict[str, float]:
        """Current for each sensor, given the change it has just experienced."""
        out: dict[str, float] = {}
        for sensor in self.sensors:
            scale = (
                float(self.per_cell_pa[sensor.cell])
                if self.per_cell_pa is not None and sensor.cell in self.per_cell_pa
                else self.response_mv
            )
            value = scale * sensor.sign * self.deviation
            if value != 0.0:
                out[sensor.cell] = value
        return out

    def reset(self) -> None:
        self.baseline = 0.0
        self.deviation = 0.0
        self._started = False
