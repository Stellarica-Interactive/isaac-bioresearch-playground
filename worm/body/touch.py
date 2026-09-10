"""Mechanosensation: a contact on the body becomes current in touch neurons.

    field = TouchField.build(plan)
    currents = field.currents(contact_fraction=0.2, ventral=False, force=1.0)
    runtime.inject_many(currents)

Why touch is worth doing before locomotion works
------------------------------------------------

``docs/model_assumptions.md`` §5C.9 establishes that our neuron model has no limit
cycle anywhere: every input drives it to a fixed point. That rules out a gait, but
it does **not** rule out a touch response, because a touch response is not a
rhythm. It is a transient -- a stimulus arrives, a signal propagates through the
circuit, motor output changes, and it decays -- and a relaxation model is exactly
the right shape for that.

So this is the one sensory experiment the oscillator problem does not block. What
it cannot deliver is the *behavioural* half: a real escape ends with the animal
reversing away, and reversing requires working locomotion. Expect to see the
circuit respond and the body fail to escape.

The circuit, which is measured
------------------------------

Gentle body touch in *C. elegans* is carried by six touch receptor neurons
(Chalfie & Sulston 1981; Chalfie et al. 1985, *J Neurosci* 5:956-964):

* **ALML, ALMR, AVM** -- anterior. Touch here and the animal reverses.
* **PLML, PLMR, PVM** -- posterior. Touch here and it accelerates forward.

They feed the command interneurons AVA/AVD/AVE (backward) and AVB/PVC (forward).
All of that wiring is in the Cook connectome already; none of it is added here.
This module only decides *which cells a contact reaches*, and how strongly.

What is measured and what is ours
---------------------------------

**Measured:** which cells are touch receptors, that ALM/AVM serve the anterior and
PLM/PVM the posterior, that AVM and PVM are ventral, and every synapse downstream.

**Ours:** the numeric extent of each receptive field along the body, and the
current per unit force. The processes' anatomical extents are documented
qualitatively -- ALM's process runs anteriorly from a cell body near mid-body,
PLM's runs anteriorly from the tail -- but turning that into fractions of body
length is a modelling choice, and turning contact into picoamps is an engineering
one. Both are ASSUMED and neither is tuned to produce a particular behaviour.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from worm.body.geometry import BodyPlan

#: Current injected into a touch neuron at the centre of its receptive field for a
#: unit contact, pA. **ASSUMED**, but not arbitrary.
#:
#: The first value here was 400 pA, chosen by analogy with the command drive and
#: the proprioceptive gain. That was wrong by more than an order of magnitude, and
#: instructively so: AVB tolerates 500 pA because it is heavily gap-junction
#: coupled and so has a large total conductance, whereas ALM has far less coupling
#: and the same current drove it to **+443 mV**. A neuron does not go to +443 mV.
#: The model has no spike mechanism and no upper bound of its own, so nothing
#: complained; the network simply carried on and eventually diverged.
#:
#: In this network ALM depolarises by roughly 1 mV per pA, so 20 pA gives a
#: receptor potential of about 20 mV. That is the right order for a graded
#: mechanoreceptor potential -- *C. elegans* touch cells respond to indentation
#: with graded receptor currents of order picoamps and depolarisations of tens of
#: millivolts (O'Hagan, Chalfie & Goodman 2005, Nat Neurosci 8:43-50) -- rather
#: than a number transplanted from a different kind of cell.
#:
#: It is still ASSUMED: the mapping from indentation depth to current is ours.
DEFAULT_TOUCH_CURRENT_PA = 20.0


@dataclass(frozen=True)
class Receptor:
    """One touch receptor neuron and the stretch of body it answers for.

    ``start`` and ``end`` are fractions of body length from the tip of the head, so
    they are independent of how many segments the body plan happens to use.
    """

    cell: str
    start: float
    end: float
    ventral_only: bool = False
    note: str = ""

    def overlap(self, fraction: float) -> float:
        """How strongly a contact at ``fraction`` falls in this field, 0 to 1.

        Graded rather than binary, peaking at the middle of the field and falling
        to zero at its edges. A step function would make the response depend on
        which side of an invented boundary a contact landed, which is a sharper
        claim than the anatomy supports.
        """
        if not (self.start <= fraction <= self.end):
            return 0.0
        mid = 0.5 * (self.start + self.end)
        half_width = 0.5 * (self.end - self.start)
        if half_width <= 0.0:
            return 0.0
        return float(max(0.0, 1.0 - abs(fraction - mid) / half_width))


#: The six gentle-touch receptors, with receptive fields along the body.
#:
#: Cell identities and their anterior/posterior division are **measured** (Chalfie
#: et al. 1985). The fractional extents are **ASSUMED**, chosen to reflect the
#: documented anatomy: ALM processes run forward from cell bodies near mid-body, so
#: the anterior field covers roughly the front half; PLM processes run forward from
#: the tail, so the posterior field covers roughly the back half; the two overlap
#: around mid-body, which is where the animal is famously least touch-sensitive
#: because a contact there recruits both and the responses oppose.
TOUCH_RECEPTORS: tuple[Receptor, ...] = (
    Receptor("ALML", 0.0, 0.55, note="anterior, lateral left"),
    Receptor("ALMR", 0.0, 0.55, note="anterior, lateral right"),
    Receptor(
        "AVM",
        0.05,
        0.50,
        ventral_only=True,
        note="anterior, ventral; born post-embryonically",
    ),
    Receptor("PLML", 0.45, 1.0, note="posterior, lateral left"),
    Receptor("PLMR", 0.45, 1.0, note="posterior, lateral right"),
    Receptor(
        "PVM",
        0.50,
        0.95,
        ventral_only=True,
        note="posterior, ventral; does not by itself drive the gentle-touch response",
    ),
)

#: Harsh touch and nose touch are different circuits with different receptors. They
#: are listed so that a future experiment does not silently reuse the gentle-touch
#: mapping for a stimulus it was never meant to describe.
HARSH_TOUCH_CELLS: tuple[str, ...] = ("PVDL", "PVDR")
NOSE_TOUCH_CELLS: tuple[str, ...] = ("ASHL", "ASHR", "FLPL", "FLPR")


@dataclass(frozen=True)
class TouchField:
    """Turns a contact somewhere on the body into current in touch neurons."""

    plan: BodyPlan
    receptors: tuple[Receptor, ...] = TOUCH_RECEPTORS
    current_pa: float = DEFAULT_TOUCH_CURRENT_PA

    @classmethod
    def build(
        cls,
        plan: BodyPlan,
        *,
        current_pa: float = DEFAULT_TOUCH_CURRENT_PA,
        cells: tuple[str, ...] | None = None,
    ) -> TouchField:
        """``cells`` restricts to receptors actually present in the network."""
        receptors = TOUCH_RECEPTORS
        if cells is not None:
            present = set(cells)
            receptors = tuple(r for r in receptors if r.cell in present)
        return cls(plan=plan, receptors=receptors, current_pa=current_pa)

    @property
    def cells(self) -> tuple[str, ...]:
        return tuple(r.cell for r in self.receptors)

    def currents(
        self, *, contact_fraction: float, ventral: bool = False, force: float = 1.0
    ) -> dict[str, float]:
        """Current for each receptor, given one contact.

        ``contact_fraction`` is position along the body, 0 at the head tip and 1 at
        the tail tip. ``force`` scales the whole response linearly, which is an
        engineering convenience rather than a claim about mechanotransduction.
        """
        out: dict[str, float] = {}
        for r in self.receptors:
            if r.ventral_only and not ventral:
                continue
            weight = r.overlap(contact_fraction)
            if weight > 0.0:
                out[r.cell] = self.current_pa * weight * float(force)
        return out

    def segment_fraction(self, segment: int) -> float:
        """Body fraction at the centre of a segment, for contacts found per segment."""
        return (segment + 0.5) / self.plan.n_segments

    def currents_from_segments(
        self, touched: np.ndarray, *, ventral: np.ndarray | None = None
    ) -> dict[str, float]:
        """Combine contacts on several segments at once.

        ``touched`` holds a contact strength per segment, zero where untouched --
        which is what a probe of finite size against a bent body actually produces.
        Each receptor takes the **strongest** contact in its field rather than the
        sum, so pressing a wide object against the body does not report a larger
        stimulus than a sharp one merely by covering more segments.
        """
        touched = np.asarray(touched, dtype=np.float64)
        if touched.shape != (self.plan.n_segments,):
            raise ValueError(
                f"expected one contact strength per segment "
                f"({self.plan.n_segments},), got {touched.shape}"
            )
        is_ventral = (
            np.zeros(self.plan.n_segments, dtype=bool)
            if ventral is None
            else np.asarray(ventral, dtype=bool)
        )

        out: dict[str, float] = {}
        for segment in np.flatnonzero(touched > 0.0):
            contribution = self.currents(
                contact_fraction=self.segment_fraction(int(segment)),
                ventral=bool(is_ventral[segment]),
                force=float(touched[segment]),
            )
            for cell, value in contribution.items():
                out[cell] = max(out.get(cell, 0.0), value)
        return out
