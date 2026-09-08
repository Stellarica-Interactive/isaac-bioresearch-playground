# Biology for software engineers

Written for someone who can read code but has not studied neuroscience. It covers
only what you need to understand the data currently in this repository, and it is
explicit about where the biology stops and our engineering starts.

---

## 1. What *C. elegans* is

A nematode worm about 1 mm long, transparent, living in rotting plant matter and
eating bacteria. It has been the workhorse of neuroscience since the 1970s for
one reason: it is small enough to be finished.

Two facts make it unique among animals:

- **The cell lineage is invariant.** Every adult hermaphrodite has exactly 959
  somatic cells, and they arise by exactly the same division pattern in every
  individual. Cell number 723 in your worm is the same cell as number 723 in mine.
  So neurons have *names*, not just types — `AVAL` is one specific cell, and it
  exists in every hermaphrodite that ever lived.
- **The nervous system is 302 neurons** in the adult hermaphrodite (plus 95 body
  wall muscles it controls). Not 302 types — 302 cells. That is small enough to
  have been reconstructed completely from electron micrographs.

Names encode anatomy, roughly. `AVAL` and `AVAR` are the left and right members
of the class `AVA`. `ASHL` is the left `ASH`. `VB7` is the seventh cell of the
`VB` ventral-cord motor neuron class. Muscles are `MDL01` (muscle, dorsal-left,
segment 1) through `MVR23`.

## 2. What a connectome is

A **connectome** is a wiring diagram: which cells physically contact which, and
how big those contacts are.

How it is made: fix and kill an animal, slice it into sections roughly 50 nm
thick, image every section with an electron microscope, then trace each cell's
processes through the whole stack by hand and mark every place a synapse is
visible. This takes years per animal.

What comes out is, per pair of cells, a count. In the Witvliet datasets that
count is the **number of presynaptic active zones** — visible release sites —
agreed on by at least two of three independent human annotators.

### What a connectome is *not*

This distinction is the backbone of the whole project, so it is worth being blunt:

| | In a connectome? |
|---|---|
| Which cells exist | ✅ yes |
| Which cells contact which | ✅ yes |
| How many synapses per contact | ✅ yes |
| Chemical vs electrical contact | ✅ yes |
| Whether a synapse **excites or inhibits** | ❌ no — see §4 |
| Synaptic **strength** in physiological units | ❌ no — a count is not a conductance |
| What neurotransmitter a cell releases | ❌ no — separate experiments |
| Whether a cell is sensory or motor | ❌ no — separate experiments |
| The animal's neural **activity** when it died | ❌ no — it was fixed and dead |
| Anything the animal had **learned** | ❌ no |
| Signalling that leaves no visible structure | ❌ no — see §6 |

A connectome constrains what a nervous system can do without determining what it
does. This is why the repository keeps anatomy and annotation in separate files
and makes every annotated field name its own source.

## 3. Chemical synapses and gap junctions

Two ways neurons connect, and they behave completely differently.

**Chemical synapse.** The presynaptic cell releases a neurotransmitter into a
small gap; the postsynaptic cell has receptors that respond. Directed: signal
flows one way. Slower, and — crucially — the *sign* is decided by the receiving
cell, not the sending one.

**Electrical synapse (gap junction).** A protein channel physically joins two
cells' cytoplasm, so ions flow directly between them. Undirected and fast, and it
tends to *equalize* the voltage of the two cells rather than push one in a
particular direction.

In code this shows up as: chemical edges are directed and stored once per ordered
pair; gap junctions are undirected and stored once per unordered pair. Getting
that convention wrong silently doubles weights — which is why
`docs/model_assumptions.md` §3 spells it out and the test suite pins it.

> **Oddity in the real data:** every Witvliet dataset contains a handful of gap
> junctions from a cell *to itself* — two processes of the same neuron joined to
> each other. We preserve these rather than dropping them as noise. They are also
> why the published weight totals do not always divide by two.

## 4. Why synaptic sign is missing, and why that matters so much

If neuron A synapses onto neuron B, does A make B more active or less?

Electron microscopy cannot tell you. The image shows vesicles and a release site;
it does not show which receptor protein sits on the other side. And that receptor
is what decides. The *same* neurotransmitter can excite one cell and inhibit
another, because they express different receptors. Acetylcholine excites body
wall muscle through one channel and inhibits other cells through another.

So to get sign you need to combine:

1. which transmitter the presynaptic cell releases (→ Wang et al. 2024, in this
   repo), and
2. which receptors the postsynaptic cell expresses (→ single-cell RNA sequencing).

Fenyves et al. (2020) did exactly this and could predict sign for roughly
two-thirds of chemical synapses. The rest remain unknown.

This is why sign is a *prediction*, never a measurement, and why the loader does
**not** apply a polarity overlay by default. A neural simulation is extremely
sensitive to whether a connection is + or −; letting a prediction slip in
unnoticed would be the fastest way for this project to produce impressive-looking
nonsense.

## 5. Graded potentials: the worm does not mostly spike

Textbook neurons fire **action potentials** — all-or-nothing voltage spikes that
travel down an axon. Digital, in a sense: a spike either happens or it does not,
and information is in the timing.

Most *C. elegans* neurons do not do this. They are **graded**: membrane voltage
varies continuously, and transmitter release varies continuously with it. Analog
rather than digital. This is possible because the worm is tiny — a neuron is
short enough that voltage is nearly uniform along it, so there is no long cable
needing a regenerating spike to carry a signal down it.

Consequences for modelling: a spiking-neural-network framework is the wrong tool
here. What you want is a continuous dynamical system — each neuron a state
variable that leaks toward rest and is pushed around by its inputs.

**But "the worm has no action potentials" is an oversimplification**, and this
repository will not repeat it as fact:

- **AWA**, an olfactory neuron, fires genuine calcium-based all-or-nothing action
  potentials (Liu et al. 2018).
- **AVL** fires compound action potentials during the defecation cycle.
- **RMD** motor neurons show **plateau potentials** — long-lasting, all-or-none,
  terminable by a negative current pulse.

So the honest statement is: *graded transmission dominates, and a graded model is
the right starting point, but named exceptions exist and are documented.*

## 6. What is invisible to every connectome

*C. elegans* neurons also signal by releasing **neuropeptides** and **monoamines**
that diffuse and act on receptors far from any synapse. There is no
ultrastructural signature for this, so it appears in no connectome at all.

The scale is not marginal. Ripoll-Sánchez et al. (2023) mapped a neuropeptide
signalling network across the worm and found it to be dense and largely
non-overlapping with the wired connectome. A model built purely from synapses is
missing a whole communication layer, and would still be missing it even if the
synaptic data were perfect.

## 7. Sensory, interneuron, motor — and why it's fuzzier than it sounds

The conventional decomposition:

- **Sensory neurons** transduce something physical — a chemical concentration, a
  touch, temperature, oxygen — into voltage.
- **Motor neurons** synapse onto muscle.
- **Interneurons** sit between.

Real cells resist the boxes. Many are **polymodal**: `ASH` senses noxious
chemicals *and* high osmolarity *and* nose touch. Some cells were classified one
way by White et al. (1986) and another by later work. `CAN` is essential for the
worm to survive at all, makes no chemical synapses, and nobody knows what it does.

This repository handles that by letting a cell carry **several** roles, and by
recording `unknown` for the seven cells (`CANL/R`, `MI`, `MCL/R`, `NSML/R`) whose
published label does not resolve to a single role — rather than picking one.

## 8. The anatomy that decides which dataset you can use

The worm's body is a tube of muscle in four quadrants running head to tail. It
crawls by propagating a wave of dorsal/ventral bending backwards along its body.

The neurons that drive that wave — the **ventral nerve cord** motor neuron
classes `DA`, `DB`, `VA`, `VB`, `VC`, `DD`, `VD`, `AS` — run the length of the
animal, not in the head.

**This is why the dataset choice is load-bearing.** Witvliet et al. reconstructed
the *brain*: the nerve ring and nearby ganglia. Superb data, eight individuals,
modern methods — and it contains no ventral-cord motor neurons and only the
front eight muscle segments of each quadrant.

So:

> **Witvliet #7 cannot make a worm crawl.** The circuitry is not in the file.

It is excellent for head circuits, sensory processing, and comparing individuals.
For locomotion, a whole-animal reconstruction (Cook et al. 2019, White et al.
1986) is required. `Scope.HEAD` is recorded in the data itself so that code can
refuse rather than silently produce a paralysed worm.

## 9. Two behaviours we are aiming at

**Chemotaxis** (milestone W3). Put a worm in a chemical gradient and it climbs it
— but not by steering. It runs roughly straight, and when concentration is
*falling* it becomes more likely to stop and turn sharply (a "pirouette"). Biased
random walk, driven by the time-derivative of concentration. `ASE` neurons are
central to salt chemotaxis, and `ASEL`/`ASER` are functionally asymmetric — the
left and right members of one class respond differently, which is unusual and
useful.

**Touch withdrawal** (milestone W4). Touch the worm near the head and it reverses;
touch near the tail and it accelerates forward. Mediated by six touch receptor
neurons — `ALML/R`, `AVM` (anterior), `PLML/R`, `PVM` (posterior) — feeding the
command interneurons `AVA`/`AVD`/`AVE` (backward) and `AVB`/`PVC` (forward). One
of the best-characterized circuits in any animal, which makes it a good test: we
know what the answer should look like.

---

## Vocabulary

| Term | Meaning |
|---|---|
| **Connectome** | Wiring diagram from electron microscopy. Anatomy only. |
| **Functional connectivity** | Which neurons influence which *in a living animal*, measured by activity. A different thing; see Randi et al. 2023. |
| **Chemical synapse** | Directed, neurotransmitter-mediated. Sign set by the postsynaptic receptor. |
| **Gap junction** | Undirected direct electrical coupling. |
| **Graded potential** | Continuous voltage signalling, no spikes. Dominant in *C. elegans*. |
| **Action potential** | All-or-nothing voltage spike. Present in a few worm neurons. |
| **Active zone** | A synaptic release site. What the Witvliet weights count. |
| **Polymodal** | A neuron responding to several distinct stimulus types. |
| **Autapse** | A synapse from a neuron onto itself. Present in this data. |
| **Ventral nerve cord** | The main longitudinal nerve tract; where locomotor motor neurons live. |
| **Nerve ring** | The dense neuropil around the pharynx — the worm's "brain". |
| **Neuropil** | A region of packed processes where synapses form. |
| **L1 … L4** | Larval stages. Witvliet datasets 1–6 are larvae; 7 and 8 are adults. |

Full citations: [references.md](references.md).
