# *Drosophila melanogaster* — deliberately not started

This directory holds a README and nothing else, on purpose.

## Why it is empty

The fly is Part II of this project, and it is not "the worm, but larger". It
differs in ways that will change the design:

| | *C. elegans* | *Drosophila* |
|---|---|---|
| Neurons | 302, individually named and invariant | ~140,000, typed rather than individually named |
| Signalling | predominantly **graded**, mostly non-spiking | substantially **spiking** — a different neuron model |
| Connectome | one file, ~2,500 edges | ~50 million synapses; subsetting is mandatory |
| Body | one soft segmented tube | six legs, two wings, flight aerodynamics |
| Neurotransmitters | measured per cell by reporter alleles | **predicted** from EM morphology, with an error rate |

Creating the eight-directory skeleton now would force `common/` into a
worm-shaped abstraction that has never been validated against a second organism.
That is exactly how a shared layer ends up wrong. The tree gets built when a fly
dataset is actually being parsed, and `common/` gets generalized then, against
real requirements rather than guessed ones.

## Prerequisite

The worm's closed sensorimotor loop working end to end — milestones W1 through
W4. Until a real connectome drives a real body in Isaac Sim and something
behaviourally interesting comes out, scaling up 500× is premature.

## Groundwork for when it starts (F0)

Not yet done; recorded so it is not forgotten.

**Choose a dataset, with reasons.** Do not default to the biggest one.

| | FAFB / FlyWire | BANC |
|---|---|---|
| Extent | female adult fly **brain**, ~140k neurons | brain **plus ventral nerve cord**, ~158k neurons |
| Motor output | descending neurons present; the motor neurons they drive are **not** | includes the nerve cord where leg and wing motor neurons live |
| Maturity | extensively annotated and published | newer |

The worm taught this lesson concretely: a brain-only reconstruction cannot
produce locomotion, however good it is. BANC is therefore the likely choice for
sensorimotor work — but compare annotation quality, sensory and motor coverage,
neurotransmitter predictions, and downloadable fields first, and write the
comparison down.

**Then, before simulating anything:** load a static download (not repeated Codex
queries), write our own normalized importer reusing `common/data/schemas.py`,
inspect the annotation system, and identify visual, olfactory, descending and
looming-sensitive pathways. Build tools that extract a *readable circuit* around
chosen cells — rendering 150,000 neurons at once produces nothing anyone can read.

**Neuron model:** research before choosing. Do not assume the worm's graded
model transfers. Fly neurons spike, and the right representation may be
event-based.

## Planned experiments

| | |
|---|---|
| **F3** | Phototaxis — orientation to light through a real visual pathway. Forbidden: `direction = normalize(light_pos - fly_pos)`. |
| **F4** | Odour navigation to fruit — static gradient, then a noisy plume, then wind, then plume loss and reacquisition. The fly never receives the source's coordinates. |
| **F5** | Looming escape — an approaching object, through experimentally identified looming-sensitive circuitry to descending neurons and takeoff. Forbidden: `if dangerous: escape()`. |
| **F6** | Competing behaviour — food approach interrupted by a looming threat. Which pathway wins, and does the other resume? |
| **F7** | A small naturalistic world with only biologically plausible sensory channels. No world coordinates, no target vectors, no navigation mesh. |

Starting points for the data: [FlyWire / Codex](https://codex.flywire.ai/),
Dorkenwald et al. *Nature* 634:124–138 (2024), Schlegel et al. *Nature*
634:139–152 (2024). See [../docs/references.md](../docs/references.md).
