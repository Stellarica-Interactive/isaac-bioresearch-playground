# *Caenorhabditis elegans*

The proof-of-concept organism: 302 neurons in the adult hermaphrodite, completely
reconstructed, with decades of behavioural work to check a model against.

Background reading: [../docs/biology.md](../docs/biology.md).

## Layout

```
worm/
  loader.py        load("witvliet_2021_7") -> annotated Connectome. Start here.
  importers/
    naming.py      cell-name canonicalization; the single source of truth
    witvliet_2021.py
    sources.py     access to sources.toml and reference_totals.toml
  annotations/
    overlays.py    functional role, neurotransmitter, neuron class
  data/
    sources.toml            where each published file comes from, + SHA-256
    reference_totals.toml   published figures our parser is checked against
    cells.csv               canonical cell registry (629 cells)
    annotations/            sim_roles, neurotransmitters, neuron_classes
    circuits/               named circuits for the visualizer, with citations
    normalized/             committed derived data
    raw/                    git-ignored; fetch with tools/fetch_datasets.py
  tests/
```

## Milestones

| | | Status |
|---|---|---|
| **W0** | Data inspection: import, verify, explore, visualize | ✅ done |
| **W1** | Neural runtime: graded dynamics, checkpointing, stability tests | ⏸ blocked on [two decisions](../docs/model_assumptions.md#6-open-decisions-before-w1) |
| **W2** | Isaac Sim body: segmented, articulated, muscle-driven | not started |
| **W3** | Chemotaxis in a concentration field | not started |
| **W4** | Touch and escape | not started |
| **W5** | Lesions, individual comparison, environment manipulation, checkpoint/restore, body swap | not started |

W1 is deliberately not started. Choosing a neuron model and deciding what to do
about synaptic sign are biological judgements, not implementation details, and
writing code first would bury them.

## The dataset trap, stated once more

**Witvliet #7 and #8 are brain-only reconstructions.** They contain no
ventral-cord motor neurons (`DA`, `DB`, `VA`, `VB`, `VC`, `DD`, `VD`, `AS`) and
only body-wall muscle segments 1–8 of each quadrant.

The crawling gait is generated in the ventral cord. So:

- **Head circuits, sensory processing, individual comparison** → Witvliet, ideal.
- **Locomotion (W2–W4)** → needs a whole-animal reconstruction. Cook et al. 2019
  is the intended source and is **not yet imported**.

`Connectome.scope` records this, and `worm/tests/test_witvliet_import.py::TestAnatomicalScope`
asserts it, so the constraint fails loudly instead of producing a paralysed worm.

## Adding a dataset

1. Add an entry to `data/sources.toml` (url, filename, SHA-256, licence, citation).
   `tools/fetch_datasets.py --print-hashes` helps with the hash.
2. Add its published totals to `data/reference_totals.toml`, with the source of
   those numbers. A dataset with no oracle cannot be trusted, and
   `build_datasets.py` refuses to write one without `--allow-unverified`.
3. Write an importer in `importers/` returning an **unannotated** `Connectome`,
   and register it. Reuse `naming.canonical_cell_id`; if the new source spells a
   cell differently, add the rule there, not in the importer.
4. Cite it in `../docs/references.md` — a test enforces this.
5. `python tools/build_datasets.py --dataset <id>` then `python -m pytest`.

Importers produce anatomy only. Anything that is not "which cells contact which,
how, how often" belongs in an annotation overlay.
