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
| **W0** | Data inspection: import, verify, explore, visualize | ✅ done — Witvliet ×8 and Cook 2019, all verified against published totals |
| **W1** | Neural runtime: graded dynamics, checkpointing, stability tests | ⏸ blocked on the [neuron-model decision](../docs/model_assumptions.md#62-neuron-model); synaptic sign is [decided](../docs/model_assumptions.md#61-synaptic-sign--decided) |
| **W2** | Isaac Sim body: segmented, articulated, muscle-driven | not started |
| **W3** | Chemotaxis in a concentration field | not started |
| **W4** | Touch and escape | not started |
| **W5** | Lesions, individual comparison, environment manipulation, checkpoint/restore, body swap | not started |

W1 is deliberately not started. Choosing a neuron model and deciding what to do
about synaptic sign are biological judgements, not implementation details, and
writing code first would bury them.

## Choosing a dataset — this matters more than it sounds

**Witvliet #7 and #8 are brain-only reconstructions.** They contain no
ventral-cord motor neurons (`DA`, `DB`, `VA`, `VB`, `VC`, `DD`, `VD`, `AS`) and
only body-wall muscle segments 1–8 of each quadrant. The crawling gait is generated
in the ventral cord, so a locomotion model built on them would be a worm that cannot
move.

**Cook et al. 2019 is the whole animal**: 302 neurons, all 95 body wall muscles, the
full ventral cord, and real neuromuscular junctions.

| Task | Dataset |
|---|---|
| Locomotion, motor circuits, whole-animal graph | `cook_2019_herm` |
| Head and sensory circuits, best modern EM | `witvliet_2021_7` |
| Comparing two individuals | `witvliet_2021_7` vs `witvliet_2021_8` |
| Development across larval stages | `witvliet_2021_1` … `_8` |

Neither dominates: Cook is whole-animal but a single older reconstruction; Witvliet is
head-only but eight individuals traced with modern methods.

`Connectome.scope` records the anatomical extent in the data itself, and
`test_witvliet_import.py::TestAnatomicalScope` asserts it, so the constraint fails
loudly rather than silently producing a paralysed worm.

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
