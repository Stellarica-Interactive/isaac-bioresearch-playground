# Architecture

## The shape of the whole thing

The end goal is a closed sensorimotor loop:

```
   Isaac Sim environment
            |
        sensors                    physical quantity -> sensory signal
            v
     sensory neurons
            |
   connectome-derived neural runtime
            |
      motor neurons
            v
      muscle model                 neural activity -> actuation
            v
        body physics
            v
   changed environment  ---------> (back to sensors)
```

Nothing in the middle may short-circuit. If a behaviour appears, it has to arrive
through sensory neurons and leave through motor neurons. A line like
`if food_detected: turn_towards(food)` defeats the entire point of the project.

**Today only the first stage of the data pipeline exists.** Everything from
"neural runtime" rightwards is unbuilt; see [Deliberately absent](#deliberately-absent).

## The data pipeline, as built

```
worm/data/raw/*.xlsx                published EM reconstruction
      |                             MEASURED. git-ignored, SHA-256 pinned.
      | importer  (worm/importers/)
      v
worm/data/normalized/<id>/          our schema, PURE ANATOMY
      |  cells.csv, connections.csv, meta.json
      |  committed, sorted, LF, hash-verified
      |
      | overlays applied at load() -- never baked in
      |  (worm/annotations/)
      v
Connectome object                   anatomy + role + neurotransmitter + class
                                    every annotated field names its own source
```

### Why anatomy and annotation are separate files

1. **Different provenance, different update cadence.** A connectome parse is
   pinned forever to a published file. A neurotransmitter atlas gets revised.
   Coupling them means a correction upstream dirties an unrelated golden file.
2. **A verification test should test one thing.** "Witvliet #7 has 2202 chemical
   edges" must not also be implicitly asserting that the role table is intact.
3. **The same annotation applies to every dataset.** `AVAL` is an interneuron in
   Witvliet and in Cook. One table, no drift between importers.
4. **Coverage reporting falls out for free.** "Which of these 181 neurons has no
   neurotransmitter assignment" is real scientific information. An overlay
   *reports* what it could not match rather than silently defaulting it.

### Why every field names its source

`Cell` and `Connection` carry `field_sources: {field_name -> source_id}`,
resolvable in `Connectome.sources`. Combined with `Confidence`, that makes the
project's central requirement machine-checkable rather than a matter of good
intentions:

```
python tools/inspect_connectome.py --dataset witvliet_2021_7 gaps
```

lists everything in the loaded graph that this reconstruction did not measure,
with a citation for each. Cell *category* is in that list — electron microscopy
shows a cell exists, not that it is a neuron.

### Why synaptic polarity is opt-out by construction

`load()` applies `classes`, `sim` and `nt`. It does **not** apply a polarity
overlay, and none is implemented. Synaptic sign is not measured by anything; it
can only be predicted. A neural simulation is extremely sensitive to whether a
connection is + or −, so a prediction arriving unnoticed is the most likely route
to impressive-looking nonsense. `validate()` rejects any connection claiming a
`MEASURED` sign.

## Swapping datasets without touching the simulator

Every importer registers under a stable `dataset_id` and returns the same
`Connectome` type. Downstream code names a dataset, never a file format:

```python
from worm.loader import load
connectome, coverage = load("witvliet_2021_7")
```

`Connectome.scope` records anatomical extent, so code that needs a whole animal
can refuse a head-only dataset instead of silently producing a worm that cannot
move. That check is load-bearing: see [datasets.md](datasets.md).

## Layout

```
common/
  data/
    schemas.py     Cell, Connection, Connectome, the enums, Provenance. No I/O.
    io.py          Deterministic CSV + JSON serialization. Hash-verified.
    registry.py    Importer protocol and dataset_id registry.
    overlay.py     Annotation overlay protocol and coverage reporting.
    validate.py    Structural invariants. Shared by tests and the CLI.
  neural/
    graph.py       Organism-agnostic graph metrics over a Connectome.

worm/
  importers/       naming.py (canonicalization), witvliet_2021.py, sources.py
  annotations/     overlays.py
  data/            sources.toml, reference_totals.toml, cells.csv,
                   annotations/, circuits/, normalized/, raw/ (git-ignored)
  loader.py        load() -- the one entry point
  tests/

tools/
  fetch_datasets.py       download + verify SHA-256
  extract_annotations.py  dev-only: regenerate annotation tables from sources
  build_datasets.py       raw -> normalized (the only writer of committed data)
  inspect_connectome.py   the CLI
  visualize_connectome.py circuit diagrams
```

### Two commands that look similar and are not

| | Question | Fails when |
|---|---|---|
| `validate` | Is the graph internally well-formed? | Dangling endpoint, duplicate edge, non-canonical gap junction, unresolvable source id |
| `verify` | Do our totals match figures published by someone else? | The parser is self-consistent but wrong |

## Determinism

Committed derived data is compared byte-for-byte, so serialization is pinned:
LF newlines everywhere (`.gitattributes`), fixed column order, rows sorted by a
fixed key, JSON with sorted keys. `meta.json` records the SHA-256 of each CSV and
`load_connectome` verifies it — a hand-edited or CRLF-mangled file fails loudly.

## Deliberately absent

Requested in the original brief, not built, and **not stubbed**. An empty module
is worse than no module: it lets `from common.neural.runtime import X` typecheck
against nothing and implies decisions that have not been made.

| Path | Milestone | Blocked on |
|---|---|---|
| `common/neural/{runtime,neuron_models,synapses}.py` | W1 | Neuron-model and synaptic-sign decisions — [model_assumptions.md](model_assumptions.md) §6 |
| `common/isaac/{sensors,actuators,simulation_loop}.py` | W2 | W1 |
| `common/visualization/neural_activity.py` | W2 | Something to visualize |
| `worm/{neural,body,isaac,experiments}/` | W1–W5 | as above |
| `fly/` beyond a README | Part II | The worm loop closing first |
| `tools/benchmark_runtime.py` | W1 | There is no runtime to benchmark |

`fly/` in particular is left as a single README on purpose. Building eight empty
directories now would force `common/` into a worm-shaped abstraction never
validated against a second organism, which is exactly how a shared layer ends up
wrong. It gets built when a FlyWire or BANC dataset is actually being parsed.

## Dependencies

Runtime: `numpy`, `networkx`, `openpyxl`. That is the whole list. Downloads use
stdlib `urllib`, CLIs use `argparse`, manifests use stdlib `tomllib`. Optional:
`matplotlib` (`[viz]`), `pytest`/`ruff`/`mypy` (`[dev]`).

No `cect` at runtime. The OpenWorm Connectome Toolbox is used as a *mirror* for
published files and as a reference for what those files contain; we parse them
ourselves, in about fifty lines, so that the parse is ours to understand and
verify.
