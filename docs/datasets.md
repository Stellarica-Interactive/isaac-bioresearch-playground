# Datasets

What is imported, what each dataset contains, and — more importantly — what it
omits.

Fetch everything with:

```
python tools/fetch_datasets.py --organism worm --all
```

Raw files are not committed (see [model_assumptions.md](model_assumptions.md)
§3.5). Derived normalized data in `worm/data/normalized/` **is** committed, each
directory carrying a `meta.json` with the source URL, licence, citation and the
SHA-256 of the bytes it was built from.

---

## Witvliet et al. 2021 — eight individual brains across development

Eight isogenic *C. elegans* hermaphrodites reconstructed by serial-section
electron microscopy, from birth to adulthood. `witvliet_2021_1` … `witvliet_2021_8`.

**Weight semantics.** The number of presynaptic active zones between a pair of
cells, counted from electron micrographs and agreed upon by at least two of three
independent annotators.

| id | stage | cells | neurons | muscles | other | chemical edges / weight | gap junctions (stored) | gap junctions (published, directed) |
|---|---|---|---|---|---|---|---|---|
| `witvliet_2021_1` | L1 | 187 | 161 | 17 | 9 | 775 / 1296 | 83 | 164 / 206 |
| `witvliet_2021_2` | L1 | 194 | 162 | 21 | 11 | 986 / 1895 | 124 | 246 / 264 |
| `witvliet_2021_3` | L1 | 198 | 162 | 27 | 9 | 1012 / 2128 | 94 | 186 / 216 |
| `witvliet_2021_4` | L1 | 204 | 168 | 26 | 10 | 1136 / 2777 | 209 | 415 / 483 |
| `witvliet_2021_5` | L2 | 211 | 174 | 28 | 9 | 1515 / 4116 | 292 | 578 / 832 |
| `witvliet_2021_6` | L3 | 216 | 175 | 32 | 9 | 1525 / 4456 | 214 | 426 / 532 |
| **`witvliet_2021_7`** | **adult** | **222** | **181** | **32** | **9** | **2202 / 7467** | **291** | **576 / 794** |
| **`witvliet_2021_8`** | **adult** | **219** | **180** | **32** | **7** | **2186 / 7970** | **310** | **612 / 851** |

Every number in this table is asserted by the test suite against an independently
published figure. See "Verification" below.

### ⚠ Anatomical extent: brain only

Witvliet et al. reconstructed the nerve ring and immediately adjacent ganglia —
the worm's brain — not the whole animal. Consequences:

- **No ventral-cord motor neurons.** None of `DA`, `DB`, `VA`, `VB`, `VC`, `DD`,
  `VD`, `AS` are present.
- **Body wall muscles limited to segments 1–8** of each quadrant: exactly
  `MDL01–08`, `MDR01–08`, `MVL01–08`, `MVR01–08`. The adult has 95 body wall
  muscles; these datasets have 32.
- **Posterior sensory neurons largely absent** (`PLM`, `PVM`, `PHA`, `PHB`).
- Some cells appear on one side only — `CANR` is present in #7, `CANL` is not.
  That is the reconstruction volume, not biology.

> **The worm's crawling gait is generated in the ventral nerve cord.** These
> datasets do not contain it. A locomotion model built on them would be a worm
> that cannot move. This is recorded in the data as `Scope.HEAD` and asserted by
> `worm/tests/test_witvliet_import.py::TestAnatomicalScope`, so the failure is
> loud rather than silent.

Good for: head and sensory circuits, nerve-ring architecture, comparing
individuals, comparing developmental stages.

### Gap-junction self-connections

Every dataset contains gap junctions from a cell to itself — plausibly two
processes of the same neuron contacting each other:

| id | self-junctions | total self weight |
|---|---|---|
| 1, 2, 3, 6 | 2 | 2 |
| 4 | 3 | 3 |
| 5, 7 | 6 | 6 |
| 8 | 8 | **9** (`AIZR` has weight 2) |

Dataset 8's weight-2 self-junction is why its published electrical weight, 851,
is odd. Preserved, not cleaned away; see
[model_assumptions.md](model_assumptions.md) §3.1–3.2.

Witvliet #7 also contains **11 chemical autapses** — neurons synapsing onto
themselves: `ADEL`, `AIMR`, `ASER`, `ASHL`, `RMDDL`, `RMDDR`, `RMED`, `RMEL`,
`RMER`, `RMGL`, `SDQL`.

### Individual comparison, with a caveat

`witvliet_2021_7` and `witvliet_2021_8` are two different adult animals prepared
and traced the same way, which makes them the natural pair for asking how much
nervous systems differ between individuals.

```
python tools/inspect_connectome.py compare witvliet_2021_7 witvliet_2021_8
```

**Be careful with the answer.** A connection present in one and absent from the
other may be genuine individual variability, or a difference in reconstruction
volume, or annotator disagreement at the detection threshold. This data cannot
distinguish between those. The CLI prints that warning with the result.

---

## Annotation sources

Applied at load time, never baked into the connectome files.

| Source | Provides | Licence | Overlay |
|---|---|---|---|
| WormAtlas + Cook et al. 2019 groupings (via OpenWorm Connectome Toolbox) | cell category; sensory/inter/motor role | MIT (compilation) | `sim` |
| Wang et al. 2024, eLife 12:RP95402 | neurotransmitter; neuron class | CC BY 4.0 | `nt`, `classes` |
| Fenyves et al. 2020, PLOS Comp Biol | predicted synaptic sign | CC BY 4.0 | **not implemented — deliberate** |

Coverage on Witvliet #7: all three overlays annotate 181/181 neurons.

```
python tools/inspect_connectome.py --dataset witvliet_2021_7 coverage
```

---

## Verification

Our parser's totals are checked against numbers we did not produce: the OpenWorm
Connectome Toolbox factsheets, an independent analysis of the same source files.
They live in `worm/data/reference_totals.toml`, each with its citation.

```
python tools/inspect_connectome.py --dataset witvliet_2021_7 verify   # against published figures
python tools/inspect_connectome.py --dataset witvliet_2021_7 validate # internal consistency
```

The two are separate on purpose: a parser can be perfectly self-consistent and
still wrong.

All eight datasets currently agree with their oracle on all eight checked fields.

---

## Not yet imported

| Dataset | Why it matters | Status |
|---|---|---|
| **Cook et al. 2019** (herm.) | Whole animal: 302 neurons, 135 muscles, full ventral cord. **Required for any locomotion work.** | Planned next. Adjacency-matrix format, so a different importer. Note its published electrical figures (2883 entries, 23313 weight) are both odd, so its gap-junction matrix is not a clean symmetric doubling — expect to spend time reconciling that. |
| **White et al. 1986** | The original connectome. Historical baseline. | Available via the same mirror. |
| **Randi et al. 2023** | *Functional* connectivity — which neurons actually influence which in a living animal, by optogenetics plus whole-brain imaging. A different kind of measurement from anatomy. | Would be a strong cross-check on any dynamics model. |
| **Ripoll-Sánchez et al. 2023** | Neuropeptide signalling network — the communication layer invisible to EM. | In scope only if extrasynaptic signalling is modelled. |
| **Cook et al. 2019** (male) | Sexually dimorphic circuits. | Out of scope for now. |
