# Data licences and attribution

**The MIT licence in [LICENSE](LICENSE) covers the code in this repository. It does
not cover the scientific data.**

That data is other people's work — in several cases years of manual electron-micrograph
tracing — and it comes with its own licences and its own citation requirements. This
file records them. Machine-readable equivalents live in
[`worm/data/sources.toml`](worm/data/sources.toml), and a test asserts the two agree.

---

## Summary

| What | Licence | Redistributed here? |
|---|---|---|
| Source code, tests, tools, documentation | MIT | yes |
| Published source datasets (`.xlsx`) | see below | **no** — fetched by script |
| Derived normalized connectomes (`worm/data/normalized/`) | see [Derived data](#derived-data) | yes |
| Derived annotation tables (`worm/data/annotations/`, `worm/data/cells.csv`) | see [Derived data](#derived-data) | yes |
| Generated figures (`docs/img/`) | MIT, but depict the datasets below | yes |

---

## Source datasets

### Witvliet et al. 2021 — eight *C. elegans* brain connectomes

> Witvliet D, Mulcahy B, Mitchell JK, Meirovitch Y, Berger DR, Wu Y, Liu Y, Koh WX,
> Parvathala R, Holmyard D, Schalek RL, Shavit N, Chisholm AD, Lichtman JW,
> Samuel ADT, Zhen M. **Connectomes across development reveal principles of brain
> maturation.** *Nature* 596:257–261 (2021). doi:10.1038/s41586-021-03778-8

- **Not redistributed here.** `tools/fetch_datasets.py` downloads the spreadsheets and
  verifies their SHA-256 against `worm/data/sources.toml`.
- Retrieved from the [openworm/ConnectomeToolbox](https://github.com/openworm/ConnectomeToolbox)
  mirror, which is MIT-licensed. The underlying scientific work is the paper above and
  must be cited as such.
- Primary distributions: [nemanode.org](http://nemanode.org),
  [bossdb.org/project/witvliet2020](https://bossdb.org/project/witvliet2020).

### Wang et al. 2024 — *C. elegans* neurotransmitter atlas

> Wang C, Vidal B, Sural S, Loer C, Aguilar GR, Merritt DM, Toker IA, Vogt MC,
> Cros CC, Hobert O. **A neurotransmitter atlas of *C. elegans* males and
> hermaphrodites.** *eLife* 12:RP95402 (2024). doi:10.7554/eLife.95402

- **Licence: CC BY 4.0.** Supplementary File 2 is the source of
  `worm/data/annotations/neurotransmitters.csv` and `neuron_classes.csv`.
- Not redistributed; fetched by script.

### Fenyves et al. 2020 — predicted synaptic polarity

> Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P.
> **Synaptic polarity and sign-balance prediction using gene expression data in the
> *Caenorhabditis elegans* chemical synapse neuronal connectome network.**
> *PLOS Computational Biology* 16(12):e1007974 (2020). doi:10.1371/journal.pcbi.1007974

- **Licence: CC BY 4.0.** S1 Data is the source of
  `worm/data/annotations/polarity_fenyves2020.csv`.
- Not redistributed; fetched by script.
- These are **predictions**, never measurements. They reach a simulation only
  through an opt-in overlay, and every row is tagged `PREDICTED`.

### Neuromuscular junction polarity — cited findings, no data

> Richmond JE, Jorgensen EM. **One GABA and two acetylcholine receptors function at
> the *C. elegans* neuromuscular junction.** *Nature Neuroscience* 2:791–797 (1999).
> doi:10.1038/12160

> McIntire SL, Jorgensen E, Kaplan J, Horvitz HR. **The GABAergic nervous system of
> *Caenorhabditis elegans*.** *Nature* 364:337–341 (1993). doi:10.1038/364337a0

Both are subscription articles. **No data from either is redistributed here.** The
`nmj` overlay encodes two facts each paper established — acetylcholine is excitatory
and GABA inhibitory at the body wall neuromuscular junction — applied to synapses we
already have from the connectome. Facts are not copyrightable; the papers are cited
because they are what makes the assignment defensible.

### WormAtlas cell listings, with Cook et al. 2019 groupings

> Cook SJ, Jarrell TA, Brittin CA, Wang Y, Bloniarz AE, Yakovlev MA, Nguyen KCQ,
> Tang LT-H, Bayer EA, Duerr JS, Bülow HE, Hobert O, Hall DH, Emmons SW.
> **Whole-animal connectomes of both *Caenorhabditis elegans* sexes.** *Nature*
> 571:63–71 (2019). doi:10.1038/s41586-019-1352-7

> [WormAtlas](https://www.wormatlas.org/) — Altun ZF, Herndon LA, Wolkow CA,
> Crocker C, Lints R, Hall DH (eds.), 2002–present.

- Source of `worm/data/cells.csv`, `annotations/sim_roles.csv` and
  `annotations/cell_descriptions.csv`, via the MIT-licensed OpenWorm
  Connectome Toolbox compilation (`cect/data/all_cell_info.csv`).
- WormAtlas content is provided for research and educational use; see
  <https://www.wormatlas.org/aboutus.htm>. Cite WormAtlas and Cook et al. 2019.

### OpenWorm *C. elegans* Connectome Toolbox

> <https://openworm.org/ConnectomeToolbox/> — [source](https://github.com/openworm/ConnectomeToolbox), MIT.

Used two ways, and **not** as a runtime dependency:

1. as a mirror for published data files, and
2. as an **independent** source of dataset totals. The figures in
   `worm/data/reference_totals.toml` come from their factsheets, which is what makes
   our verification meaningful — the oracle did not come from our own parser.

---

## Derived data

`worm/data/normalized/`, `worm/data/annotations/` and `worm/data/cells.csv` are
committed to this repository and are **derived** from the sources above: reformatted,
canonicalized and re-tabulated by the scripts in `tools/`, with no scientific content
added.

We make these available under the MIT licence **as far as our own contribution goes**,
which is the schema, the formatting and the transformation code. The underlying facts
belong to the cited works, and in some jurisdictions (notably the EU) database rights
may attach to the original compilations. If you redistribute the derived files:

- **cite the original publications**, not this repository, for the data itself;
- keep the `meta.json` provenance blocks and `source_ref` columns intact — every row
  and every field names where it came from, and that is the point;
- check the licence of each source above for your use case.

If you are an author of any of these datasets and would prefer we not redistribute a
derived form, please open an issue and we will remove it.

---

## How to cite

If this repository is useful to you, cite **the underlying data first**. The connectomes
represent years of painstaking manual reconstruction; the code here is a few weeks of
engineering on top of them.

For the software itself, see [CITATION.cff](CITATION.cff).
