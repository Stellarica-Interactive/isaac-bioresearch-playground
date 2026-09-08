# References

Every dataset used, every paper a claim in this repository rests on, and the
licence each is distributed under.

Machine-readable provenance for the datasets lives in `worm/data/sources.toml`;
`worm/tests/test_docs.py` asserts that every source id there appears here, so the
two cannot drift apart.

---

## Datasets in use

### `witvliet_2021_1`, `witvliet_2021_2`, `witvliet_2021_3`, `witvliet_2021_4`, `witvliet_2021_5`, `witvliet_2021_6`, `witvliet_2021_7`, `witvliet_2021_8`

Witvliet D, Mulcahy B, Mitchell JK, Meirovitch Y, Berger DR, Wu Y, Liu Y, Koh WX,
Parvathala R, Holmyard D, Schalek RL, Shavit N, Chisholm AD, Lichtman JW,
Samuel ADT, Zhen M. **Connectomes across development reveal principles of brain
maturation.** *Nature* 596:257–261 (2021). doi:[10.1038/s41586-021-03778-8](https://doi.org/10.1038/s41586-021-03778-8)

Eight isogenic *C. elegans* brains reconstructed by serial-section electron
microscopy, L1 through adult. Datasets 7 and 8 are adults.

- Retrieved from the MIT-licensed [openworm/ConnectomeToolbox](https://github.com/openworm/ConnectomeToolbox) mirror of the published supplementary spreadsheets.
- Original data also at [nemanode.org](http://nemanode.org) and [bossdb.org/project/witvliet2020](https://bossdb.org/project/witvliet2020).
- **Brain-only reconstruction.** See [datasets.md](datasets.md).

### `cook_2019_herm`

Cook SJ, Jarrell TA, Brittin CA, Wang Y, Bloniarz AE, Yakovlev MA, Nguyen KCQ,
Tang LT-H, Bayer EA, Duerr JS, Bülow HE, Hobert O, Hall DH, Emmons SW.
**Whole-animal connectomes of both *Caenorhabditis elegans* sexes.** *Nature*
571:63–71 (2019). doi:[10.1038/s41586-019-1352-7](https://doi.org/10.1038/s41586-019-1352-7)

The whole-animal hermaphrodite connectome: 302 neurons, 135 muscles, 36 other
cells, **including the ventral nerve cord** and all 95 body wall muscles. The first
dataset here that could support a locomotion model.

- Supplementary Information 5, July 2020 correction. Adjacency matrices rather than
  an edge list, hence a separate importer.
- Retrieved from the MIT-licensed [openworm/ConnectomeToolbox](https://github.com/openworm/ConnectomeToolbox) mirror.
- Only the hermaphrodite sheets are imported; the male sheets are in the file but out
  of scope.

### `openworm_cect_all_cell_info` / `wormatlas_cook_2019_via_cect`

OpenWorm *C. elegans* Connectome Toolbox, `cect/data/all_cell_info.csv` — a
compilation of cell listings from [WormAtlas](https://www.wormatlas.org/NeuronNames.htm)
with the cell groupings of:

Cook SJ, Jarrell TA, Brittin CA, Wang Y, Bloniarz AE, Yakovlev MA, Nguyen KCQ,
Tang LT-H, Bayer EA, Duerr JS, Bülow HE, Hobert O, Hall DH, Emmons SW.
**Whole-animal connectomes of both *Caenorhabditis elegans* sexes.** *Nature*
571:63–71 (2019). doi:[10.1038/s41586-019-1352-7](https://doi.org/10.1038/s41586-019-1352-7)

Licence: MIT for the compilation. Provides cell category and coarse
sensory/interneuron/motor role. Derived tables: `worm/data/cells.csv`,
`worm/data/annotations/sim_roles.csv`.

### `wang_2024_nt_atlas`

Wang C, Vidal B, Sural S, Loer C, Aguilar GR, Merritt DM, Toker IA, Vogt MC,
Cros CC, Hobert O. **A neurotransmitter atlas of *C. elegans* males and
hermaphrodites.** *eLife* 12:RP95402 (2024). doi:[10.7554/eLife.95402](https://doi.org/10.7554/eLife.95402)

CRISPR/Cas9 knock-in reporter alleles for every neurotransmitter pathway gene,
scored across all 302 hermaphrodite neurons. Licence: CC BY 4.0. Supplementary
File 2 (hermaphrodite) provides neurotransmitter assignment and neuron class.
Derived tables: `worm/data/annotations/neurotransmitters.csv`,
`neuron_classes.csv`.

---

## Cited but not yet imported

**Fenyves BG, Arnold C, Gerencsér VG, Vörös P, Pollner P, Csermely P, Korcsmáros T.
Synaptic polarity and sign-balance prediction using gene expression data in the
*Caenorhabditis elegans* chemical synapse neuronal connectome network.**
*PLOS Computational Biology* 16(12):e1007974 (2020).
doi:[10.1371/journal.pcbi.1007974](https://doi.org/10.1371/journal.pcbi.1007974)
Licence CC BY 4.0. Predicts excitatory/inhibitory sign for ~2/3 of chemical
synapses from presynaptic transmitter and postsynaptic receptor gene expression.
Would be the source for an opt-in polarity overlay —
[model_assumptions.md](model_assumptions.md) §6.1.

**Randi F, Sharma AK, Dvali S, Leifer AM. Neural signal propagation atlas of
*Caenorhabditis elegans*.** *Nature* 623:406–414 (2023).
doi:[10.1038/s41586-023-06683-4](https://doi.org/10.1038/s41586-023-06683-4)
*Functional* connectivity — which neurons actually influence which in a living
animal, by targeted optogenetics plus whole-brain calcium imaging. A different
kind of measurement from anatomy, and a strong future cross-check on any dynamics
model.

**Ripoll-Sánchez L, Watteyne J, Sun H, Fernandez R, Taylor SR, Weinreb A,
Bentley BL, Hammarlund M, Miller DM 3rd, Hobert O, Beets I, Vértes PE,
Schafer WR. The neuropeptidergic connectome of *C. elegans*.** *Neuron*
111:3570–3589 (2023). doi:[10.1016/j.neuron.2023.09.043](https://doi.org/10.1016/j.neuron.2023.09.043)
The extrasynaptic signalling layer, invisible to electron microscopy.

**White JG, Southgate E, Thomson JN, Brenner S. The structure of the nervous
system of the nematode *Caenorhabditis elegans*.** *Phil. Trans. R. Soc. Lond. B*
314:1–340 (1986). doi:[10.1098/rstb.1986.0056](https://doi.org/10.1098/rstb.1986.0056)
The original connectome. Still the baseline every later reconstruction is
compared against.

---

## Electrophysiology and neuron models

**Goodman MB, Hall DH, Avery L, Lockery SR. Active currents regulate sensitivity
and dynamic range in *C. elegans* neurons.** *Neuron* 20:763–772 (1998).
doi:[10.1016/S0896-6273(00)81014-4](https://doi.org/10.1016/S0896-6273(00)81014-4)
The source of essentially all measured passive properties: whole-cell capacitance
0.5–3 pF, input resistance in the gigaohm range, near-isopotentiality, and the
absence of classical action potentials in ASER and 42 further neurons.

**Liu Q, Kidd PB, Dobosiewicz M, Bargmann CI. *C. elegans* AWA olfactory neurons
fire calcium-mediated all-or-none action potentials.** *Cell* 175:57–70 (2018).
doi:[10.1016/j.cell.2018.08.018](https://doi.org/10.1016/j.cell.2018.08.018)
The counterexample to "worm neurons do not spike".

**Lockery SR, Goodman MB. The quest for action potentials in *C. elegans* neurons
hits a plateau.** *Nature Neuroscience* 12:377–378 (2009).
doi:[10.1038/nn0409-377](https://doi.org/10.1038/nn0409-377)
Review of graded versus plateau versus spiking signalling, including RMD plateau
potentials.

**Jiang J, Su Y, Zhang R, Li H, Tao L, Liu Q. *C. elegans* enteric motor neurons
fire synchronized action potentials underlying the defecation motor program.**
*Nature Communications* 13:2783 (2022).
doi:[10.1038/s41467-022-30452-y](https://doi.org/10.1038/s41467-022-30452-y)
AVL and DVB compound action potentials.

**Wicks SR, Roehrig CJ, Rankin CH. A dynamic network simulation of the nematode
tap withdrawal circuit: predictions concerning synaptic function using
behavioral criteria.** *Journal of Neuroscience* 16:4017–4031 (1996).
doi:[10.1523/JNEUROSCI.16-12-04017.1996](https://doi.org/10.1523/JNEUROSCI.16-12-04017.1996)
The origin of the graded leaky-integrator formulation used by most later
whole-network *C. elegans* models.

**Kunert J, Shlizerman E, Kutz JN. Low-dimensional functionality of complex
network dynamics: neurosensory integration in the *Caenorhabditis elegans*
connectome.** *Physical Review E* 89:052805 (2014).
doi:[10.1103/PhysRevE.89.052805](https://doi.org/10.1103/PhysRevE.89.052805)
The whole-connectome form of the Wicks model, and the source of the parameter set
quoted in [model_assumptions.md](model_assumptions.md) §6.2 — every constant of
which is an order-of-magnitude assumption rather than a measurement.

**Kunert-Graf JM, Shlizerman E, Walker A, Kutz JN. Multistability and
long-timescale transients encoded by network structure in a model of *C. elegans*
connectome dynamics.** *Frontiers in Computational Neuroscience* 11:53 (2017).
doi:[10.3389/fncom.2017.00053](https://doi.org/10.3389/fncom.2017.00053)

**Nicoletti M, Loppini A, Chiodo L, Folli V, Ruocco G, Filippi S. Biophysical
modeling of *C. elegans* neurons: single ion currents and whole-cell dynamics of
AWC^on and RMD.** *PLOS ONE* 14:e0218738 (2019).
doi:[10.1371/journal.pone.0218738](https://doi.org/10.1371/journal.pone.0218738)
Conductance-based models — the most faithful option, and characterized for only
two neuron types.

---

## Behaviour and circuits

**Chalfie M, Sulston JE, White JG, Southgate E, Thomson JN, Brenner S. The neural
circuit for touch sensitivity in *Caenorhabditis elegans*.** *Journal of
Neuroscience* 5:956–964 (1985).
doi:[10.1523/JNEUROSCI.05-04-00956.1985](https://doi.org/10.1523/JNEUROSCI.05-04-00956.1985)

**Kaplan JM, Horvitz HR. A dual mechanosensory and chemosensory neuron in
*Caenorhabditis elegans*.** *PNAS* 90:2227–2231 (1993).
doi:[10.1073/pnas.90.6.2227](https://doi.org/10.1073/pnas.90.6.2227)
ASH as a polymodal nociceptor.

**Bargmann CI, Hartwieg E, Horvitz HR. Odorant-selective genes and neurons
mediate olfaction in *C. elegans*.** *Cell* 74:515–527 (1993).
doi:[10.1016/0092-8674(93)80053-H](https://doi.org/10.1016/0092-8674(93)80053-H)

**Pierce-Shimomura JT, Morse TM, Lockery SR. The fundamental role of pirouettes in
*Caenorhabditis elegans* chemotaxis.** *Journal of Neuroscience* 19:9557–9569 (1999).
doi:[10.1523/JNEUROSCI.19-21-09557.1999](https://doi.org/10.1523/JNEUROSCI.19-21-09557.1999)
Chemotaxis as a biased random walk rather than steering.

**Suzuki H, Thiele TR, Faumont S, Ezcurra M, Lockery SR, Schafer WR. Functional
asymmetry in *Caenorhabditis elegans* taste neurons and its computational role in
chemotaxis.** *Nature* 454:114–117 (2008).
doi:[10.1038/nature06927](https://doi.org/10.1038/nature06927)
ASEL/ASER asymmetry.

**Macosko EZ, Pokala N, Feinberg EH, Chalasani SH, Butcher RA, Clardy J,
Bargmann CI. A hub-and-spoke circuit drives pheromone attraction and social
behaviour in *C. elegans*.** *Nature* 458:1171–1175 (2009).
doi:[10.1038/nature07886](https://doi.org/10.1038/nature07886)
The RMG gap-junction hub.

**Kawano T, Po MD, Gao S, Leung G, Ryu WS, Zhen M. An imbalancing act: gap
junctions reduce the backward motor circuit activity to bias *C. elegans* for
forward locomotion.** *Neuron* 72:572–586 (2011).
doi:[10.1016/j.neuron.2011.09.005](https://doi.org/10.1016/j.neuron.2011.09.005)

---

## Alternative modelling approaches

**Morrison M, Young L-S.** A data-driven biophysical network model reproduces
*C. elegans* premotor neural dynamics. arXiv:[2501.00278](https://arxiv.org/abs/2501.00278) (2024).
A phenomenological model whose state variable is calcium-imaging brightness
rather than voltage, with weights fitted to imaging data — discussed as option B
in [model_assumptions.md](model_assumptions.md) §6.2.

---

## Tools and resources

- **OpenWorm *C. elegans* Connectome Toolbox** — <https://openworm.org/ConnectomeToolbox/>, [source](https://github.com/openworm/ConnectomeToolbox) (MIT). Used as a mirror for published data files and as an independent source of the dataset totals our tests verify against. **Not a runtime dependency.**
- **WormAtlas** — <https://www.wormatlas.org/>
- **WormBase** — <https://wormbase.org/>
- **NemaNode** — <http://nemanode.org>
- **CeNGEN** (single-cell RNA-seq) — <https://www.cengen.org/>
- **BossDB Witvliet 2020 EM volumes** — <https://bossdb.org/project/witvliet2020>

## For Part II (*Drosophila*), not yet used

- **FlyWire / Codex** — <https://codex.flywire.ai/>
- Dorkenwald S, et al. **Neuronal wiring diagram of an adult brain.** *Nature* 634:124–138 (2024). doi:[10.1038/s41586-024-07558-y](https://doi.org/10.1038/s41586-024-07558-y)
- Schlegel P, et al. **Whole-brain annotation and multi-connectome cell typing of *Drosophila*.** *Nature* 634:139–152 (2024). doi:[10.1038/s41586-024-07686-5](https://doi.org/10.1038/s41586-024-07686-5)
