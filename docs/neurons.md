# Neuron reference

Every cell in **`witvliet_2021_7`** (Caenorhabditis elegans, hermaphrodite, adult), with what its name stands for, what it does, and how heavily it is wired.

> **This file is generated.** Run `python tools/generate_neuron_reference.py` to rebuild it. Edit the generator, not this file.

## How to read a name

*C. elegans* neuron names are acronyms, and the same cell has the same name in every animal — the worm's cell lineage is invariant, so `AVAL` is one specific cell, not a cell type. `ADEL` is the **A**nterior **DE**irid neuron, **L**eft; `ADER` is its right-hand partner. Both belong to class `ADE`.

The trailing letters are usually position:

| Suffix | Meaning |
| --- | --- |
| `L` / `R` | left / right of the bilateral pair |
| `D` / `V` | dorsal / ventral |
| `DL`, `DR`, `VL`, `VR` | one of four radial quadrants |
| a number | position along the body, head to tail (`VB1` … `VB11`) |

The **lineage** column is the cell-division path from the fertilized egg. It is the reason named cells exist at all: division is invariant, so the same cell appears in the same place in every hermaphrodite.

## Notable cells

Hand-written notes, each with its own citation, for the cells you are most likely to meet first — the outliers in the figures and the ones the planned experiments depend on. Everything below this section is generated from data.

### ADE

Dopaminergic mechanosensory neuron of the anterior deirid, involved in sensing bacterial lawn texture and the resulting slowing response.

- **In `witvliet_2021_7`:** `ADEL`, `ADER`
- **Source:** Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000).

### AFD

Principal thermosensory neuron, setting the temperature the animal prefers based on prior experience.

- **In `witvliet_2021_7`:** `AFDL`, `AFDR`
- **Source:** Mori I, Ohshima Y. Nature 376:344-348 (1995).

### AIA

First-layer amphid interneuron receiving convergent chemosensory input.

- **In `witvliet_2021_7`:** `AIAL`, `AIAR`
- **Source:** White JG, et al. Phil Trans R Soc Lond B 314:1-340 (1986).

### AIB

First-layer amphid interneuron acting largely in opposition to AIY in the chemotaxis turning decision.

- **In `witvliet_2021_7`:** `AIBL`, `AIBR`
- **Source:** Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).

### AIY

First-layer amphid interneuron; a major integration point downstream of chemosensory and thermosensory input, implicated in the turn/run decision underlying chemotaxis.

- **In `witvliet_2021_7`:** `AIYL`, `AIYR`
- **Source:** Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).

### AIZ

First-layer amphid interneuron in the chemotaxis pathway.

- **In `witvliet_2021_7`:** `AIZL`, `AIZR`
- **Source:** Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).

### ALA

Interneuron mediating a sleep-like quiescent state after cellular stress.

- **In `witvliet_2021_7`:** `ALA`
- **Source:** Hill AJ, et al. Curr Biol 24:2399-2405 (2014).

### ALM

Anterior gentle-touch receptor. One of six touch receptor neurons; ablating them abolishes the response to gentle body touch. Anterior touch is associated with reversal.

- **In `witvliet_2021_7`:** `ALML`, `ALMR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### AQR

Oxygen-sensing neuron with a process exposed to the body cavity, part of the RMG hub circuit.

- **In `witvliet_2021_7`:** `AQR`
- **Source:** Macosko EZ, et al. Nature 458:1171-1175 (2009).

### ASE

Principal salt-sensing chemosensory pair, and functionally left/right asymmetric: ASEL responds to salt increases and ASER to decreases. The two members of one class are not interchangeable.

- **In `witvliet_2021_7`:** `ASEL`, `ASER`
- **Source:** Suzuki H, et al. Nature 454:114-117 (2008).

### ASH

Polymodal nociceptor: responds to noxious chemicals, high osmolarity and nose touch, and drives avoidance. A good example of why a single sensory/inter/motor label is a simplification.

- **In `witvliet_2021_7`:** `ASHL`, `ASHR`
- **Source:** Kaplan JM, Horvitz HR. PNAS 90:2227-2231 (1993).

### AVA

Premotor command interneuron associated with backward locomotion. Among the most heavily connected cells in the animal.

- **In `witvliet_2021_7`:** `AVAL`, `AVAR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985); Kawano T, et al. Neuron 72:572-586 (2011).

### AVB

Premotor command interneuron associated with forward locomotion.

- **In `witvliet_2021_7`:** `AVBL`, `AVBR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985); Kawano T, et al. Neuron 72:572-586 (2011).

### AVD

Premotor interneuron in the backward-locomotion pathway, downstream of the anterior touch receptors.

- **In `witvliet_2021_7`:** `AVDL`, `AVDR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### AVE

Premotor interneuron associated with backward locomotion, acting on the anterior body.

- **In `witvliet_2021_7`:** `AVEL`, `AVER`
- **Source:** Kawano T, et al. Neuron 72:572-586 (2011).

### AVM

Anterior ventral gentle-touch receptor, born post-embryonically. Works with ALM for anterior touch.

- **In `witvliet_2021_7`:** `AVM`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### AWA

Olfactory neuron for attractive volatile odours. Notable for the model: AWA fires genuine calcium-mediated all-or-none action potentials, which the common 'C. elegans neurons are graded, not spiking' summary omits.

- **In `witvliet_2021_7`:** `AWAL`, `AWAR`
- **Source:** Bargmann CI, et al. Cell 74:515-527 (1993); Liu Q, et al. Cell 175:57-70 (2018).

### AWB

Olfactory neuron mediating avoidance of repulsive volatile odours.

- **In `witvliet_2021_7`:** `AWBL`, `AWBR`
- **Source:** Troemel ER, et al. Cell 91:161-169 (1997).

### AWC

Olfactory neuron for attractive volatile odours, and asymmetric between left and right in the odorants it detects.

- **In `witvliet_2021_7`:** `AWCL`, `AWCR`
- **Source:** Bargmann CI, et al. Cell 74:515-527 (1993); Wes PD, Bargmann CI. Nature 410:698-701 (2001).

### BAG

Sensory neuron responding to carbon dioxide and to falling oxygen.

- **In `witvliet_2021_7`:** `BAGL`, `BAGR`
- **Source:** Hallem EA, Sternberg PW. PNAS 105:8038-8043 (2008).

### CAN

Makes no chemical synapses at all, yet the animal dies without it. Its function remains unclear, which is why this repository records its role as unknown rather than guessing.

- **In `witvliet_2021_7`:** `CANR`
- **Source:** Forrester WC, Garriga G. Development 124:1831-1843 (1997).

### CEP

Dopaminergic mechanosensory neuron of the head, acting with ADE and PDE in the food-induced slowing response.

- **In `witvliet_2021_7`:** `CEPDL`, `CEPDR`, `CEPVL`, `CEPVR`
- **Source:** Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000).

### DVA

Interneuron with a stretch-sensitive, proprioceptive role: it reports the body's own bending back to the nervous system, so it is both an interneuron by anatomy and a sensor by function.

- **In `witvliet_2021_7`:** `DVA`
- **Source:** Li W, Feng Z, Sternberg PW, Xu XZS. Nature 440:684-687 (2006).

### HSN

Hermaphrodite-specific serotonergic motor neuron driving egg laying.

- **In `witvliet_2021_7`:** `HSNL`, `HSNR`
- **Source:** Desai C, Garriga G, McIntire SL, Horvitz HR. Nature 336:638-646 (1988).

### NSM

Serotonergic neurosecretory-motor neuron of the pharynx; detects bacterial food and slows the animal. Combines sensory, motor and neurosecretory function, so a single SIM role does not apply.

- **In `witvliet_2021_7`:** _not in this dataset_
- **Source:** Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000); Rhoades JL, et al. Cell 176:85-97 (2019).

### PLM

Posterior gentle-touch receptor. Posterior touch is associated with forward acceleration rather than reversal.

- **In `witvliet_2021_7`:** _not in this dataset_
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### PVC

Premotor command interneuron associated with forward locomotion, downstream of the posterior touch receptors.

- **In `witvliet_2021_7`:** `PVCL`, `PVCR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### PVD

Highly branched nociceptor responding to harsh touch and cold.

- **In `witvliet_2021_7`:** _not in this dataset_
- **Source:** Way JC, Chalfie M. Genes Dev 3:1823-1833 (1989); Chatzigeorgiou M, et al. Nat Neurosci 13:861-868 (2010).

### RIA

Highly connected ring interneuron integrating sensory input with head motor output; compartmentalized calcium dynamics within one cell.

- **In `witvliet_2021_7`:** `RIAL`, `RIAR`
- **Source:** Hendricks M, et al. Nature 487:99-103 (2012).

### RIM

Motor/interneuron coupled to the backward-locomotion circuit; its published label records disagreement between studies about whether it is an interneuron or a motor neuron.

- **In `witvliet_2021_7`:** `RIML`, `RIMR`
- **Source:** Kawano T, et al. Neuron 72:572-586 (2011).

### RIP

The only direct connection between the pharyngeal nervous system and the rest of the animal.

- **In `witvliet_2021_7`:** `RIPL`, `RIPR`
- **Source:** Albertson DG, Thomson JN. Phil Trans R Soc Lond B 275:299-325 (1976).

### RIS

Interneuron that induces developmentally timed sleep.

- **In `witvliet_2021_7`:** `RIS`
- **Source:** Turek M, Lewandrowski I, Bringmann H. Curr Biol 23:2215-2223 (2013).

### RMG

Hub of a gap-junction 'hub-and-spoke' circuit: several sensory neurons are electrically coupled to RMG, which aggregates them.

- **In `witvliet_2021_7`:** `RMGL`, `RMGR`
- **Source:** Macosko EZ, et al. Nature 458:1171-1175 (2009).

### URX

Oxygen-sensing neuron, electrically coupled into the RMG hub circuit.

- **In `witvliet_2021_7`:** `URXL`, `URXR`
- **Source:** Macosko EZ, et al. Nature 458:1171-1175 (2009).

> Not present in `witvliet_2021_7`: `NSM`, `PLM`, `PVD`. This is a head reconstruction — see [datasets.md](datasets.md).

## Most heavily connected neurons

Ranked by total number of partners in this dataset. Wiring counts are **measured**; role and neurotransmitter are annotations from other work.

| Neuron | Class | Role | NT | Chem in | Chem out | Gap partners | Description |
| --- | --- | --- | --- | ---: | ---: | ---: | --- |
| `RIML` | RIM | inter/motor | Glu | 24 | 27 | 9 | Ring motor neuron |
| `AIBR` | AIB | inter | Glu | 33 | 21 | 5 | Amphid interneuron |
| `RIAL` | RIA | inter | Glu | 32 | 21 | 2 | Ring interneuron, many synapses |
| `RIH` | RIH | inter | ACh | 14 | 36 | 4 | Ring interneuron |
| `RMDR` | RMD | motor | ACh | 26 | 20 | 6 | Ring motor neuron/interneuron, many synapses |
| `RIAR` | RIA | inter | Glu | 27 | 22 | 1 | Ring interneuron, many synapses |
| `RIMR` | RIM | inter/motor | Glu | 21 | 23 | 6 | Ring motor neuron |
| `AIBL` | AIB | inter | Glu | 30 | 15 | 1 | Amphid interneuron |
| `RIBL` | RIB | inter | GABA | 17 | 19 | 10 | Ring interneuron |
| `RIGR` | RIG | inter | Glu | 18 | 17 | 10 | Ring interneuron |
| `RMDDL` | RMD | motor | ACh | 21 | 16 | 8 | Ring motor neuron/interneuron, many synapses |
| `RMDL` | RMD | motor | ACh | 22 | 19 | 4 | Ring motor neuron/interneuron, many synapses |
| `AVBL` | AVB | inter | ACh | 40 | 0 | 4 | Ventral cord interneuron |
| `AVBR` | AVB | inter | ACh | 38 | 0 | 6 | Ventral cord interneuron |
| `DVA` | DVA | sensory | ACh | 15 | 28 | 1 | Ring interneurons, cell bodies in dorsorectal ganglion |
| `RICR` | RIC | inter | OA | 21 | 19 | 4 | Ring interneuron |
| `RMGR` | RMG | inter | — | 12 | 27 | 5 | Ring interneuron |
| `AVEL` | AVE | inter | ACh | 32 | 0 | 11 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord |
| `AVER` | AVE | inter | ACh | 35 | 0 | 8 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord |
| `RIBR` | RIB | inter | GABA | 17 | 15 | 11 | Ring interneuron |
| `RMGL` | RMG | inter | — | 9 | 29 | 5 | Ring interneuron |
| `ADER` | ADE | sensory | DA | 6 | 34 | 2 | Anterior deirid, sensory neuron |
| `RIGL` | RIG | inter | Glu | 8 | 22 | 10 | Ring interneuron |
| `RMDDR` | RMD | motor | ACh | 21 | 12 | 7 | Ring motor neuron/interneuron, many synapses |
| `URXL` | URX | sensory | ACh | 15 | 18 | 7 | Ring interneuron |

## All neurons in this dataset

181 neurons. Sorted by name.

**Where the columns disagree, that disagreement is real and is left visible.** `DVA` is listed with role *sensory* but described as a ring interneuron: it is an interneuron by anatomy that acts as a stretch receptor, reporting the body's own bending. `RIM` carries both *inter* and *motor* because the published label records a difference between studies. `RMG` is a major hub with no known transmitter at all. Flattening any of these would hide something true.

Neurotransmitter marks: `*` dim and variable reporter expression; 
`†` taken up from neighbouring cells rather than synthesized; 
`—` no known transmitter pathway gene detected (a published result, not a gap).

| Neuron | Class | Name stands for | Role | NT | Partners | Description | Lineage |
| --- | --- | --- | --- | --- | ---: | --- | --- |
| `ADAL` | ADA | Anterior Process from Deirid Commissure A Left | inter | Glu | 35 | Ring interneuron | `AB plapaaaapp` |
| `ADAR` | ADA | Anterior Process from Deirid Commissure A Right | inter | Glu | 34 | Ring interneuron | `AB prapaaaapp` |
| `ADEL` | ADE | Anterior DEirid Neuron Left | sensory | DA | 38 | Anterior deirid, sensory neuron | `AB plapaaaapa` |
| `ADER` | ADE | Anterior DEirid Neuron Right | sensory | DA | 42 | Anterior deirid, sensory neuron | `AB prapaaaapa` |
| `ADFL` | ADF | Amphid Dual Ciliated Ending F Left | sensory | ACh | 27 | Amphid neuron | `AB alpppppaa` |
| `ADFR` | ADF | Amphid Dual Ciliated Ending F Right | sensory | ACh | 21 | Amphid neuron | `AB praaappaa` |
| `ADLL` | ADL | Amphid Dual Ciliated Ending L Left | sensory | Glu | 29 | Amphid neuron | `AB alppppaad` |
| `ADLR` | ADL | Amphid Dual Ciliated Ending L Right | sensory | Glu | 18 | Amphid neuron | `AB praaapaad` |
| `AFDL` | AFD | Amphid Finger-like Endings D Left | sensory | Glu | 12 | Amphid finger cell | `AB alpppapav` |
| `AFDR` | AFD | Amphid Finger-like Endings D Right | sensory | Glu | 9 | Amphid finger cell | `AB praaaapav` |
| `AIAL` | AIA | Anterior Interneuron A Left | inter | ACh | 35 | Amphid interneuron | `AB plppaappa` |
| `AIAR` | AIA | Anterior Interneuron A Right | inter | ACh | 35 | Amphid interneuron | `AB prppaappa` |
| `AIBL` | AIB | Anterior Interneuron B Left | inter | Glu | 46 | Amphid interneuron | `AB plaapappa` |
| `AIBR` | AIB | Anterior Interneuron B Right | inter | Glu | 59 | Amphid interneuron | `AB praapappa` |
| `AIML` | AIM | Anterior Interneuron M Left | inter | Glu | 27 | Ring interneuron | `AB plpaapppa` |
| `AIMR` | AIM | Anterior Interneuron M Right | inter | Glu | 29 | Ring interneuron | `AB prpaapppa` |
| `AINL` | AIN | Anterior Interneuron N Left | inter | ACh | 19 | Ring interneuron | `AB alaaaalal` |
| `AINR` | AIN | Anterior Interneuron N Right | inter | ACh | 17 | Ring interneuron | `AB alaapaaar` |
| `AIYL` | AIY | Anterior Interneuron Y Left | inter | ACh | 22 | Amphid interneuron | `AB plpapaaap` |
| `AIYR` | AIY | Anterior Interneuron Y Right | inter | ACh | 24 | Amphid interneuron | `AB prpapaaap` |
| `AIZL` | AIZ | Anterior Interneuron Z Left | inter | Glu | 32 | Amphid interneuron | `AB plapaaapav` |
| `AIZR` | AIZ | Anterior Interneuron Z Right | inter | Glu | 38 | Amphid interneuron | `AB prapaaapav` |
| `ALA` | ALA | Anterior Lateral Neuron A | inter | GABA | 8 | Neuron, sends processes laterally and along dorsal cord | `AB alapppaaa` |
| `ALML` | ALM | Anterior Lateral Microtubule Neuron Left | sensory | Glu | 15 | Anterior lateral microtubule cell | `AB arppaappa` |
| `ALMR` | ALM | Anterior Lateral Microtubule Neuron Right | sensory | Glu | 19 | Anterior lateral microtubule cell | `AB arpppappa` |
| `ALNL` | ALN | Anterior Lateral Neuron N Left | sensory | ACh | 9 | Neuron associated with ALM | `AB plapappppap` |
| `ALNR` | ALN | Anterior Lateral Neuron N Right | sensory | ACh | 13 | Neuron associated with ALM | `AB prapappppap` |
| `AQR` | AQR | Anterior, Q-cell Derived Receptor | sensory | Glu | 33 | Neuron, basal body. Not part of a sensillum, projects into ring | `QR.ap` |
| `ASEL` | ASE | Amphid Single Cilium E Left | sensory | Glu | 21 | Amphid neurons, single ciliated endings | `AB alppppppaa` |
| `ASER` | ASE | Amphid Single Cilium E Right | sensory | Glu | 33 | Amphid neurons, single ciliated endings | `AB praaapppaa` |
| `ASGL` | ASG | Amphid Single Cilium G Left | sensory | Glu | 9 | Amphid neurons, single ciliated endings | `AB plaapapap` |
| `ASGR` | ASG | Amphid Single Cilium G Right | sensory | Glu | 11 | Amphid neurons, single ciliated endings | `AB praapapap` |
| `ASHL` | ASH | Amphid Single Cilium H Left | sensory | Glu | 31 | Amphid neurons, single ciliated endings | `AB plpaappaa` |
| `ASHR` | ASH | Amphid Single Cilium H Right | sensory | Glu | 30 | Amphid neurons, single ciliated endings | `AB prpaappaa` |
| `ASIL` | ASI | Amphid Single Cilium I Left | sensory | betaine † | 16 | Amphid neurons, single ciliated endings | `AB plaapapppa` |
| `ASIR` | ASI | Amphid Single Cilium I Right | sensory | betaine † | 13 | Amphid neurons, single ciliated endings | `AB praapapppa` |
| `ASJL` | ASJ | Amphid Single Cilium J Left | sensory | ACh | 13 | Amphid neurons, single ciliated endings | `AB alpppppppa` |
| `ASJR` | ASJ | Amphid Single Cilium J Right | sensory | ACh | 7 | Amphid neurons, single ciliated endings | `AB praaappppa` |
| `ASKL` | ASK | Amphid Single Cilium K Left | sensory | Glu | 17 | Amphid neurons, single ciliated endings | `AB alpppapppa` |
| `ASKR` | ASK | Amphid Single Cilium K Right | sensory | Glu | 18 | Amphid neurons, single ciliated endings | `AB praaaapppa` |
| `AUAL` | AUA | Amphid-associated Unknown Receptor A Left | inter | Glu | 17 | Neuron, process runs with amphid processes but lacks ciliated ending | `AB alpppppppp` |
| `AUAR` | AUA | Amphid-associated Unknown Receptor A Right | inter | Glu | 18 | Neuron, process runs with amphid processes but lacks ciliated ending | `AB praaappppp` |
| `AVAL` | AVA | Anterior Ventral Process A Left | inter | ACh | 32 | Ventral cord interneuron | `AB alppaaapa` |
| `AVAR` | AVA | Anterior Ventral Process A Right | inter | ACh | 33 | Ventral cord interneuron | `AB alaappapa` |
| `AVBL` | AVB | Anterior Ventral Process B Left | inter | ACh | 44 | Ventral cord interneuron | `AB plpaapaap` |
| `AVBR` | AVB | Anterior Ventral Process B Right | inter | ACh | 44 | Ventral cord interneuron | `AB prpaapaap` |
| `AVDL` | AVD | Anterior Ventral Process D Left | inter | ACh | 21 | Ventral cord interneuron | `AB alaaapalr` |
| `AVDR` | AVD | Anterior Ventral Process D Right | inter | ACh | 24 | Ventral cord interneuron | `AB alaaapprl` |
| `AVEL` | AVE | Anterior Ventral Process E Left | inter | ACh | 43 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord | `AB alpppaaaa` |
| `AVER` | AVE | Anterior Ventral Process E Right | inter | ACh | 43 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord | `AB praaaaaaa` |
| `AVFL` | AVF | Anterior Ventral Process F Left | inter | GABA † | 18 | Interneuron | `P1.aaaa/ W.aaa` |
| `AVFR` | AVF | Anterior Ventral Process F Right | inter | GABA † | 25 | Interneuron | `P1.aaaa/ W.aaa` |
| `AVHL` | AVH | Anterior Ventral Process H Left | inter | — | 23 | Neuron, mainly postsynaptic in ventral cord and presynaptic in the ring | `AB alapaaaaa` |
| `AVHR` | AVH | Anterior Ventral Process H Right | inter | — | 30 | Neuron, mainly postsynaptic in ventral cord and presynaptic in the ring | `AB alappapaa` |
| `AVJL` | AVJ | Anterior Ventral Process J Left | inter | GABA * | 25 | Neuron, synapses like AVHL/R | `AB alapapppa` |
| `AVJR` | AVJ | Anterior Ventral Process J Right | inter | GABA * | 21 | Neuron, synapses like AVHL/R | `AB alapppppa` |
| `AVKL` | AVK | Anterior Ventral Process K Left | inter | — | 35 | Ring and ventral cord interneuron | `AB plpapapap` |
| `AVKR` | AVK | Anterior Ventral Process K Right | inter | — | 35 | Ring and ventral cord interneuron | `AB prpapapap` |
| `AVL` | AVL | Anterior Ventral Process L | inter | GABA | 8 | Ring and ventral cord interneuron and an excitatory GABAergic motor neuron for rectal muscles. Few synapses | `AB prpappaap` |
| `AVM` | AVM | Anterior Ventral Microtubule Neuron | sensory | Glu | 21 | Anterior ventral microtubule cell, touch receptor | `QR.paa` |
| `AWAL` | AWA | Amphid Wing Neuron A Left | sensory | ACh * | 21 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB plaapapaa` |
| `AWAR` | AWA | Amphid Wing Neuron A Right | sensory | ACh * | 20 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB praapapaa` |
| `AWBL` | AWB | Amphid Wing Neuron B Left | sensory | ACh | 23 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB alpppppap` |
| `AWBR` | AWB | Amphid Wing Neuron B Right | sensory | ACh | 18 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB praaappap` |
| `AWCL` | AWC | Amphid Wing Neuron C Left | sensory | Glu | 21 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB plpaaaaap` |
| `AWCR` | AWC | Amphid Wing Neuron C Right | sensory | Glu | 23 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB prpaaaaap` |
| `BAGL` | BAG | BAG-like Dendritic Ending Left | sensory | Glu | 20 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB alppappap` |
| `BAGR` | BAG | BAG-like Dendritic Ending Right | sensory | Glu | 31 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB arappppap` |
| `BDUL` | BDU | Black Vesicles, Deirid Locale, Unknown Function Left | inter | — | 22 | Neuron, process runs along excretory canal and into ring, unique darkly staining synaptic vesicles | `AB arppaappp` |
| `BDUR` | BDU | Black Vesicles, Deirid Locale, Unknown Function Right | inter | — | 17 | Neuron, process runs along excretory canal and into ring, unique darkly staining synaptic vesicles | `AB arpppappp` |
| `CANR` | CAN | Excretory CANal-associated Neuron Right | unknown | betaine † | 1 | Process runs along excretory canal, no synapses, essential for survival | `AB alappappa` |
| `CEPDL` | CEP | CEPhalic Sensory Neuron Dorsal Left | sensory | DA | 39 | Cephalic neurons, contain dopamine | `AB plaaaaappa` |
| `CEPDR` | CEP | CEPhalic Sensory Neuron Dorsal Right | sensory | DA | 36 | Cephalic neurons, contain dopamine | `AB arpapaappa` |
| `CEPVL` | CEP | CEPhalic Sensory Neuron Ventral Left | sensory | DA | 34 | Cephalic neurons, contain dopamine | `AB plpaappppa` |
| `CEPVR` | CEP | CEPhalic Sensory Neuron Ventral Right | sensory | DA | 32 | Cephalic neurons, contain dopamine | `AB prpaappppa` |
| `DVA` | DVA | Dorsorectal Ganglion Ventral Process A | sensory | ACh | 44 | Ring interneurons, cell bodies in dorsorectal ganglion | `AB prppppapp` |
| `DVC` | DVC | Dorsorectal Ganglion Ventral Process C | inter | Glu | 24 | Ring interneurons, cell bodies in dorsorectal ganglion | `C aapaa` |
| `FLPL` | FLP | FLaP-like Dendritic Ending Left | sensory | Glu | 19 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB plapaaapad` |
| `FLPR` | FLP | FLaP-like Dendritic Ending Right | sensory | Glu | 14 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB prapaaapad` |
| `HSNL` | HSN | Hermaphrodite-Specific Neuron Left | motor | ACh | 32 | Hermaphrodite specific motor neurons (die in male embryo), innervate vulval muscles, serotonergic | `AB plapppappa` |
| `HSNR` | HSN | Hermaphrodite-Specific Neuron Right | motor | ACh | 34 | Hermaphrodite specific motor neurons (die in male embryo), innervate vulval muscles, serotonergic | `AB prapppappa` |
| `IL1DL` | IL1 | Inner Labial 1 Dorsal Left | sensory | Glu | 22 | Inner labial neuron | `AB alapappaaa` |
| `IL1DR` | IL1 | Inner Labial 1 Dorsal Right | sensory | Glu | 19 | Inner labial neuron | `AB alappppaaa` |
| `IL1L` | IL1 | Inner Labial 1 Left | sensory | Glu | 30 | Inner labial neuron | `AB alapaappaa` |
| `IL1R` | IL1 | Inner Labial 1 Right | sensory | Glu | 31 | Inner labial neuron | `AB alaappppaa` |
| `IL1VL` | IL1 | Inner Labial 1 Ventral Left | sensory | Glu | 20 | Inner labial neuron | `AB alppapppaa` |
| `IL1VR` | IL1 | Inner Labial 1 Ventral Right | sensory | Glu | 19 | Inner labial neuron | `AB arapppppaa` |
| `IL2DL` | IL2 | Inner Labial 2 Dorsal Left | sensory | ACh | 18 | Inner labial neuron | `AB alapappap` |
| `IL2DR` | IL2 | Inner Labial 2 Dorsal Right | sensory | ACh | 14 | Inner labial neuron | `AB alappppap` |
| `IL2L` | IL2 | Inner Labial 2 Left | sensory | ACh | 22 | Inner labial neuron | `AB alapaappp` |
| `IL2R` | IL2 | Inner Labial 2 Right | sensory | ACh | 31 | Inner labial neuron | `AB alaappppp` |
| `IL2VL` | IL2 | Inner Labial 2 Ventral Left | sensory | ACh | 15 | Inner labial neuron | `AB alppapppp` |
| `IL2VR` | IL2 | Inner Labial 2 Ventral Right | sensory | ACh | 17 | Inner labial neuron | `AB arapppppp` |
| `OLLL` | OLL | Outer Labial Lateral Dendrite Left | sensory | Glu | 34 | Lateral outer labial neurons | `AB alppppapaa` |
| `OLLR` | OLL | Outer Labial Lateral Dendrite Right | sensory | Glu | 31 | Lateral outer labial neurons | `AB praaapapaa` |
| `OLQDL` | OLQ | Outer Labial Quadrant Dendrite Dorsal Left | sensory | Glu | 21 | Quadrant outer labial neuron | `AB alapapapaa` |
| `OLQDR` | OLQ | Outer Labial Quadrant Dendrite Dorsal Right | sensory | Glu | 23 | Quadrant outer labial neuron | `AB alapppapaa` |
| `OLQVL` | OLQ | Outer Labial Quadrant Dendrite Ventral Left | sensory | Glu | 19 | Quadrant outer labial neuron | `AB plpaaappaa` |
| `OLQVR` | OLQ | Outer Labial Quadrant Dendrite Ventral Right | sensory | Glu | 21 | Quadrant outer labial neuron | `AB prpaaappaa` |
| `PLNL` | PLN | Posterior Lateral N Left | sensory | ACh | 9 | Interneuron, associated with PLM | `TL.pppap` |
| `PLNR` | PLN | Posterior Lateral N Right | sensory | ACh | 8 | Interneuron, associated with PLM | `TR.pppap` |
| `PVCL` | PVC | Posterior Ventral Process C Left | inter | ACh | 24 | Ventral cord interneuron, cell body in lumbar ganglion, synapses onto VB and DB motor neurons, formerly called delta | `AB plpppaapaa` |
| `PVCR` | PVC | Posterior Ventral Process C Right | inter | ACh | 26 | Ventral cord interneuron, cell body in lumbar ganglion, synapses onto VB and DB motor neurons, formerly called delta | `AB prpppaapaa` |
| `PVNL` | PVN | Posterior Ventral Process N Left | inter | ACh | 28 | Interneuron/motor neuron, post. vent. cord, few synapses | `TL.appp` |
| `PVNR` | PVN | Posterior Ventral Process N Right | inter | ACh | 20 | Interneuron/motor neuron, post. vent. cord, few synapses | `TR.appp` |
| `PVPL` | PVP | Posterior Ventral Process P Left | inter | ACh | 22 | Interneuron, cell body in preanal ganglion, projects along ventral cord to nerve ring | `AB plppppaaa` |
| `PVPR` | PVP | Posterior Ventral Process P Right | inter | ACh | 21 | Interneuron, cell body in preanal ganglion, projects along ventral cord to nerve ring | `AB prppppaaa` |
| `PVQL` | PVQ | Posterior Ventral Process Q Left | inter | — | 15 | Interneuron, projects along ventral cord to ring | `AB plapppaaa` |
| `PVQR` | PVQ | Posterior Ventral Process Q Right | inter | — | 15 | Interneuron, projects along ventral cord to ring | `AB prapppaaa` |
| `PVR` | PVR | Posterior Ventral Process R | inter | Glu | 21 | Interneuron, projects along ventral cord to ring | `C aappa` |
| `PVT` | PVT | Posterior Ventral Process T | inter | — | 23 | Interneuron, projects along ventral cord to ring | `AB plpappppa` |
| `RIAL` | RIA | Ring Interneuron A Left | inter | Glu | 55 | Ring interneuron, many synapses | `AB alapaapaa` |
| `RIAR` | RIA | Ring Interneuron A Right | inter | Glu | 50 | Ring interneuron, many synapses | `AB alaapppaa` |
| `RIBL` | RIB | Ring Interneuron B Left | inter | GABA | 46 | Ring interneuron | `AB plpaappap` |
| `RIBR` | RIB | Ring Interneuron B Right | inter | GABA | 43 | Ring interneuron | `AB prpaappap` |
| `RICL` | RIC | Ring Interneuron C Left | inter | OA | 34 | Ring interneuron | `AB plppaaaapp` |
| `RICR` | RIC | Ring Interneuron C Right | inter | OA | 44 | Ring interneuron | `AB prppaaaapp` |
| `RID` | RID | Ring Interneuron D | inter | — | 9 | Ring interneuron, projects along dorsal cord | `AB alappaapa` |
| `RIFL` | RIF | Ring Interneuron F Left | inter | ACh | 21 | Ring interneuron | `AB plppapaaap` |
| `RIFR` | RIF | Ring Interneuron F Right | inter | ACh | 28 | Ring interneuron | `AB prppapaaap` |
| `RIGL` | RIG | Ring Interneuron G Left | inter | Glu | 40 | Ring interneuron | `AB plppappaa` |
| `RIGR` | RIG | Ring Interneuron G Right | inter | Glu | 45 | Ring interneuron | `AB prppappaa` |
| `RIH` | RIH | Ring Interneuron H | inter | ACh | 54 | Ring interneuron | `AB prpappaaa` |
| `RIML` | RIM | Ring Interneuron M Left | inter/motor | Glu | 60 | Ring motor neuron | `AB plppaapap` |
| `RIMR` | RIM | Ring Interneuron M Right | inter/motor | Glu | 50 | Ring motor neuron | `AB prppaapap` |
| `RIPL` | RIP | Ring Interneuron P Left | inter | ACh | 11 | Ring/pharynx interneuron, only direct connection between pharynx and ring | `AB alpapaaaa` |
| `RIPR` | RIP | Ring Interneuron P Right | inter | ACh | 14 | Ring/pharynx interneuron, only direct connection between pharynx and ring | `AB arappaaaa` |
| `RIR` | RIR | Ring Interneuron R | inter | ACh | 36 | Ring interneuron | `AB prpapppaa` |
| `RIS` | RIS | Ring Interneuron S | inter | GABA | 39 | Ring interneuron | `AB prpappapa` |
| `RIVL` | RIV | Ring Interneuron V Left | motor | ACh | 24 | Ring interneuron | `AB plpaapaaa` |
| `RIVR` | RIV | Ring Interneuron V Right | motor | ACh | 32 | Ring interneuron | `AB prpaapaaa` |
| `RMDDL` | RMD | Ring Motor Neuron D Dorsal Left | motor | ACh | 45 | Ring motor neuron/interneuron, many synapses | `AB alpapapaa` |
| `RMDDR` | RMD | Ring Motor Neuron D Dorsal Right | motor | ACh | 40 | Ring motor neuron/interneuron, many synapses | `AB arappapaa` |
| `RMDL` | RMD | Ring Motor Neuron D Left | motor | ACh | 45 | Ring motor neuron/interneuron, many synapses | `AB alpppapad` |
| `RMDR` | RMD | Ring Motor Neuron D Right | motor | ACh | 52 | Ring motor neuron/interneuron, many synapses | `AB praaaapad` |
| `RMDVL` | RMD | Ring Motor Neuron D Ventral Left | motor | ACh | 38 | Ring motor neuron/interneuron, many synapses | `AB alppapaaa` |
| `RMDVR` | RMD | Ring Motor Neuron D Ventral Right | motor | ACh | 39 | Ring motor neuron/interneuron, many synapses | `AB arapppaaa` |
| `RMED` | RME | Ring Motor Neuron E Dorsal | motor | GABA | 31 | Ring motor neuron | `AB alapppaap` |
| `RMEL` | RME | Ring Motor Neuron E Left | motor | GABA | 25 | Ring motor neuron | `AB alaaaarlp` |
| `RMER` | RME | Ring Motor Neuron E Right | motor | GABA | 31 | Ring motor neuron | `AB alaaaarrp` |
| `RMEV` | RME | Ring Motor Neuron E Ventral | motor | GABA | 26 | Ring motor neuron | `AB plpappaaa` |
| `RMFL` | RMF | Ring Motor Neuron F Left | inter | ACh | 23 | Ring motor neuron/interneuron | `G2.al` |
| `RMFR` | RMF | Ring Motor Neuron F Right | inter | ACh | 36 | Ring motor neuron/interneuron | `G2.ar` |
| `RMGL` | RMG | Ring Motor Neuron G Left | inter | — | 43 | Ring interneuron | `AB plapaaapp` |
| `RMGR` | RMG | Ring Motor Neuron G Right | inter | — | 44 | Ring interneuron | `AB prapaaapp` |
| `RMHL` | RMH | Ring Motor Neuron H Left | motor | ACh | 27 | Ring motor neuron/interneuron | `G1.l` |
| `RMHR` | RMH | Ring Motor Neuron H Right | motor | ACh | 33 | Ring motor neuron/interneuron | `G1.r` |
| `SAADL` | SAA | Sublateral Anterior A Dorsal Left | inter | ACh | 18 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB alppapapa` |
| `SAADR` | SAA | Sublateral Anterior A Dorsal Right | inter | ACh | 17 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB arapppapa` |
| `SAAVL` | SAA | Sublateral Anterior A Ventral Left | inter | ACh | 16 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB plpaaaaaa` |
| `SAAVR` | SAA | Sublateral Anterior A Ventral Right | inter | ACh | 17 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB prpaaaaaa` |
| `SDQL` | SDQ | Sublateral Dorsal Q-cell Derived Left | sensory | ACh | 24 | Post. lateral interneuron, process projects into ring | `QL.pap` |
| `SDQR` | SDQ | Sublateral Dorsal Q-cell Derived Right | sensory | ACh | 19 | Ant. lateral interneuron, process projects into ring | `QR.pap` |
| `SIADL` | SIA | Sublateral Interneuron A Dorsal Left | inter/motor | ACh | 9 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB plpapaapa` |
| `SIADR` | SIA | Sublateral Interneuron A Dorsal Right | inter/motor | ACh | 5 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB prpapaapa` |
| `SIAVL` | SIA | Sublateral Interneuron A Ventral Left | inter/motor | ACh | 10 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB plpapappa` |
| `SIAVR` | SIA | Sublateral Interneuron A Ventral Right | inter/motor | ACh | 7 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB prpapappa` |
| `SIBDL` | SIB | Sublateral Interneuron B Dorsal Left | inter/motor | ACh | 12 | Similar to SIA | `AB plppaaaaa` |
| `SIBDR` | SIB | Sublateral Interneuron B Dorsal Right | inter/motor | ACh | 8 | Similar to SIA | `AB prppaaaaa` |
| `SIBVL` | SIB | Sublateral Interneuron B Ventral Left | inter/motor | ACh | 9 | Similar to SIA | `AB plpapaapp` |
| `SIBVR` | SIB | Sublateral Interneuron B Ventral Right | inter/motor | ACh | 12 | Similar to SIA | `AB prpapaapp` |
| `SMBDL` | SMB | Sublateral Motor Neuron B Dorsal Left | motor | ACh | 32 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB alpapapapp` |
| `SMBDR` | SMB | Sublateral Motor Neuron B Dorsal Right | motor | ACh | 28 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB arappapapp` |
| `SMBVL` | SMB | Sublateral Motor Neuron B Ventral Left | motor | ACh | 28 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB alpapappp` |
| `SMBVR` | SMB | Sublateral Motor Neuron B Ventral Right | motor | ACh | 30 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB arappappp` |
| `SMDDL` | SMD | Sublateral Motor Neuron D Dorsal Left | motor | ACh | 34 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB plpapaaaa` |
| `SMDDR` | SMD | Sublateral Motor Neuron D Dorsal Right | motor | ACh | 31 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB prpapaaaa` |
| `SMDVL` | SMD | Sublateral Motor Neuron D Ventral Left | motor | ACh | 32 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB alppappaa` |
| `SMDVR` | SMD | Sublateral Motor Neuron D Ventral Right | motor | ACh | 33 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB arappppaa` |
| `URADL` | URA | Unknown Receptor, not Ciliated A Dorsal Left | motor | ACh | 12 | Ring motor neuron | `AB plaaaaaaa` |
| `URADR` | URA | Unknown Receptor, not Ciliated A Dorsal Right | motor | ACh | 12 | Ring motor neuron | `AB arpapaaaa` |
| `URAVL` | URA | Unknown Receptor, not Ciliated A Ventral Left | motor | ACh | 16 | Ring motor neuron | `AB plpaaapaa` |
| `URAVR` | URA | Unknown Receptor, not Ciliated A Ventral Right | motor | ACh | 16 | Ring motor neuron | `AB prpaaapaa` |
| `URBL` | URB | Unknown Receptor, not Ciliated B Left | inter | ACh | 27 | Neuron, presynaptic in ring, ending in head | `AB plaapaapa` |
| `URBR` | URB | Unknown Receptor, not Ciliated B Right | inter | ACh | 22 | Neuron, presynaptic in ring, ending in head | `AB praapaapa` |
| `URXL` | URX | Unknown Receptor, not Ciliated X Left | sensory | ACh | 40 | Ring interneuron | `AB plaaaaappp` |
| `URXR` | URX | Unknown Receptor, not Ciliated X Right | sensory | ACh | 31 | Ring interneuron | `AB arpapaappp` |
| `URYDL` | URY | Unknown Receptor, not Ciliated Y Dorsal Left | sensory | Glu | 23 | Neuron, presynaptic in ring, ending in head | `AB alapapapp` |
| `URYDR` | URY | Unknown Receptor, not Ciliated Y Dorsal Right | sensory | Glu | 18 | Neuron, presynaptic in ring, ending in head | `AB alapppapp` |
| `URYVL` | URY | Unknown Receptor, not Ciliated Y Ventral Left | sensory | Glu | 22 | Neuron, presynaptic in ring, ending in head | `AB plpaaappp` |
| `URYVR` | URY | Unknown Receptor, not Ciliated Y Ventral Right | sensory | Glu | 19 | Neuron, presynaptic in ring, ending in head | `AB prpaaappp` |

## Non-neuronal cells in this dataset

41 cells that are not neurons. Muscles are the nervous system's actual output, and glia are increasingly understood to shape neural function, so neither is filtered out.

| Cell | Category | Partners | Description |
| --- | --- | ---: | --- |
| `CEPshDL` | glia | 2 | Sheath cell other than amphid sheath and phasmid |
| `CEPshDR` | glia | 4 | Sheath cell other than amphid sheath and phasmid |
| `CEPshVL` | glia | 8 | Sheath cell other than amphid sheath and phasmid |
| `CEPshVR` | glia | 9 | Sheath cell other than amphid sheath and phasmid |
| `MDL01` | muscle | 9 | Head muscle |
| `MDL02` | muscle | 10 | Head muscle |
| `MDL03` | muscle | 11 | Head muscle |
| `MDL04` | muscle | 10 | Head muscle |
| `MDL05` | muscle | 7 | Head muscle |
| `MDL06` | muscle | 6 | Head muscle |
| `MDL07` | muscle | 6 | Head muscle |
| `MDL08` | muscle | 4 | Main body muscle |
| `MDR01` | muscle | 8 | Head muscle |
| `MDR02` | muscle | 12 | Head muscle |
| `MDR03` | muscle | 8 | Head muscle |
| `MDR04` | muscle | 11 | Head muscle |
| `MDR05` | muscle | 7 | Head muscle |
| `MDR06` | muscle | 9 | Head muscle |
| `MDR07` | muscle | 3 | Head muscle |
| `MDR08` | muscle | 5 | Main body muscle |
| `MVL01` | muscle | 8 | Head muscle |
| `MVL02` | muscle | 9 | Head muscle |
| `MVL03` | muscle | 11 | Head muscle |
| `MVL04` | muscle | 11 | Head muscle |
| `MVL05` | muscle | 9 | Head muscle |
| `MVL06` | muscle | 7 | Head muscle |
| `MVL07` | muscle | 4 | Head muscle |
| `MVL08` | muscle | 2 | Main body muscle |
| `MVR01` | muscle | 10 | Head muscle |
| `MVR02` | muscle | 12 | Head muscle |
| `MVR03` | muscle | 8 | Head muscle |
| `MVR04` | muscle | 14 | Head muscle |
| `MVR05` | muscle | 5 | Head muscle |
| `MVR06` | muscle | 6 | Head muscle |
| `MVR07` | muscle | 4 | Head muscle |
| `MVR08` | muscle | 3 | Main body muscle |
| `GLRDL` | other | 3 | GLR cell |
| `GLRDR` | other | 2 | GLR cell |
| `GLRL` | other | 5 | GLR cell |
| `GLRR` | other | 6 | GLR cell |
| `GLRVR` | other | 2 | GLR cell |

## Sources

| Column | Source | Confidence |
| --- | --- | --- |
| name expansion, lineage, description | WormAtlas cell listings | published annotation |
| role | WormAtlas + Cook et al. 2019 groupings | published annotation |
| class, neurotransmitter | Wang et al. 2024, eLife 12:RP95402 | published annotation |
| partner counts | Witvliet D, et al. Connectomes across development reveal principles of brain maturation. Nature 596:257-261 (2021). | **measured** |

Full citations and licences: [references.md](references.md). What each source does and does not establish: [model_assumptions.md](model_assumptions.md).
