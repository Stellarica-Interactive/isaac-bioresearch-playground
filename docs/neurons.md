# Neuron reference

Every cell in **`cook_2019_herm`** (Caenorhabditis elegans, hermaphrodite, adult), with what its name stands for, what it does, and how heavily it is wired.

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

- **In `cook_2019_herm`:** `ADEL`, `ADER`
- **Source:** Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000).

### AFD

Principal thermosensory neuron, setting the temperature the animal prefers based on prior experience.

- **In `cook_2019_herm`:** `AFDL`, `AFDR`
- **Source:** Mori I, Ohshima Y. Nature 376:344-348 (1995).

### AIA

First-layer amphid interneuron receiving convergent chemosensory input.

- **In `cook_2019_herm`:** `AIAL`, `AIAR`
- **Source:** White JG, et al. Phil Trans R Soc Lond B 314:1-340 (1986).

### AIB

First-layer amphid interneuron acting largely in opposition to AIY in the chemotaxis turning decision.

- **In `cook_2019_herm`:** `AIBL`, `AIBR`
- **Source:** Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).

### AIY

First-layer amphid interneuron; a major integration point downstream of chemosensory and thermosensory input, implicated in the turn/run decision underlying chemotaxis.

- **In `cook_2019_herm`:** `AIYL`, `AIYR`
- **Source:** Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).

### AIZ

First-layer amphid interneuron in the chemotaxis pathway.

- **In `cook_2019_herm`:** `AIZL`, `AIZR`
- **Source:** Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).

### ALA

Interneuron mediating a sleep-like quiescent state after cellular stress.

- **In `cook_2019_herm`:** `ALA`
- **Source:** Hill AJ, et al. Curr Biol 24:2399-2405 (2014).

### ALM

Anterior gentle-touch receptor. One of six touch receptor neurons; ablating them abolishes the response to gentle body touch. Anterior touch is associated with reversal.

- **In `cook_2019_herm`:** `ALML`, `ALMR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### AQR

Oxygen-sensing neuron with a process exposed to the body cavity, part of the RMG hub circuit.

- **In `cook_2019_herm`:** `AQR`
- **Source:** Macosko EZ, et al. Nature 458:1171-1175 (2009).

### ASE

Principal salt-sensing chemosensory pair, and functionally left/right asymmetric: ASEL responds to salt increases and ASER to decreases. The two members of one class are not interchangeable.

- **In `cook_2019_herm`:** `ASEL`, `ASER`
- **Source:** Suzuki H, et al. Nature 454:114-117 (2008).

### ASH

Polymodal nociceptor: responds to noxious chemicals, high osmolarity and nose touch, and drives avoidance. A good example of why a single sensory/inter/motor label is a simplification.

- **In `cook_2019_herm`:** `ASHL`, `ASHR`
- **Source:** Kaplan JM, Horvitz HR. PNAS 90:2227-2231 (1993).

### AVA

Premotor command interneuron associated with backward locomotion. Among the most heavily connected cells in the animal.

- **In `cook_2019_herm`:** `AVAL`, `AVAR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985); Kawano T, et al. Neuron 72:572-586 (2011).

### AVB

Premotor command interneuron associated with forward locomotion.

- **In `cook_2019_herm`:** `AVBL`, `AVBR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985); Kawano T, et al. Neuron 72:572-586 (2011).

### AVD

Premotor interneuron in the backward-locomotion pathway, downstream of the anterior touch receptors.

- **In `cook_2019_herm`:** `AVDL`, `AVDR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### AVE

Premotor interneuron associated with backward locomotion, acting on the anterior body.

- **In `cook_2019_herm`:** `AVEL`, `AVER`
- **Source:** Kawano T, et al. Neuron 72:572-586 (2011).

### AVM

Anterior ventral gentle-touch receptor, born post-embryonically. Works with ALM for anterior touch.

- **In `cook_2019_herm`:** `AVM`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### AWA

Olfactory neuron for attractive volatile odours. Notable for the model: AWA fires genuine calcium-mediated all-or-none action potentials, which the common 'C. elegans neurons are graded, not spiking' summary omits.

- **In `cook_2019_herm`:** `AWAL`, `AWAR`
- **Source:** Bargmann CI, et al. Cell 74:515-527 (1993); Liu Q, et al. Cell 175:57-70 (2018).

### AWB

Olfactory neuron mediating avoidance of repulsive volatile odours.

- **In `cook_2019_herm`:** `AWBL`, `AWBR`
- **Source:** Troemel ER, et al. Cell 91:161-169 (1997).

### AWC

Olfactory neuron for attractive volatile odours, and asymmetric between left and right in the odorants it detects.

- **In `cook_2019_herm`:** `AWCL`, `AWCR`
- **Source:** Bargmann CI, et al. Cell 74:515-527 (1993); Wes PD, Bargmann CI. Nature 410:698-701 (2001).

### BAG

Sensory neuron responding to carbon dioxide and to falling oxygen.

- **In `cook_2019_herm`:** `BAGL`, `BAGR`
- **Source:** Hallem EA, Sternberg PW. PNAS 105:8038-8043 (2008).

### CAN

Makes no chemical synapses at all, yet the animal dies without it. Its function remains unclear, which is why this repository records its role as unknown rather than guessing.

- **In `cook_2019_herm`:** `CANL`, `CANR`
- **Source:** Forrester WC, Garriga G. Development 124:1831-1843 (1997).

### CEP

Dopaminergic mechanosensory neuron of the head, acting with ADE and PDE in the food-induced slowing response.

- **In `cook_2019_herm`:** `CEPDL`, `CEPDR`, `CEPVL`, `CEPVR`
- **Source:** Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000).

### DVA

Interneuron with a stretch-sensitive, proprioceptive role: it reports the body's own bending back to the nervous system, so it is both an interneuron by anatomy and a sensor by function.

- **In `cook_2019_herm`:** `DVA`
- **Source:** Li W, Feng Z, Sternberg PW, Xu XZS. Nature 440:684-687 (2006).

### HSN

Hermaphrodite-specific serotonergic motor neuron driving egg laying.

- **In `cook_2019_herm`:** `HSNL`, `HSNR`
- **Source:** Desai C, Garriga G, McIntire SL, Horvitz HR. Nature 336:638-646 (1988).

### NSM

Serotonergic neurosecretory-motor neuron of the pharynx; detects bacterial food and slows the animal. Combines sensory, motor and neurosecretory function, so a single SIM role does not apply.

- **In `cook_2019_herm`:** `NSML`, `NSMR`
- **Source:** Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000); Rhoades JL, et al. Cell 176:85-97 (2019).

### PLM

Posterior gentle-touch receptor. Posterior touch is associated with forward acceleration rather than reversal.

- **In `cook_2019_herm`:** `PLML`, `PLMR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### PVC

Premotor command interneuron associated with forward locomotion, downstream of the posterior touch receptors.

- **In `cook_2019_herm`:** `PVCL`, `PVCR`
- **Source:** Chalfie M, et al. J Neurosci 5:956-964 (1985).

### PVD

Highly branched nociceptor responding to harsh touch and cold.

- **In `cook_2019_herm`:** `PVDL`, `PVDR`
- **Source:** Way JC, Chalfie M. Genes Dev 3:1823-1833 (1989); Chatzigeorgiou M, et al. Nat Neurosci 13:861-868 (2010).

### RIA

Highly connected ring interneuron integrating sensory input with head motor output; compartmentalized calcium dynamics within one cell.

- **In `cook_2019_herm`:** `RIAL`, `RIAR`
- **Source:** Hendricks M, et al. Nature 487:99-103 (2012).

### RIM

Motor/interneuron coupled to the backward-locomotion circuit; its published label records disagreement between studies about whether it is an interneuron or a motor neuron.

- **In `cook_2019_herm`:** `RIML`, `RIMR`
- **Source:** Kawano T, et al. Neuron 72:572-586 (2011).

### RIP

The only direct connection between the pharyngeal nervous system and the rest of the animal.

- **In `cook_2019_herm`:** `RIPL`, `RIPR`
- **Source:** Albertson DG, Thomson JN. Phil Trans R Soc Lond B 275:299-325 (1976).

### RIS

Interneuron that induces developmentally timed sleep.

- **In `cook_2019_herm`:** `RIS`
- **Source:** Turek M, Lewandrowski I, Bringmann H. Curr Biol 23:2215-2223 (2013).

### RMG

Hub of a gap-junction 'hub-and-spoke' circuit: several sensory neurons are electrically coupled to RMG, which aggregates them.

- **In `cook_2019_herm`:** `RMGL`, `RMGR`
- **Source:** Macosko EZ, et al. Nature 458:1171-1175 (2009).

### URX

Oxygen-sensing neuron, electrically coupled into the RMG hub circuit.

- **In `cook_2019_herm`:** `URXL`, `URXR`
- **Source:** Macosko EZ, et al. Nature 458:1171-1175 (2009).

## Most heavily connected neurons

Ranked by total number of partners in this dataset. Wiring counts are **measured**; role and neurotransmitter are annotations from other work.

| Neuron | Class | Role | NT | Chem in | Chem out | Gap partners | Description |
| --- | --- | --- | --- | ---: | ---: | ---: | --- |
| `AVAR` | AVA | inter | ACh | 59 | 49 | 51 | Ventral cord interneuron |
| `AVAL` | AVA | inter | ACh | 65 | 44 | 45 | Ventral cord interneuron |
| `DVA` | DVA | sensory | ACh | 30 | 46 | 24 | Ring interneurons, cell bodies in dorsorectal ganglion |
| `PVCR` | PVC | inter | ACh | 33 | 38 | 28 | Ventral cord interneuron, cell body in lumbar ganglion, synapses onto VB and DB motor neurons, formerly called delta |
| `AVBL` | AVB | inter | ACh | 45 | 24 | 26 | Ventral cord interneuron |
| `PVCL` | PVC | inter | ACh | 34 | 36 | 25 | Ventral cord interneuron, cell body in lumbar ganglion, synapses onto VB and DB motor neurons, formerly called delta |
| `AVBR` | AVB | inter | ACh | 45 | 21 | 28 | Ventral cord interneuron |
| `AVDR` | AVD | inter | ACh | 46 | 31 | 16 | Ventral cord interneuron |
| `AVEL` | AVE | inter | ACh | 46 | 24 | 15 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord |
| `AVER` | AVE | inter | ACh | 44 | 25 | 16 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord |
| `AVDL` | AVD | inter | ACh | 39 | 25 | 16 | Ventral cord interneuron |
| `SMDDR` | SMD | motor | ACh | 27 | 35 | 13 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally |
| `SMDDL` | SMD | motor | ACh | 23 | 38 | 11 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally |
| `AIBR` | AIB | inter | Glu | 32 | 26 | 13 | Amphid interneuron |
| `AIBL` | AIB | inter | Glu | 33 | 23 | 14 | Amphid interneuron |
| `HSNR` | HSN | motor | ACh | 19 | 44 | 7 | Hermaphrodite specific motor neurons (die in male embryo), innervate vulval muscles, serotonergic |
| `PVNR` | PVN | inter | ACh | 16 | 37 | 17 | Interneuron/motor neuron, post. vent. cord, few synapses |
| `PVR` | PVR | inter | Glu | 17 | 31 | 22 | Interneuron, projects along ventral cord to ring |
| `RIAR` | RIA | inter | Glu | 35 | 27 | 7 | Ring interneuron, many synapses |
| `SMDVL` | SMD | motor | ACh | 27 | 32 | 9 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally |
| `RIAL` | RIA | inter | Glu | 36 | 27 | 4 | Ring interneuron, many synapses |
| `PVT` | PVT | inter | — | 20 | 13 | 31 | Interneuron, projects along ventral cord to ring |
| `RMGL` | RMG | inter | — | 20 | 33 | 10 | Ring interneuron |
| `PVPR` | PVP | inter | ACh | 17 | 23 | 22 | Interneuron, cell body in preanal ganglion, projects along ventral cord to nerve ring |
| `RIML` | RIM | inter/motor | Glu | 26 | 27 | 9 | Ring motor neuron |

## All neurons in this dataset

302 neurons. Sorted by name.

**Where the columns disagree, that disagreement is real and is left visible.** `DVA` is listed with role *sensory* but described as a ring interneuron: it is an interneuron by anatomy that acts as a stretch receptor, reporting the body's own bending. `RIM` carries both *inter* and *motor* because the published label records a difference between studies. `RMG` is a major hub with no known transmitter at all. Flattening any of these would hide something true.

Neurotransmitter marks: `*` dim and variable reporter expression; 
`†` taken up from neighbouring cells rather than synthesized; 
`—` no known transmitter pathway gene detected (a published result, not a gap).

| Neuron | Class | Name stands for | Role | NT | Partners | Description | Lineage |
| --- | --- | --- | --- | --- | ---: | --- | --- |
| `ADAL` | ADA | Anterior Process from Deirid Commissure A Left | inter | Glu | 43 | Ring interneuron | `AB plapaaaapp` |
| `ADAR` | ADA | Anterior Process from Deirid Commissure A Right | inter | Glu | 39 | Ring interneuron | `AB prapaaaapp` |
| `ADEL` | ADE | Anterior DEirid Neuron Left | sensory | DA | 55 | Anterior deirid, sensory neuron | `AB plapaaaapa` |
| `ADER` | ADE | Anterior DEirid Neuron Right | sensory | DA | 45 | Anterior deirid, sensory neuron | `AB prapaaaapa` |
| `ADFL` | ADF | Amphid Dual Ciliated Ending F Left | sensory | ACh | 40 | Amphid neuron | `AB alpppppaa` |
| `ADFR` | ADF | Amphid Dual Ciliated Ending F Right | sensory | ACh | 49 | Amphid neuron | `AB praaappaa` |
| `ADLL` | ADL | Amphid Dual Ciliated Ending L Left | sensory | Glu | 38 | Amphid neuron | `AB alppppaad` |
| `ADLR` | ADL | Amphid Dual Ciliated Ending L Right | sensory | Glu | 46 | Amphid neuron | `AB praaapaad` |
| `AFDL` | AFD | Amphid Finger-like Endings D Left | sensory | Glu | 24 | Amphid finger cell | `AB alpppapav` |
| `AFDR` | AFD | Amphid Finger-like Endings D Right | sensory | Glu | 22 | Amphid finger cell | `AB praaaapav` |
| `AIAL` | AIA | Anterior Interneuron A Left | inter | ACh | 49 | Amphid interneuron | `AB plppaappa` |
| `AIAR` | AIA | Anterior Interneuron A Right | inter | ACh | 37 | Amphid interneuron | `AB prppaappa` |
| `AIBL` | AIB | Anterior Interneuron B Left | inter | Glu | 70 | Amphid interneuron | `AB plaapappa` |
| `AIBR` | AIB | Anterior Interneuron B Right | inter | Glu | 71 | Amphid interneuron | `AB praapappa` |
| `AIML` | AIM | Anterior Interneuron M Left | inter | Glu | 30 | Ring interneuron | `AB plpaapppa` |
| `AIMR` | AIM | Anterior Interneuron M Right | inter | Glu | 35 | Ring interneuron | `AB prpaapppa` |
| `AINL` | AIN | Anterior Interneuron N Left | inter | ACh | 26 | Ring interneuron | `AB alaaaalal` |
| `AINR` | AIN | Anterior Interneuron N Right | inter | ACh | 30 | Ring interneuron | `AB alaapaaar` |
| `AIYL` | AIY | Anterior Interneuron Y Left | inter | ACh | 28 | Amphid interneuron | `AB plpapaaap` |
| `AIYR` | AIY | Anterior Interneuron Y Right | inter | ACh | 30 | Amphid interneuron | `AB prpapaaap` |
| `AIZL` | AIZ | Anterior Interneuron Z Left | inter | Glu | 52 | Amphid interneuron | `AB plapaaapav` |
| `AIZR` | AIZ | Anterior Interneuron Z Right | inter | Glu | 39 | Amphid interneuron | `AB prapaaapav` |
| `ALA` | ALA | Anterior Lateral Neuron A | inter | GABA | 30 | Neuron, sends processes laterally and along dorsal cord | `AB alapppaaa` |
| `ALML` | ALM | Anterior Lateral Microtubule Neuron Left | sensory | Glu | 34 | Anterior lateral microtubule cell | `AB arppaappa` |
| `ALMR` | ALM | Anterior Lateral Microtubule Neuron Right | sensory | Glu | 16 | Anterior lateral microtubule cell | `AB arpppappa` |
| `ALNL` | ALN | Anterior Lateral Neuron N Left | sensory | ACh | 16 | Neuron associated with ALM | `AB plapappppap` |
| `ALNR` | ALN | Anterior Lateral Neuron N Right | sensory | ACh | 20 | Neuron associated with ALM | `AB prapappppap` |
| `AQR` | AQR | Anterior, Q-cell Derived Receptor | sensory | Glu | 34 | Neuron, basal body. Not part of a sensillum, projects into ring | `QR.ap` |
| `AS1` | AS | A-type Short Motor Neuron 1 | motor | ACh | 25 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P1.apa` |
| `AS10` | AS | A-type Short Motor Neuron 10 | motor | ACh | 23 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P10.apa` |
| `AS11` | AS | A-type Short Motor Neuron 11 | motor | ACh | 32 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P11.apa` |
| `AS2` | AS | A-type Short Motor Neuron 2 | motor | ACh | 22 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P2.apa` |
| `AS3` | AS | A-type Short Motor Neuron 3 | motor | ACh | 16 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P3.apa` |
| `AS4` | AS | A-type Short Motor Neuron 4 | motor | ACh | 20 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P4.apa` |
| `AS5` | AS | A-type Short Motor Neuron 5 | motor | ACh | 18 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P5.apa` |
| `AS6` | AS | A-type Short Motor Neuron 6 | motor | ACh | 16 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P6.apa` |
| `AS7` | AS | A-type Short Motor Neuron 7 | motor | ACh | 11 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P7.apa` |
| `AS8` | AS | A-type Short Motor Neuron 8 | motor | ACh | 11 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P8.apa` |
| `AS9` | AS | A-type Short Motor Neuron 9 | motor | ACh | 9 | Ventral cord motor neuron, innervates dorsal muscles, no ventral counterpart | `P9.apa` |
| `ASEL` | ASE | Amphid Single Cilium E Left | sensory | Glu | 35 | Amphid neurons, single ciliated endings | `AB alppppppaa` |
| `ASER` | ASE | Amphid Single Cilium E Right | sensory | Glu | 44 | Amphid neurons, single ciliated endings | `AB praaapppaa` |
| `ASGL` | ASG | Amphid Single Cilium G Left | sensory | Glu | 17 | Amphid neurons, single ciliated endings | `AB plaapapap` |
| `ASGR` | ASG | Amphid Single Cilium G Right | sensory | Glu | 23 | Amphid neurons, single ciliated endings | `AB praapapap` |
| `ASHL` | ASH | Amphid Single Cilium H Left | sensory | Glu | 37 | Amphid neurons, single ciliated endings | `AB plpaappaa` |
| `ASHR` | ASH | Amphid Single Cilium H Right | sensory | Glu | 50 | Amphid neurons, single ciliated endings | `AB prpaappaa` |
| `ASIL` | ASI | Amphid Single Cilium I Left | sensory | betaine † | 26 | Amphid neurons, single ciliated endings | `AB plaapapppa` |
| `ASIR` | ASI | Amphid Single Cilium I Right | sensory | betaine † | 23 | Amphid neurons, single ciliated endings | `AB praapapppa` |
| `ASJL` | ASJ | Amphid Single Cilium J Left | sensory | ACh | 20 | Amphid neurons, single ciliated endings | `AB alpppppppa` |
| `ASJR` | ASJ | Amphid Single Cilium J Right | sensory | ACh | 10 | Amphid neurons, single ciliated endings | `AB praaappppa` |
| `ASKL` | ASK | Amphid Single Cilium K Left | sensory | Glu | 21 | Amphid neurons, single ciliated endings | `AB alpppapppa` |
| `ASKR` | ASK | Amphid Single Cilium K Right | sensory | Glu | 31 | Amphid neurons, single ciliated endings | `AB praaaapppa` |
| `AUAL` | AUA | Amphid-associated Unknown Receptor A Left | inter | Glu | 25 | Neuron, process runs with amphid processes but lacks ciliated ending | `AB alpppppppp` |
| `AUAR` | AUA | Amphid-associated Unknown Receptor A Right | inter | Glu | 25 | Neuron, process runs with amphid processes but lacks ciliated ending | `AB praaappppp` |
| `AVAL` | AVA | Anterior Ventral Process A Left | inter | ACh | 154 | Ventral cord interneuron | `AB alppaaapa` |
| `AVAR` | AVA | Anterior Ventral Process A Right | inter | ACh | 159 | Ventral cord interneuron | `AB alaappapa` |
| `AVBL` | AVB | Anterior Ventral Process B Left | inter | ACh | 95 | Ventral cord interneuron | `AB plpaapaap` |
| `AVBR` | AVB | Anterior Ventral Process B Right | inter | ACh | 94 | Ventral cord interneuron | `AB prpaapaap` |
| `AVDL` | AVD | Anterior Ventral Process D Left | inter | ACh | 80 | Ventral cord interneuron | `AB alaaapalr` |
| `AVDR` | AVD | Anterior Ventral Process D Right | inter | ACh | 93 | Ventral cord interneuron | `AB alaaapprl` |
| `AVEL` | AVE | Anterior Ventral Process E Left | inter | ACh | 85 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord | `AB alpppaaaa` |
| `AVER` | AVE | Anterior Ventral Process E Right | inter | ACh | 85 | Ventral cord interneuron, like AVD but outputs restricted to anterior cord | `AB praaaaaaa` |
| `AVFL` | AVF | Anterior Ventral Process F Left | inter | GABA † | 57 | Interneuron | `P1.aaaa/ W.aaa` |
| `AVFR` | AVF | Anterior Ventral Process F Right | inter | GABA † | 49 | Interneuron | `P1.aaaa/ W.aaa` |
| `AVG` | AVG | Anterior Ventral Process G | inter | ACh * | 56 | Ventral cord interneuron | `AB prpapppap` |
| `AVHL` | AVH | Anterior Ventral Process H Left | inter | — | 60 | Neuron, mainly postsynaptic in ventral cord and presynaptic in the ring | `AB alapaaaaa` |
| `AVHR` | AVH | Anterior Ventral Process H Right | inter | — | 58 | Neuron, mainly postsynaptic in ventral cord and presynaptic in the ring | `AB alappapaa` |
| `AVJL` | AVJ | Anterior Ventral Process J Left | inter | GABA * | 60 | Neuron, synapses like AVHL/R | `AB alapapppa` |
| `AVJR` | AVJ | Anterior Ventral Process J Right | inter | GABA * | 61 | Neuron, synapses like AVHL/R | `AB alapppppa` |
| `AVKL` | AVK | Anterior Ventral Process K Left | inter | — | 57 | Ring and ventral cord interneuron | `AB plpapapap` |
| `AVKR` | AVK | Anterior Ventral Process K Right | inter | — | 50 | Ring and ventral cord interneuron | `AB prpapapap` |
| `AVL` | AVL | Anterior Ventral Process L | inter | GABA | 57 | Ring and ventral cord interneuron and an excitatory GABAergic motor neuron for rectal muscles. Few synapses | `AB prpappaap` |
| `AVM` | AVM | Anterior Ventral Microtubule Neuron | sensory | Glu | 28 | Anterior ventral microtubule cell, touch receptor | `QR.paa` |
| `AWAL` | AWA | Amphid Wing Neuron A Left | sensory | ACh * | 16 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB plaapapaa` |
| `AWAR` | AWA | Amphid Wing Neuron A Right | sensory | ACh * | 28 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB praapapaa` |
| `AWBL` | AWB | Amphid Wing Neuron B Left | sensory | ACh | 30 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB alpppppap` |
| `AWBR` | AWB | Amphid Wing Neuron B Right | sensory | ACh | 29 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB praaappap` |
| `AWCL` | AWC | Amphid Wing Neuron C Left | sensory | Glu | 29 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB plpaaaaap` |
| `AWCR` | AWC | Amphid Wing Neuron C Right | sensory | Glu | 36 | Amphid wing cells, neurons having ciliated sheet-like sensory endings closely associated with amphid sheath | `AB prpaaaaap` |
| `BAGL` | BAG | BAG-like Dendritic Ending Left | sensory | Glu | 35 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB alppappap` |
| `BAGR` | BAG | BAG-like Dendritic Ending Right | sensory | Glu | 27 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB arappppap` |
| `BDUL` | BDU | Black Vesicles, Deirid Locale, Unknown Function Left | inter | — | 22 | Neuron, process runs along excretory canal and into ring, unique darkly staining synaptic vesicles | `AB arppaappp` |
| `BDUR` | BDU | Black Vesicles, Deirid Locale, Unknown Function Right | inter | — | 25 | Neuron, process runs along excretory canal and into ring, unique darkly staining synaptic vesicles | `AB arpppappp` |
| `CANL` | CAN | Excretory CANal-associated Neuron Left | unknown | betaine † | 5 | Process runs along excretory canal, no synapses, essential for survival | `AB alapaaapa` |
| `CANR` | CAN | Excretory CANal-associated Neuron Right | unknown | betaine † | 5 | Process runs along excretory canal, no synapses, essential for survival | `AB alappappa` |
| `CEPDL` | CEP | CEPhalic Sensory Neuron Dorsal Left | sensory | DA | 47 | Cephalic neurons, contain dopamine | `AB plaaaaappa` |
| `CEPDR` | CEP | CEPhalic Sensory Neuron Dorsal Right | sensory | DA | 41 | Cephalic neurons, contain dopamine | `AB arpapaappa` |
| `CEPVL` | CEP | CEPhalic Sensory Neuron Ventral Left | sensory | DA | 42 | Cephalic neurons, contain dopamine | `AB plpaappppa` |
| `CEPVR` | CEP | CEPhalic Sensory Neuron Ventral Right | sensory | DA | 47 | Cephalic neurons, contain dopamine | `AB prpaappppa` |
| `DA1` | DA | Dorsal A-type Motor Neuron 1 | motor | ACh | 26 | Ventral cord motor neurons, innervate dorsal muscles | `AB prppapaap` |
| `DA2` | DA | Dorsal A-type Motor Neuron 2 | motor | ACh | 35 | Ventral cord motor neurons, innervate dorsal muscles | `AB plppapapa` |
| `DA3` | DA | Dorsal A-type Motor Neuron 3 | motor | ACh | 34 | Ventral cord motor neurons, innervate dorsal muscles | `AB prppapapa` |
| `DA4` | DA | Dorsal A-type Motor Neuron 4 | motor | ACh | 26 | Ventral cord motor neurons, innervate dorsal muscles | `AB plppapapp` |
| `DA5` | DA | Dorsal A-type Motor Neuron 5 | motor | ACh | 31 | Ventral cord motor neurons, innervate dorsal muscles | `AB prppapapp` |
| `DA6` | DA | Dorsal A-type Motor Neuron 6 | motor | ACh | 30 | Ventral cord motor neurons, innervate dorsal muscles | `AB plpppaaap` |
| `DA7` | DA | Dorsal A-type Motor Neuron 7 | motor | ACh | 30 | Ventral cord motor neurons, innervate dorsal muscles | `AB prpppaaap` |
| `DA8` | DA | Dorsal A-type Motor Neuron 8 | motor | ACh | 52 | Ventral cord motor neurons, innervate dorsal muscles | `AB prpapappp` |
| `DA9` | DA | Dorsal A-type Motor Neuron 9 | motor | ACh | 50 | Ventral cord motor neurons, innervate dorsal muscles | `AB plpppaaaa` |
| `DB1` | DB | Dorsal B-type Motor Neuron 1 | motor | ACh | 29 | Ventral cord motor neurons, innervate dorsal muscles, reciprocal inhibitor (aka DB1/3) | `AB plpaaaapp` |
| `DB2` | DB | Dorsal B-type Motor Neuron 2 | motor | ACh | 28 | Ventral cord motor neurons, innervate dorsal muscles, reciprocal inhibitor (aka DB3/1) | `AB arappappa` |
| `DB3` | DB | Dorsal B-type Motor Neuron 3/1 | motor | ACh | 34 | Ventral cord motor neurons, innervate dorsal muscles, reciprocal inhibitor | `AB prpaaaapp` |
| `DB4` | DB | Dorsal B-type Motor Neuron 4 | motor | ACh | 30 | Ventral cord motor neurons, innervate dorsal muscles, reciprocal inhibitor | `AB prpappapp` |
| `DB5` | DB | Dorsal B-type Motor Neuron 5 | motor | ACh | 31 | Ventral cord motor neurons, innervate dorsal muscles, reciprocal inhibitor | `AB plpapappp` |
| `DB6` | DB | Dorsal B-type Motor Neuron 6 | motor | ACh | 24 | Ventral cord motor neurons, innervate dorsal muscles, reciprocal inhibitor | `AB plppaappp` |
| `DB7` | DB | Dorsal B-type Motor Neuron 7 | motor | ACh | 45 | Ventral cord motor neurons, innervate dorsal muscles, reciprocal inhibitor | `AB prppaappp` |
| `DD1` | DD | Dorsal D-type Motor Neuron 1 | motor | GABA | 47 | Ventral cord motor neurons, reciprocal inhibitors, change synaptic pattern during L1 | `AB plppappap` |
| `DD2` | DD | Dorsal D-type Motor Neuron 2 | motor | GABA | 31 | Ventral cord motor neurons, reciprocal inhibitors, change synaptic pattern during L1 | `AB prppappap` |
| `DD3` | DD | Dorsal D-type Motor Neuron 3 | motor | GABA | 35 | Ventral cord motor neurons, reciprocal inhibitors, change synaptic pattern during L1 | `AB plppapppa` |
| `DD4` | DD | Dorsal D-type Motor Neuron 4 | motor | GABA | 27 | Ventral cord motor neurons, reciprocal inhibitors, change synaptic pattern during L1 | `AB prppapppa` |
| `DD5` | DD | Dorsal D-type Motor Neuron 5 | motor | GABA | 21 | Ventral cord motor neurons, reciprocal inhibitors, change synaptic pattern during L1 | `AB plppapppp` |
| `DD6` | DD | Dorsal D-type Motor Neuron 6 | motor | GABA | 47 | Ventral cord motor neurons, reciprocal inhibitors, change synaptic pattern during L1 | `AB prppapppp` |
| `DVA` | DVA | Dorsorectal Ganglion Ventral Process A | sensory | ACh | 100 | Ring interneurons, cell bodies in dorsorectal ganglion | `AB prppppapp` |
| `DVB` | DVB | Dorsorectal Ganglion Ventral Process B | inter | GABA | 25 | An excitatory GABAergic motor neuron/interneuron located in dorso-rectal ganglion. Innervates rectal muscles | `K.p` |
| `DVC` | DVC | Dorsorectal Ganglion Ventral Process C | inter | Glu | 56 | Ring interneurons, cell bodies in dorsorectal ganglion | `C aapaa` |
| `FLPL` | FLP | FLaP-like Dendritic Ending Left | sensory | Glu | 31 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB plapaaapad` |
| `FLPR` | FLP | FLaP-like Dendritic Ending Right | sensory | Glu | 31 | Neuron, ciliated ending in head, no supporting cells, associated with ILso | `AB prapaaapad` |
| `HSNL` | HSN | Hermaphrodite-Specific Neuron Left | motor | ACh | 58 | Hermaphrodite specific motor neurons (die in male embryo), innervate vulval muscles, serotonergic | `AB plapppappa` |
| `HSNR` | HSN | Hermaphrodite-Specific Neuron Right | motor | ACh | 70 | Hermaphrodite specific motor neurons (die in male embryo), innervate vulval muscles, serotonergic | `AB prapppappa` |
| `I1L` | I1 | Interneuron 1 (pharynx) Left | inter | ACh | 25 | Pharyngeal interneurons: ant sensory, input from RIP | `AB alpapppaa` |
| `I1R` | I1 | Interneuron 1 (pharynx) Right | inter | ACh | 31 | Pharyngeal interneurons: ant sensory, input from RIP | `AB arapappaa` |
| `I2L` | I2 | Interneuron 2 (pharynx) Left | inter | Glu | 26 | Pharyngeal interneurons, ant sensory | `AB alpappaapa` |
| `I2R` | I2 | Interneuron 2 (pharynx) Right | inter | Glu | 32 | Pharyngeal interneurons, ant sensory | `AB arapapaapa` |
| `I3` | I3 | Interneuron 3 (pharynx) | inter | ACh | 28 | Pharyngeal interneuron, ant sensory | `MS aaaaapaa` |
| `I4` | I4 | Interneuron 4 (pharynx) | inter | Glu * | 25 | Pharyngeal interneuron | `MS aaaapaa` |
| `I5` | I5 | Interneuron 5 (pharynx) | inter | Glu | 35 | Pharyngeal interneuron, post sensory | `AB arapapapp` |
| `I6` | I6 | Interneuron 6 (pharynx) | inter | — | 18 | Pharyngeal interneuron, post sensory | `MS paaapaa` |
| `IL1DL` | IL1 | Inner Labial 1 Dorsal Left | sensory | Glu | 20 | Inner labial neuron | `AB alapappaaa` |
| `IL1DR` | IL1 | Inner Labial 1 Dorsal Right | sensory | Glu | 28 | Inner labial neuron | `AB alappppaaa` |
| `IL1L` | IL1 | Inner Labial 1 Left | sensory | Glu | 32 | Inner labial neuron | `AB alapaappaa` |
| `IL1R` | IL1 | Inner Labial 1 Right | sensory | Glu | 32 | Inner labial neuron | `AB alaappppaa` |
| `IL1VL` | IL1 | Inner Labial 1 Ventral Left | sensory | Glu | 23 | Inner labial neuron | `AB alppapppaa` |
| `IL1VR` | IL1 | Inner Labial 1 Ventral Right | sensory | Glu | 22 | Inner labial neuron | `AB arapppppaa` |
| `IL2DL` | IL2 | Inner Labial 2 Dorsal Left | sensory | ACh | 18 | Inner labial neuron | `AB alapappap` |
| `IL2DR` | IL2 | Inner Labial 2 Dorsal Right | sensory | ACh | 15 | Inner labial neuron | `AB alappppap` |
| `IL2L` | IL2 | Inner Labial 2 Left | sensory | ACh | 30 | Inner labial neuron | `AB alapaappp` |
| `IL2R` | IL2 | Inner Labial 2 Right | sensory | ACh | 24 | Inner labial neuron | `AB alaappppp` |
| `IL2VL` | IL2 | Inner Labial 2 Ventral Left | sensory | ACh | 14 | Inner labial neuron | `AB alppapppp` |
| `IL2VR` | IL2 | Inner Labial 2 Ventral Right | sensory | ACh | 22 | Inner labial neuron | `AB arapppppp` |
| `LUAL` | LUA | LUmbar Ganglion A Left | inter | Glu | 23 | Interneuron, short process in post ventral cord | `AB plpppaapap` |
| `LUAR` | LUA | LUmbar Ganglion A Right | inter | Glu | 29 | Interneuron, short process in post ventral cord | `AB prpppaapap` |
| `M1` | M1 | Motor Neuron 1 (pharynx) | motor | ACh | 29 | Pharyngeal motorneuron | `MS paapaaa` |
| `M2L` | M2 | Motor Neuron 2 (pharynx) Left | motor | ACh | 24 | Pharyngeal motorneurons | `AB araapappa` |
| `M2R` | M2 | Motor Neuron 2 (pharynx) Right | motor | ACh | 21 | Pharyngeal motorneurons | `AB araappppa` |
| `M3L` | M3 | Motor Neuron 3 (pharynx) Left | motor | Glu | 22 | Pharyngeal sensory-motorneurons | `AB araapappp` |
| `M3R` | M3 | Motor Neuron 3 (pharynx) Right | motor | Glu | 23 | Pharyngeal sensory-motorneurons | `AB araappppp` |
| `M4` | M4 | Motor Neuron 4 (pharynx) | motor | ACh | 33 | Pharyngeal motorneuron | `MS paaaaaa` |
| `M5` | M5 | Motor Neuron 5 (pharynx) | motor | ACh | 22 | Pharyngeal motorneuron | `MS paaapap` |
| `MCL` | MC | Marginal Cell Neuron (pharynx) Left | unknown | ACh | 9 | Pharyngeal neurons that synapse onto marginal cells | `AB alpaaappp` |
| `MCR` | MC | Marginal Cell Neuron (pharynx) Right | unknown | ACh | 10 | Pharyngeal neurons that synapse onto marginal cells | `AB arapaappp` |
| `MI` | MI | Motor/Interneuron (pharynx) | unknown | Glu | 24 | Pharyngeal motor neuron/interneuron | `AB araappaaa` |
| `NSML` | NSM | NeuroSecretory Motor Neuron (pharynx) Left | unknown | 5-HT | 23 | Pharyngeal neurosecretory motorneuron, contain serotonin | `AB araapapaav` |
| `NSMR` | NSM | NeuroSecretory Motor Neuron (pharynx) Right | unknown | 5-HT | 26 | Pharyngeal neurosecretory motorneuron, contain serotonin | `AB araapppaav` |
| `OLLL` | OLL | Outer Labial Lateral Dendrite Left | sensory | Glu | 30 | Lateral outer labial neurons | `AB alppppapaa` |
| `OLLR` | OLL | Outer Labial Lateral Dendrite Right | sensory | Glu | 42 | Lateral outer labial neurons | `AB praaapapaa` |
| `OLQDL` | OLQ | Outer Labial Quadrant Dendrite Dorsal Left | sensory | Glu | 21 | Quadrant outer labial neuron | `AB alapapapaa` |
| `OLQDR` | OLQ | Outer Labial Quadrant Dendrite Dorsal Right | sensory | Glu | 24 | Quadrant outer labial neuron | `AB alapppapaa` |
| `OLQVL` | OLQ | Outer Labial Quadrant Dendrite Ventral Left | sensory | Glu | 26 | Quadrant outer labial neuron | `AB plpaaappaa` |
| `OLQVR` | OLQ | Outer Labial Quadrant Dendrite Ventral Right | sensory | Glu | 27 | Quadrant outer labial neuron | `AB prpaaappaa` |
| `PDA` | PDA | Preanal Cell Body Dorsal Axon A | motor | ACh | 25 | Motor neuron, process in dorsal cord, same as Y cell in hermaphrodite, Y.a in male | `AB prpppaaaa` |
| `PDB` | PDB | Preanal Cell Body Dorsal Axon B | motor | ACh | 19 | Motor neuron, process in dorsal cord, cell body in pre-anal ganglion | `P12.apa` |
| `PDEL` | PDE | Posterior DEirid Left | sensory | DA | 24 | Neuron, dopaminergic of postderid sensillum | `V5L.paaa` |
| `PDER` | PDE | Posterior DEirid Right | sensory | DA | 28 | Neuron, dopaminergic of postderid sensillum | `V5R.paaa` |
| `PHAL` | PHA | PHasmid Neuron A Left | sensory | Glu * | 20 | Phasmid neurons, chemosensory | `AB plpppaapp` |
| `PHAR` | PHA | PHasmid Neuron A Right | sensory | Glu * | 24 | Phasmid neurons, chemosensory | `AB prpppaapp` |
| `PHBL` | PHB | PHasmid Neuron B Left | sensory | Glu | 26 | Phasmid neurons, chemosensory | `AB plapppappp` |
| `PHBR` | PHB | PHasmid Neuron B Right | sensory | Glu | 20 | Phasmid neurons, chemosensory | `AB prapppappp` |
| `PHCL` | PHC | PHasmid Neuron C Left | sensory | Glu | 19 | Neuron, striated rootlet in male, possibly sensory in tail spike | `TL.pppaa` |
| `PHCR` | PHC | PHasmid Neuron C Right | sensory | Glu | 21 | Neuron, striated rootlet in male, possibly sensory in tail spike | `TR.pppaa` |
| `PLML` | PLM | Posterior Lateral Microtubule Neuron Left | sensory | Glu | 18 | Posterior lateral microtubule cells, touch receptor neurons | `AB plapappppaa` |
| `PLMR` | PLM | Posterior Lateral Microtubule Neuron Right | sensory | Glu | 14 | Posterior lateral microtubule cells, touch receptor neurons | `AB prapappppaa` |
| `PLNL` | PLN | Posterior Lateral N Left | sensory | ACh | 10 | Interneuron, associated with PLM | `TL.pppap` |
| `PLNR` | PLN | Posterior Lateral N Right | sensory | ACh | 13 | Interneuron, associated with PLM | `TR.pppap` |
| `PQR` | PQR | Posterior Q-cell Derived Receptor | sensory | Glu | 31 | Neuron, basal body, not part of a sensillum, projects into preanal ganglion | `QL.ap` |
| `PVCL` | PVC | Posterior Ventral Process C Left | inter | ACh | 95 | Ventral cord interneuron, cell body in lumbar ganglion, synapses onto VB and DB motor neurons, formerly called delta | `AB plpppaapaa` |
| `PVCR` | PVC | Posterior Ventral Process C Right | inter | ACh | 99 | Ventral cord interneuron, cell body in lumbar ganglion, synapses onto VB and DB motor neurons, formerly called delta | `AB prpppaapaa` |
| `PVDL` | PVD | Posterior Ventral Process D Left | sensory | Glu | 32 | Neuron, lateral process adjacent to excretory canal | `V5L.paapa` |
| `PVDR` | PVD | Posterior Ventral Process D Right | sensory | Glu | 29 | Neuron, lateral process adjacent to excretory canal | `V5R.paapa` |
| `PVM` | PVM | Posterior Ventral Microtubule Neuron | sensory | — | 25 | Posterior ventral microtubule cell, touch receptor | `QL.paa` |
| `PVNL` | PVN | Posterior Ventral Process N Left | inter | ACh | 51 | Interneuron/motor neuron, post. vent. cord, few synapses | `TL.appp` |
| `PVNR` | PVN | Posterior Ventral Process N Right | inter | ACh | 70 | Interneuron/motor neuron, post. vent. cord, few synapses | `TR.appp` |
| `PVPL` | PVP | Posterior Ventral Process P Left | inter | ACh | 60 | Interneuron, cell body in preanal ganglion, projects along ventral cord to nerve ring | `AB plppppaaa` |
| `PVPR` | PVP | Posterior Ventral Process P Right | inter | ACh | 62 | Interneuron, cell body in preanal ganglion, projects along ventral cord to nerve ring | `AB prppppaaa` |
| `PVQL` | PVQ | Posterior Ventral Process Q Left | inter | — | 33 | Interneuron, projects along ventral cord to ring | `AB plapppaaa` |
| `PVQR` | PVQ | Posterior Ventral Process Q Right | inter | — | 45 | Interneuron, projects along ventral cord to ring | `AB prapppaaa` |
| `PVR` | PVR | Posterior Ventral Process R | inter | Glu | 70 | Interneuron, projects along ventral cord to ring | `C aappa` |
| `PVT` | PVT | Posterior Ventral Process T | inter | — | 64 | Interneuron, projects along ventral cord to ring | `AB plpappppa` |
| `PVWL` | PVW | Posterior Ventral Process W Left | inter | — | 21 | Interneuron, posterior ventral cord, few synapses | `TL.ppa` |
| `PVWR` | PVW | Posterior Ventral Process W Right | inter | — | 26 | Interneuron, posterior ventral cord, few synapses | `TR.ppa` |
| `RIAL` | RIA | Ring Interneuron A Left | inter | Glu | 67 | Ring interneuron, many synapses | `AB alapaapaa` |
| `RIAR` | RIA | Ring Interneuron A Right | inter | Glu | 69 | Ring interneuron, many synapses | `AB alaapppaa` |
| `RIBL` | RIB | Ring Interneuron B Left | inter | GABA | 54 | Ring interneuron | `AB plpaappap` |
| `RIBR` | RIB | Ring Interneuron B Right | inter | GABA | 55 | Ring interneuron | `AB prpaappap` |
| `RICL` | RIC | Ring Interneuron C Left | inter | OA | 40 | Ring interneuron | `AB plppaaaapp` |
| `RICR` | RIC | Ring Interneuron C Right | inter | OA | 48 | Ring interneuron | `AB prppaaaapp` |
| `RID` | RID | Ring Interneuron D | inter | — | 40 | Ring interneuron, projects along dorsal cord | `AB alappaapa` |
| `RIFL` | RIF | Ring Interneuron F Left | inter | ACh | 22 | Ring interneuron | `AB plppapaaap` |
| `RIFR` | RIF | Ring Interneuron F Right | inter | ACh | 38 | Ring interneuron | `AB prppapaaap` |
| `RIGL` | RIG | Ring Interneuron G Left | inter | Glu | 54 | Ring interneuron | `AB plppappaa` |
| `RIGR` | RIG | Ring Interneuron G Right | inter | Glu | 48 | Ring interneuron | `AB prppappaa` |
| `RIH` | RIH | Ring Interneuron H | inter | ACh | 59 | Ring interneuron | `AB prpappaaa` |
| `RIML` | RIM | Ring Interneuron M Left | inter/motor | Glu | 62 | Ring motor neuron | `AB plppaapap` |
| `RIMR` | RIM | Ring Interneuron M Right | inter/motor | Glu | 58 | Ring motor neuron | `AB prppaapap` |
| `RIPL` | RIP | Ring Interneuron P Left | inter | ACh | 34 | Ring/pharynx interneuron, only direct connection between pharynx and ring | `AB alpapaaaa` |
| `RIPR` | RIP | Ring Interneuron P Right | inter | ACh | 34 | Ring/pharynx interneuron, only direct connection between pharynx and ring | `AB arappaaaa` |
| `RIR` | RIR | Ring Interneuron R | inter | ACh | 49 | Ring interneuron | `AB prpapppaa` |
| `RIS` | RIS | Ring Interneuron S | inter | GABA | 56 | Ring interneuron | `AB prpappapa` |
| `RIVL` | RIV | Ring Interneuron V Left | motor | ACh | 31 | Ring interneuron | `AB plpaapaaa` |
| `RIVR` | RIV | Ring Interneuron V Right | motor | ACh | 37 | Ring interneuron | `AB prpaapaaa` |
| `RMDDL` | RMD | Ring Motor Neuron D Dorsal Left | motor | ACh | 50 | Ring motor neuron/interneuron, many synapses | `AB alpapapaa` |
| `RMDDR` | RMD | Ring Motor Neuron D Dorsal Right | motor | ACh | 44 | Ring motor neuron/interneuron, many synapses | `AB arappapaa` |
| `RMDL` | RMD | Ring Motor Neuron D Left | motor | ACh | 62 | Ring motor neuron/interneuron, many synapses | `AB alpppapad` |
| `RMDR` | RMD | Ring Motor Neuron D Right | motor | ACh | 49 | Ring motor neuron/interneuron, many synapses | `AB praaaapad` |
| `RMDVL` | RMD | Ring Motor Neuron D Ventral Left | motor | ACh | 46 | Ring motor neuron/interneuron, many synapses | `AB alppapaaa` |
| `RMDVR` | RMD | Ring Motor Neuron D Ventral Right | motor | ACh | 59 | Ring motor neuron/interneuron, many synapses | `AB arapppaaa` |
| `RMED` | RME | Ring Motor Neuron E Dorsal | motor | GABA | 32 | Ring motor neuron | `AB alapppaap` |
| `RMEL` | RME | Ring Motor Neuron E Left | motor | GABA | 24 | Ring motor neuron | `AB alaaaarlp` |
| `RMER` | RME | Ring Motor Neuron E Right | motor | GABA | 30 | Ring motor neuron | `AB alaaaarrp` |
| `RMEV` | RME | Ring Motor Neuron E Ventral | motor | GABA | 44 | Ring motor neuron | `AB plpappaaa` |
| `RMFL` | RMF | Ring Motor Neuron F Left | inter | ACh | 30 | Ring motor neuron/interneuron | `G2.al` |
| `RMFR` | RMF | Ring Motor Neuron F Right | inter | ACh | 16 | Ring motor neuron/interneuron | `G2.ar` |
| `RMGL` | RMG | Ring Motor Neuron G Left | inter | — | 63 | Ring interneuron | `AB plapaaapp` |
| `RMGR` | RMG | Ring Motor Neuron G Right | inter | — | 53 | Ring interneuron | `AB prapaaapp` |
| `RMHL` | RMH | Ring Motor Neuron H Left | motor | ACh | 33 | Ring motor neuron/interneuron | `G1.l` |
| `RMHR` | RMH | Ring Motor Neuron H Right | motor | ACh | 29 | Ring motor neuron/interneuron | `G1.r` |
| `SAADL` | SAA | Sublateral Anterior A Dorsal Left | inter | ACh | 33 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB alppapapa` |
| `SAADR` | SAA | Sublateral Anterior A Dorsal Right | inter | ACh | 35 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB arapppapa` |
| `SAAVL` | SAA | Sublateral Anterior A Ventral Left | inter | ACh | 31 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB plpaaaaaa` |
| `SAAVR` | SAA | Sublateral Anterior A Ventral Right | inter | ACh | 40 | Ring interneuron, anteriorly projecting process that runs sublaterally | `AB prpaaaaaa` |
| `SABD` | SAB | Sublateral Anterior B Dorsal | inter/motor | ACh | 29 | Ring interneuron, anteriorly projecting process that runs sublaterally, synapses to anterior body muscles in L1 | `AB plppapaap` |
| `SABVL` | SAB | Sublateral Anterior B Ventral Left | inter/motor | ACh | 15 | Ring interneuron, anteriorly projecting process that runs sublaterally, synapses to anterior body muscles in L1 | `AB plppapaaaa` |
| `SABVR` | SAB | Sublateral Anterior B Ventral Right | inter/motor | ACh | 16 | Ring interneuron, anteriorly projecting process that runs sublaterally, synapses to anterior body muscles in L1 | `AB prppapaaaa` |
| `SDQL` | SDQ | Sublateral Dorsal Q-cell Derived Left | sensory | ACh | 30 | Post. lateral interneuron, process projects into ring | `QL.pap` |
| `SDQR` | SDQ | Sublateral Dorsal Q-cell Derived Right | sensory | ACh | 27 | Ant. lateral interneuron, process projects into ring | `QR.pap` |
| `SIADL` | SIA | Sublateral Interneuron A Dorsal Left | inter/motor | ACh | 36 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB plpapaapa` |
| `SIADR` | SIA | Sublateral Interneuron A Dorsal Right | inter/motor | ACh | 42 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB prpapaapa` |
| `SIAVL` | SIA | Sublateral Interneuron A Ventral Left | inter/motor | ACh | 38 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB plpapappa` |
| `SIAVR` | SIA | Sublateral Interneuron A Ventral Right | inter/motor | ACh | 35 | Receive a few synapses in the ring, have posteriorly directed processes that run sublaterally | `AB prpapappa` |
| `SIBDL` | SIB | Sublateral Interneuron B Dorsal Left | inter/motor | ACh | 28 | Similar to SIA | `AB plppaaaaa` |
| `SIBDR` | SIB | Sublateral Interneuron B Dorsal Right | inter/motor | ACh | 32 | Similar to SIA | `AB prppaaaaa` |
| `SIBVL` | SIB | Sublateral Interneuron B Ventral Left | inter/motor | ACh | 37 | Similar to SIA | `AB plpapaapp` |
| `SIBVR` | SIB | Sublateral Interneuron B Ventral Right | inter/motor | ACh | 38 | Similar to SIA | `AB prpapaapp` |
| `SMBDL` | SMB | Sublateral Motor Neuron B Dorsal Left | motor | ACh | 48 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB alpapapapp` |
| `SMBDR` | SMB | Sublateral Motor Neuron B Dorsal Right | motor | ACh | 59 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB arappapapp` |
| `SMBVL` | SMB | Sublateral Motor Neuron B Ventral Left | motor | ACh | 50 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB alpapappp` |
| `SMBVR` | SMB | Sublateral Motor Neuron B Ventral Right | motor | ACh | 54 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB arappappp` |
| `SMDDL` | SMD | Sublateral Motor Neuron D Dorsal Left | motor | ACh | 72 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB plpapaaaa` |
| `SMDDR` | SMD | Sublateral Motor Neuron D Dorsal Right | motor | ACh | 75 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB prpapaaaa` |
| `SMDVL` | SMD | Sublateral Motor Neuron D Ventral Left | motor | ACh | 68 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB alppappaa` |
| `SMDVR` | SMD | Sublateral Motor Neuron D Ventral Right | motor | ACh | 58 | Ring motor neuron/interneuron, has a posteriorly directed process that runs sublaterally | `AB arappppaa` |
| `URADL` | URA | Unknown Receptor, not Ciliated A Dorsal Left | motor | ACh | 13 | Ring motor neuron | `AB plaaaaaaa` |
| `URADR` | URA | Unknown Receptor, not Ciliated A Dorsal Right | motor | ACh | 17 | Ring motor neuron | `AB arpapaaaa` |
| `URAVL` | URA | Unknown Receptor, not Ciliated A Ventral Left | motor | ACh | 28 | Ring motor neuron | `AB plpaaapaa` |
| `URAVR` | URA | Unknown Receptor, not Ciliated A Ventral Right | motor | ACh | 25 | Ring motor neuron | `AB prpaaapaa` |
| `URBL` | URB | Unknown Receptor, not Ciliated B Left | inter | ACh | 26 | Neuron, presynaptic in ring, ending in head | `AB plaapaapa` |
| `URBR` | URB | Unknown Receptor, not Ciliated B Right | inter | ACh | 38 | Neuron, presynaptic in ring, ending in head | `AB praapaapa` |
| `URXL` | URX | Unknown Receptor, not Ciliated X Left | sensory | ACh | 45 | Ring interneuron | `AB plaaaaappp` |
| `URXR` | URX | Unknown Receptor, not Ciliated X Right | sensory | ACh | 52 | Ring interneuron | `AB arpapaappp` |
| `URYDL` | URY | Unknown Receptor, not Ciliated Y Dorsal Left | sensory | Glu | 21 | Neuron, presynaptic in ring, ending in head | `AB alapapapp` |
| `URYDR` | URY | Unknown Receptor, not Ciliated Y Dorsal Right | sensory | Glu | 24 | Neuron, presynaptic in ring, ending in head | `AB alapppapp` |
| `URYVL` | URY | Unknown Receptor, not Ciliated Y Ventral Left | sensory | Glu | 26 | Neuron, presynaptic in ring, ending in head | `AB plpaaappp` |
| `URYVR` | URY | Unknown Receptor, not Ciliated Y Ventral Right | sensory | Glu | 26 | Neuron, presynaptic in ring, ending in head | `AB prpaaappp` |
| `VA1` | VA | Ventral A-type Motor Neuron 1 | motor | ACh | 18 | Ventral cord motor neuron, innervates vent. body muscles | `W.pa` |
| `VA10` | VA | Ventral A-type Motor Neuron 10 | motor | ACh | 31 | Ventral cord motor neuron, innervates vent. body muscles | `P10.aaaa` |
| `VA11` | VA | Ventral A-type Motor Neuron 11 | motor | ACh | 35 | Ventral cord motor neuron, innervates vent. body muscles | `P11.aaaa` |
| `VA12` | VA | Ventral A-type Motor Neuron 12 | motor | ACh | 52 | Ventral cord motor neuron, innervates vent. body muscles, but also interneuron in preanal ganglion | `P12.aaaa` |
| `VA2` | VA | Ventral A-type Motor Neuron 2 | motor | ACh | 32 | Ventral cord motor neuron, innervates vent. body muscles | `P2.aaaa` |
| `VA3` | VA | Ventral A-type Motor Neuron 3 | motor | ACh | 31 | Ventral cord motor neuron, innervates vent. body muscles | `P3.aaaa` |
| `VA4` | VA | Ventral A-type Motor Neuron 4 | motor | ACh | 27 | Ventral cord motor neuron, innervates vent. body muscles | `P4.aaaa` |
| `VA5` | VA | Ventral A-type Motor Neuron 5 | motor | ACh | 22 | Ventral cord motor neuron, innervates vent. body muscles | `P5.aaaa` |
| `VA6` | VA | Ventral A-type Motor Neuron 6 | motor | ACh | 28 | Ventral cord motor neuron, innervates vent. body muscles | `P6.aaaa` |
| `VA7` | VA | Ventral A-type Motor Neuron 7 | motor | ACh | 23 | Ventral cord motor neuron, innervates vent. body muscles | `P7.aaaa` |
| `VA8` | VA | Ventral A-type Motor Neuron 8 | motor | ACh | 25 | Ventral cord motor neuron, innervates vent. body muscles | `P8.aaaa` |
| `VA9` | VA | Ventral A-type Motor Neuron 9 | motor | ACh | 27 | Ventral cord motor neuron, innervates vent. body muscles | `P9.aaaa` |
| `VB1` | VB | Ventral B-type Motor Neuron 1 | motor | ACh | 40 | Ventral cord motor neuron, innervates vent. body muscles, also interneuron in ring | `P1.aaap` |
| `VB10` | VB | Ventral B-type Motor Neuron 10 | motor | ACh | 24 | Ventral cord motor neuron, innervates vent. body muscles | `P9.aaap` |
| `VB11` | VB | Ventral B-type Motor Neuron 11 | motor | ACh | 29 | Ventral cord motor neuron, innervates vent. body muscles | `P10.aaap` |
| `VB2` | VB | Ventral B-type Motor Neuron 2 | motor | ACh | 31 | Ventral cord motor neuron, innervates vent. body muscles | `W.aap` |
| `VB3` | VB | Ventral B-type Motor Neuron 3 | motor | ACh | 25 | Ventral cord motor neuron, innervates vent. body muscles | `P2.aaap` |
| `VB4` | VB | Ventral B-type Motor Neuron 4 | motor | ACh | 26 | Ventral cord motor neuron, innervates vent. body muscles | `P3.aaap` |
| `VB5` | VB | Ventral B-type Motor Neuron 5 | motor | ACh | 19 | Ventral cord motor neuron, innervates vent. body muscles | `P4.aaap` |
| `VB6` | VB | Ventral B-type Motor Neuron 6 | motor | ACh | 27 | Ventral cord motor neuron, innervates vent. body muscles | `P5.aaap` |
| `VB7` | VB | Ventral B-type Motor Neuron 7 | motor | ACh | 26 | Ventral cord motor neuron, innervates vent. body muscles | `P6.aaap` |
| `VB8` | VB | Ventral B-type Motor Neuron 8 | motor | ACh | 21 | Ventral cord motor neuron, innervates vent. body muscles | `P7.aaap` |
| `VB9` | VB | Ventral B-type Motor Neuron 9 | motor | ACh | 18 | Ventral cord motor neuron, innervates vent. body muscles | `P8.aaap` |
| `VC1` | VC | Ventral C-type Motor Neuron 1 | motor | ACh | 31 | Hermaphrodite specific ventral cord motor neuron innervates vulval muscles and ventral body muscles | `P3.aap` |
| `VC2` | VC | Ventral C-type Motor Neuron 2 | motor | ACh | 39 | Hermaphrodite specific ventral cord motor neuron innervates vulval muscles and ventral body muscles | `P4.aap` |
| `VC3` | VC | Ventral C-type Motor Neuron 3 | motor | ACh | 47 | Hermaphrodite specific ventral cord motor neuron innervates vulval muscles and ventral body muscles | `P5.aap` |
| `VC4` | VC | Ventral C-type Motor Neuron 4 | motor | ACh | 16 | Hermaphrodite specific ventral cord motor neuron innervates vulval muscles and ventral body muscles | `P6.aap` |
| `VC5` | VC | Ventral C-type Motor Neuron 5 | motor | ACh | 26 | Hermaphrodite specific ventral cord motor neuron innervates vulval muscles and ventral body muscles | `P7.aap` |
| `VC6` | VC | Ventral C-type Motor Neuron 6 | motor | ACh | 14 | Hermaphrodite specific ventral cord motor neuron innervates vulval muscles and ventral body muscles | `P8.aap` |
| `VD1` | VD | Ventral D-type Motor Neuron 1 | motor | GABA | 25 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `W.pp` |
| `VD10` | VD | Ventral D-type Motor Neuron 10 | motor | GABA | 20 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P9.app` |
| `VD11` | VD | Ventral D-type Motor Neuron 11 | motor | GABA | 32 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P10.app` |
| `VD12` | VD | Ventral D-type Motor Neuron 12 | motor | GABA | 42 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P11.app` |
| `VD13` | VD | Ventral D-type Motor Neuron 13 | motor | GABA | 44 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P12.app` |
| `VD2` | VD | Ventral D-type Motor Neuron 2 | motor | GABA | 30 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P1.app` |
| `VD3` | VD | Ventral D-type Motor Neuron 3 | motor | GABA | 27 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P2.app` |
| `VD4` | VD | Ventral D-type Motor Neuron 4 | motor | GABA | 25 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P3.app` |
| `VD5` | VD | Ventral D-type Motor Neuron 5 | motor | GABA | 28 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P4.app` |
| `VD6` | VD | Ventral D-type Motor Neuron 6 | motor | GABA | 30 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P5.app` |
| `VD7` | VD | Ventral D-type Motor Neuron 7 | motor | GABA | 25 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P6.app` |
| `VD8` | VD | Ventral D-type Motor Neuron 8 | motor | GABA | 21 | Ventralcord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P7.app` |
| `VD9` | VD | Ventral D-type Motor Neuron 9 | motor | GABA | 26 | Ventral cord motor neuron, innervates vent body muscles, reciprocal inhibitor | `P8.app` |

## Non-neuronal cells in this dataset

171 cells that are not neurons. Muscles are the nervous system's actual output, and glia are increasingly understood to shape neural function, so neither is filtered out.

| Cell | Category | Partners | Description |
| --- | --- | ---: | --- |
| `CEPshDL` | glia | 3 | Sheath cell other than amphid sheath and phasmid |
| `CEPshDR` | glia | 2 | Sheath cell other than amphid sheath and phasmid |
| `CEPshVL` | glia | 6 | Sheath cell other than amphid sheath and phasmid |
| `CEPshVR` | glia | 13 | Sheath cell other than amphid sheath and phasmid |
| `g1AL` | glia | 8 | Pharyngeal glial cell |
| `g1AR` | glia | 9 | Pharyngeal glial cell |
| `g1P` | glia | 3 | Pharyngeal glial cell |
| `g2L` | glia | 2 | Pharyngeal glial cell |
| `g2R` | glia | 1 | Pharyngeal glial cell |
| `MDL01` | muscle | 12 | Head muscle |
| `MDL02` | muscle | 12 | Head muscle |
| `MDL03` | muscle | 13 | Head muscle |
| `MDL04` | muscle | 10 | Head muscle |
| `MDL05` | muscle | 16 | Head muscle |
| `MDL06` | muscle | 19 | Head muscle |
| `MDL07` | muscle | 13 | Head muscle |
| `MDL08` | muscle | 12 | Main body muscle |
| `MDL09` | muscle | 15 | Main body muscle |
| `MDL10` | muscle | 12 | Main body muscle |
| `MDL11` | muscle | 14 | Main body muscle |
| `MDL12` | muscle | 15 | Main body muscle |
| `MDL13` | muscle | 14 | Main body muscle |
| `MDL14` | muscle | 15 | Main body muscle |
| `MDL15` | muscle | 12 | Main body muscle |
| `MDL16` | muscle | 14 | Main body muscle |
| `MDL17` | muscle | 8 | Main body muscle |
| `MDL18` | muscle | 8 | Main body muscle |
| `MDL19` | muscle | 10 | Main body muscle |
| `MDL20` | muscle | 9 | Main body muscle |
| `MDL21` | muscle | 10 | Main body muscle |
| `MDL22` | muscle | 8 | Main body muscle |
| `MDL23` | muscle | 7 | Main body muscle |
| `MDL24` | muscle | 6 | Main body muscle |
| `MDR01` | muscle | 13 | Head muscle |
| `MDR02` | muscle | 10 | Head muscle |
| `MDR03` | muscle | 14 | Head muscle |
| `MDR04` | muscle | 17 | Head muscle |
| `MDR05` | muscle | 17 | Head muscle |
| `MDR06` | muscle | 17 | Head muscle |
| `MDR07` | muscle | 19 | Head muscle |
| `MDR08` | muscle | 11 | Main body muscle |
| `MDR09` | muscle | 13 | Main body muscle |
| `MDR10` | muscle | 13 | Main body muscle |
| `MDR11` | muscle | 14 | Main body muscle |
| `MDR12` | muscle | 15 | Main body muscle |
| `MDR13` | muscle | 14 | Main body muscle |
| `MDR14` | muscle | 14 | Main body muscle |
| `MDR15` | muscle | 13 | Main body muscle |
| `MDR16` | muscle | 14 | Main body muscle |
| `MDR17` | muscle | 8 | Main body muscle |
| `MDR18` | muscle | 8 | Main body muscle |
| `MDR19` | muscle | 10 | Main body muscle |
| `MDR20` | muscle | 8 | Main body muscle |
| `MDR21` | muscle | 9 | Main body muscle |
| `MDR22` | muscle | 9 | Main body muscle |
| `MDR23` | muscle | 8 | Main body muscle |
| `MDR24` | muscle | 6 | Main body muscle |
| `MVL01` | muscle | 16 | Head muscle |
| `MVL02` | muscle | 12 | Head muscle |
| `MVL03` | muscle | 12 | Head muscle |
| `MVL04` | muscle | 12 | Head muscle |
| `MVL05` | muscle | 11 | Head muscle |
| `MVL06` | muscle | 10 | Head muscle |
| `MVL07` | muscle | 17 | Head muscle |
| `MVL08` | muscle | 15 | Main body muscle |
| `MVL09` | muscle | 15 | Main body muscle |
| `MVL10` | muscle | 14 | Main body muscle |
| `MVL11` | muscle | 21 | Main body muscle |
| `MVL12` | muscle | 16 | Main body muscle |
| `MVL13` | muscle | 19 | Main body muscle |
| `MVL14` | muscle | 14 | Main body muscle |
| `MVL15` | muscle | 17 | Main body muscle |
| `MVL16` | muscle | 19 | Main body muscle |
| `MVL17` | muscle | 9 | Main body muscle |
| `MVL18` | muscle | 9 | Main body muscle |
| `MVL19` | muscle | 14 | Main body muscle |
| `MVL20` | muscle | 11 | Main body muscle |
| `MVL21` | muscle | 11 | Main body muscle |
| `MVL22` | muscle | 8 | Main body muscle |
| `MVL23` | muscle | 4 | Main body muscle |
| `MVR01` | muscle | 4 | Head muscle |
| `MVR02` | muscle | 9 | Head muscle |
| `MVR03` | muscle | 12 | Head muscle |
| `MVR04` | muscle | 12 | Head muscle |
| `MVR05` | muscle | 9 | Head muscle |
| `MVR06` | muscle | 16 | Head muscle |
| `MVR07` | muscle | 16 | Head muscle |
| `MVR08` | muscle | 18 | Main body muscle |
| `MVR09` | muscle | 14 | Main body muscle |
| `MVR10` | muscle | 15 | Main body muscle |
| `MVR11` | muscle | 20 | Main body muscle |
| `MVR12` | muscle | 15 | Main body muscle |
| `MVR13` | muscle | 18 | Main body muscle |
| `MVR14` | muscle | 16 | Main body muscle |
| `MVR15` | muscle | 15 | Main body muscle |
| `MVR16` | muscle | 17 | Main body muscle |
| `MVR17` | muscle | 8 | Main body muscle |
| `MVR18` | muscle | 9 | Main body muscle |
| `MVR19` | muscle | 8 | Main body muscle |
| `MVR20` | muscle | 15 | Main body muscle |
| `MVR21` | muscle | 13 | Main body muscle |
| `MVR22` | muscle | 10 | Main body muscle |
| `MVR23` | muscle | 7 | Main body muscle |
| `MVR24` | muscle | 3 | Main body muscle |
| `mu_anal` | muscle | 4 | Anal/sphincter muscle |
| `mu_intL` | muscle | 5 | Intestinal muscles |
| `mu_intR` | muscle | 6 | Intestinal muscles |
| `mu_sph` | muscle | 3 | Anal/sphincter muscle |
| `pm1` | muscle | 7 | Pharyngeal muscle |
| `pm2D` | muscle | 3 | Pharyngeal muscle |
| `pm2VL` | muscle | 3 | Pharyngeal muscle |
| `pm2VR` | muscle | 2 | Pharyngeal muscle |
| `pm3D` | muscle | 5 | Pharyngeal muscle |
| `pm3VL` | muscle | 9 | Pharyngeal muscle |
| `pm3VR` | muscle | 7 | Pharyngeal muscle |
| `pm4D` | muscle | 19 | Pharyngeal muscle |
| `pm4VL` | muscle | 13 | Pharyngeal muscle |
| `pm4VR` | muscle | 11 | Pharyngeal muscle |
| `pm5D` | muscle | 22 | Pharyngeal muscle |
| `pm5VL` | muscle | 19 | Pharyngeal muscle |
| `pm5VR` | muscle | 16 | Pharyngeal muscle |
| `pm6D` | muscle | 4 | Pharyngeal muscle |
| `pm6VL` | muscle | 6 | Pharyngeal muscle |
| `pm6VR` | muscle | 5 | Pharyngeal muscle |
| `pm7D` | muscle | 4 | Pharyngeal muscle |
| `pm7VL` | muscle | 3 | Pharyngeal muscle |
| `pm7VR` | muscle | 4 | Pharyngeal muscle |
| `pm8` | muscle | 6 | Pharyngeal muscle |
| `um1AL` | muscle | 4 | Uterine muscle |
| `um1AR` | muscle | 4 | Uterine muscle |
| `um1PL` | muscle | 4 | Uterine muscle |
| `um1PR` | muscle | 4 | Uterine muscle |
| `um2AL` | muscle | 2 | Uterine muscle |
| `um2AR` | muscle | 2 | Uterine muscle |
| `um2PL` | muscle | 2 | Uterine muscle |
| `um2PR` | muscle | 2 | Uterine muscle |
| `vm1AL` | muscle | 3 | Vulval muscle |
| `vm1AR` | muscle | 6 | Vulval muscle |
| `vm1PL` | muscle | 3 | Vulval muscle |
| `vm1PR` | muscle | 3 | Vulval muscle |
| `vm2AL` | muscle | 8 | Vulval muscle |
| `vm2AR` | muscle | 8 | Vulval muscle |
| `vm2PL` | muscle | 8 | Vulval muscle |
| `vm2PR` | muscle | 8 | Vulval muscle |
| `GLRDL` | other | 8 | GLR cell |
| `GLRDR` | other | 8 | GLR cell |
| `GLRL` | other | 6 | GLR cell |
| `GLRR` | other | 6 | GLR cell |
| `GLRVL` | other | 3 | GLR cell |
| `GLRVR` | other | 3 | GLR cell |
| `bm` | other | 4 | Pharyngeal basement membrane |
| `e2D` | other | 1 | Pharyngeal epithelium |
| `e2VL` | other | 1 | Pharyngeal epithelium |
| `e2VR` | other | 2 | Pharyngeal epithelium |
| `e3D` | other | 4 | Pharyngeal epithelium |
| `e3VL` | other | 5 | Pharyngeal epithelium |
| `e3VR` | other | 3 | Pharyngeal epithelium |
| `exc_cell` | other | 4 | Excretory cell |
| `exc_gl` | other | 6 | Excretory gland |
| `hmc` | other | 14 | Head mesodermal cell |
| `hyp` | other | 118 | Hypodermis |
| `int` | other | 3 | Intestine |
| `mc1DL` | other | 5 | Marginal cell of the pharynx |
| `mc1DR` | other | 4 | Marginal cell of the pharynx |
| `mc1V` | other | 4 | Marginal cell of the pharynx |
| `mc2DL` | other | 6 | Marginal cell of the pharynx |
| `mc2DR` | other | 7 | Marginal cell of the pharynx |
| `mc2V` | other | 4 | Marginal cell of the pharynx |
| `mc3DL` | other | 6 | Marginal cell of the pharynx |
| `mc3DR` | other | 6 | Marginal cell of the pharynx |
| `mc3V` | other | 6 | Marginal cell of the pharynx |

## Sources

| Column | Source | Confidence |
| --- | --- | --- |
| name expansion, lineage, description | WormAtlas cell listings | published annotation |
| role | WormAtlas + Cook et al. 2019 groupings | published annotation |
| class, neurotransmitter | Wang et al. 2024, eLife 12:RP95402 | published annotation |
| partner counts | Cook SJ, Jarrell TA, Brittin CA, Wang Y, Bloniarz AE, Yakovlev MA, Nguyen KCQ, Tang LT-H, Bayer EA, Duerr JS, Bulow HE, Hobert O, Hall DH, Emmons SW. Whole-animal connectomes of both Caenorhabditis elegans sexes. Nature 571:63-71 (2019). | **measured** |

Full citations and licences: [references.md](references.md). What each source does and does not establish: [model_assumptions.md](model_assumptions.md).
