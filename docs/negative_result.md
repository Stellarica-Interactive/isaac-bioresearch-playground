# What a graded connectome model does not do

**A negative result.** An embodied *C. elegans* model built from the measured
whole-animal connectome, measured neuromuscular polarity, and the standard graded
leaky-integrator dynamics does not produce a travelling undulatory wave. Nor does
one in which a published conductance-based motor neuron replaces the graded cell.
This document states that claim precisely, gives the measurements behind it, lists
the nine explanations tested and rejected, and is explicit about what it does
*not* establish.

It is not an isolated failure. **No published model obtains *C. elegans*
locomotion from measured connectivity with measured parameters** — Boyle et al.
fit their oscillator's components, OpenWorm's c302 prescribes the phase gradient
as an explicit per-segment delay ladder, and MetaWorm optimises connection
weights. §6 gives the comparison. The interesting question is therefore not why
this implementation fails, but what the field is actually supplying when it
succeeds.

It is written for people building neuromechanical models of *C. elegans*, and
assumes familiarity with the Wicks/Kunert lineage and with Boyle, Berri & Cohen
2012. The supporting detail, including every retraction along the way, is in
[`model_assumptions.md`](model_assumptions.md) §5C–§5X.

---

## 1. The claim

With the components below, forward locomotion does not emerge:

| component | source | status |
|---|---|---|
| connectivity | Cook et al. 2019, hermaphrodite, 302 neurons + 135 muscles | measured |
| neuromuscular sign | McIntire et al. 1993 and others, per-synapse citations | measured physiology |
| interneuronal sign | Fenyves et al. 2020, receptor expression | predicted |
| neurotransmitters | Wang et al. 2024 | published annotation |
| dynamics | graded leaky integrator, `C dV/dt = leak + gap + chemical + I_ext` | standard |
| proprioception | B-type motor neurons sense anterior curvature (Wen et al. 2012) | measured mechanism, assumed gain |
| body | 24-segment articulated chain, anisotropic drag | engineering |

The body is capable of the target behaviour: driven by a scripted sinusoidal wave
it covers **7.2 body lengths in 60 s** with an inter-joint phase of **+23.0°**.
Driven by the connectome, the same body oscillates at comparable amplitude and
goes essentially nowhere as a crawler.

## 2. The measurement

Undulation is characterised by three numbers, each derived from the signal's own
timescale rather than a fixed window — a correction that mattered, see §5.3.

| | amplitude | period | **inter-joint phase** | distance |
|---|---|---|---|---|
| scripted wave (positive control) | 20.99° | 2.0 s | **+23.0°** | 7.216 BL |
| connectome-driven | 16.67° | 12.0 s | **−0.8°** | 0.994 BL |

Inter-joint phase is the quantity that separates crawling from flexing in place.
The connectome-driven body has **amplitude without propagation**: every segment
does approximately the same thing at approximately the same time.

The proximate cause is measured directly. Across the eighteen B-type motor
neurons:

- **99.3%** of variance in synaptic activation lies in the first principal
  component
- **95.9%** of variance in *membrane potential* likewise
- all eighteen operating points lie within **6.7 mV** of each other

They are, to a good approximation, one signal. A travelling wave requires
neighbouring segments to differ at a given moment; these do not.

## 3. What was tested and rejected

Each of these was a live hypothesis with an argument behind it. All were measured,
not assumed away.

| # | hypothesis | test | result |
|---|---|---|---|
| 1 | Gap junctions synchronise the population | `g_gap` swept to 1 pS, a hundredth of the assumed value | **74.7%** still shared. Rejected |
| 2 | The motor neurons need bistability (Boyle et al.) | Schmitt-trigger output, dead band swept, gain-matched | Phase **+0.0°** at every width; amplitude *falls* from 16.67° to ~3.2° |
| 3 | The D-class antagonist is missing | Checked the wiring | Present and correctly signed: 121 inhibitory D→muscle edges. Not missing |
| 4 | Mutual inhibition between B-type cells is missing | Counted it | **Three** B-to-B chemical synapses exist in the animal, weight 1 each. Nothing to restore |
| 5 | The proprioceptive delay is wrong | Sensing offset swept 0–12 joints | Phase within ±0.8° of zero throughout. Oscillation exists only for offsets 1–4; propagation at none |
| 6 | Dropped synapses (29% of weight is unsigned) | All four `UnknownSignPolicy` modes | No mode produces sustained propagation; see §3.2 |
| 7 | The body cannot express a wave | Scripted wave in the same body | 7.216 BL, phase +23.0°. The body is fine |
| 8 | Coupling needs a conduction delay (c302 uses one) | One uniform delay, 5–100 ms, chemical transmission only | Phase within ±0.8° of zero at every value; shared variance 99.3% → 90–93% |
| 9 | The cells need real intrinsic dynamics | VB6's published conductance model promoted into the live network | A 22 mV transient decaying over ~10 s to a fixed point. No limit cycle |

Two further controls bound the result from both ends: with zero muscle drive the
body moves **exactly 0.000 BL** (so nothing is drift), and with isotropic drag a
real gait loses a factor of **6.1** in distance (so the drag model is doing real
work).

### 3.1 The strongest single comparison

A scripted wave with drag anisotropy **removed** (1.187 BL) still outruns the
connectome-driven body with every mechanical advantage the model offers
(0.994 BL). Whatever is missing is not in the body.

### 3.2 The unsigned synapses, and a transient that looked like a result

29% of synaptic weight carries no measured or predicted polarity. Under the
default `EXCLUDE` it is dropped. All four policies, 60 s:

| policy | amplitude | phase | distance |
|---|---|---|---|
| `exclude` (default) | 16.67° | −0.8° | 0.994 BL |
| `excitatory` | 19.57° | −0.2° | 0.979 BL |
| `neutral` | 0.07° | — | 0.020 BL |
| `inhibitory` | 3.11° | **−11.0°** | 0.192 BL |

`neutral` — keeping the connections as pure shunts — does not perturb the model,
it extinguishes it.

The `inhibitory` row read as the first propagation this model had ever produced,
tail-to-head, and it is not. Run for 120 s and windowed in ten-second blocks, the
body oscillates for the first ten seconds and then **stops completely**: joint
angles at 20 s and at 110 s are identical to 0.000°, correlated at +1.0000, with
amplitude of exactly zero thereafter. The animal locks into a rigid 15° coil, and
the distance it covers is that coil being rotated and dragged.

The −11.0° was real but transient, measured over a window three quarters of which
was a frozen body. It is recorded because the near-miss is instructive: a
sustained-behaviour question was answered with a whole-recording statistic, and a
decaying transient is exactly what that conflates. See §5.3.

### 3.3 The conductance-based test, which is the strongest of the nine

Hypothesis 9 is the one the others pointed at. The graded leaky integrator is a
relaxation system with no limit cycle under any drive, so the obvious reading of
every null above is that the *cells* are too simple — a view this project held for
a dozen sections.

Nicoletti et al. 2024 publish a conductance-based model of **VB6**, a B-type motor
neuron inside the proprioceptive loop, with thirteen channels at measured
densities. Promoted into the live 302-cell network with no external drive:

| window | mean V | swing |
|---|---|---|
| 0–5 s | −48.66 mV | 22.145 mV |
| 5–10 s | −34.58 mV | 0.573 mV |
| 10–15 s | −34.53 mV | **0.000 mV** |
| 55–60 s | −34.53 mV | **0.000 mV** |

A twenty-two millivolt transient lasting ten seconds, then a fixed point held to
three decimals for fifty more. In isolation the cell settles likewise, its
currents summing to zero.

So the failure is not that the cells are too simple. A measured cell with
voltage-gated conductances, inside the measured connectome, still relaxes. That
is a considerably stronger statement than the graded result alone, and it is the
reason this document is worth writing: the obvious fix was available, published,
and does not work.

## 4. The remaining explanation

Every rejected hypothesis is a mechanism that would let a population *hold* or
*express* a phase difference. What has not been ruled out — and what the
measurements increasingly implicate — is that **there is no phase difference to
hold, because the input is common**.

The eighteen B-type cells receive:

- the same constant AVB command drive, gap-junction coupled, carrying no rhythm
  and no spatial pattern; and
- a proprioceptive signal that is **97.1% one signal**, because the body bends in
  a single mode.

This is circular in a way the model may not escape without an additional
mechanism: no wave in the body means no difference between segments, which means
no wave in the neurons. Boyle et al. break the circularity with fitted components
(a half-body receptive field, a dorsoventral gain asymmetry, binary neurons with
hysteresis); imported individually into this model, each was measured to make it
*worse* (§5M, §5Q).

A separate structural finding supports the same reading: the graded leaky
integrator has **no limit cycle anywhere** in the parameter ranges explored. It is
a relaxation system. Any oscillation it shows is driven by the sensorimotor loop,
not intrinsic.

## 5. What this does not establish

Stated plainly, because a negative result is easy to over-read.

### 5.1 It is not a claim about the animal

*C. elegans* crawls. This is a statement about a model, and a failure here is
evidence about the model's components, not about the biology.

### 5.2 It is not exhaustive over parameters

Seven of eleven biophysical parameters remain `ASSUMED` — they come from the
published graded-model lineage, where they are order-of-magnitude choices rather
than per-neuron measurements. The proprioceptive gain in particular (20 mV per
radian) is ours. A parameter regime that produces propagation may exist and not
have been found. What can be said is that propagation did not appear anywhere in
the ranges swept, and that the failure is insensitive to the parameters most
likely to matter.

### 5.3 Five numbers here were first reported wrongly

Documented rather than quietly fixed, because each inflates confidence in the
wrong direction if not known, and because the pattern matters more than any one
of them.

Two were defects in the metrics themselves:

- The propagation metric used a fixed correlation lag. In closed form its reading
  is `2·sin(φ)·sin(ω·lag)`, so it returned the true value scaled by a constant
  nobody had chosen — **0.156** at the committed lag. A genuine gait read +0.12 on
  a scale documented as reaching +1.
- The amplitude metric used a one-second window on a body whose dominant period is
  twelve seconds, reporting 0.96° where whole cycles give 14.75°.

A third limitation is live rather than fixed: the metric summarises the *whole*
recording, so a run that behaves differently at the start and the end returns a
blend of both. That is the right answer to "what did this recording contain" and
the wrong one to "what is it doing now". §3.2 is what that looks like in practice.
Any run whose transient is a large fraction of its length should be windowed
before its numbers are quoted.

Neither changes the conclusion — re-measured correctly, the connectome body's
phase is −0.8° against the control's +23.0° — but every magnitude in earlier
analyses was on a compressed scale.

Three more were results reported from a window shorter than the transient being
measured, each later reversed:

- **Ablating the B-type motor neurons appeared to *raise* oscillation amplitude**
  (1.82° against 0.96°), which would have meant the cells supposed to drive
  locomotion were suppressing it. At sixty seconds it is 16.67° intact against
  6.69° ablated. The twenty-second window was transient.
- **`--unknown-sign inhibitory` appeared to produce the first propagation in the
  project**, phase −11.0°. At 120 s the body is frozen: joint angles at 20 s and
  110 s identical to 0.000°, correlated at +1.0000, a rigid coil being dragged.
- **A conductance-based VB6 appeared to oscillate where the graded cell is
  identically static** — 4.627 mV against 0.000 mV, measured over three seconds.
  Over sixty it is a ten-second transient decaying to a fixed point.

The last is the instructive one. The control was correct and the comparison was
sound: the graded cell really is static to three decimals. It was written up with
the explicit argument that no window-length objection could turn 0.000 into
4.627 — and that argument was wrong, because the graded cell settles instantly
while the conductance cell takes ten seconds, so any window inside those ten
seconds shows a difference that vanishes outside them.

**A control establishes what a difference is attributable to, not that the
difference persists.** The general rule, learned five times: measure the settling
time before choosing a window.

### 5.4 Missing biology that could matter

- **Extrasynaptic and neuropeptidergic signalling** (Ripoll-Sánchez et al. 2023)
  is entirely absent. It is not a small system in *C. elegans*.
- **Intrinsic conductances, for the other 301 cells.** VB6 has been tested
  (§3.3) and the rest have not. Nicoletti et al. 2024 publish seven cells, all
  seven of which now run here; the remaining 295 neurons are still graded.
  Promoting a whole class would mean applying one cell's fitted densities to its
  seventeen siblings, which is not a measurement.
- **Sign coverage.** 29% of synaptic weight has no measured or predicted polarity.
  All four handling policies were tested, but none of them is the truth.
- **Muscle model.** Activation is read as contraction, and the neuromuscular
  transform is ours.

### 5.5 The negative result is conditional on the class of model

The claim is about a *graded, non-spiking, connectome-constrained* model with
proprioceptive feedback. It says nothing about whether some other dynamics over
the same connectome would work.

### 5.6 The simulation is not converged with respect to timestep

Identical configuration at different physics rates gives materially different
trajectories — `bend` 56.4°, 12.0° and 22.9° at 240, 120 and 60 Hz. **No
trajectory quoted here is reproducible at another timestep**, and the body spends
most of its time coiled and near its joint limits, which is where a physics
engine's contact handling contributes most.

The locomotion conclusion survives it: inter-joint phase reads within 2° of zero
at every rate, against the scripted control's +23.0°. Comparative results — one
condition against another, all taken at 240 Hz — remain internally comparable.
The absolute numbers were never predictions and should not be read as any.

Whether this is chaos or stiffness is not established. If chaos, no timestep is
"correct" and the honest form is distributions over seeds rather than single
trajectories.

## 6. What the published models supply instead

Three models reproduce *C. elegans* locomotion. None does so from measured
connectivity with measured parameters.

**Boyle, Berri & Cohen 2012** undulate with no central pattern generator, from
proprioceptive feedback — the same architecture as this project. The oscillator is
binary B-class motor neurons with hysteresis, switching at 0.75 and 0.25 on a unit
scale, with half-body stretch receptive fields and a dorsoventral gain asymmetry.
Each is a fitted modelling choice. Imported individually here, each was measured
to make this model *worse*.

**c302 (OpenWorm)** produces head-to-tail waves in muscle, from a hand-selected
39-neuron forward-locomotion network. The wave travels because the delays say so.
From `c302/parameters_C2.py`:

    AVBR_to_DB1_elec_syn_delay      0 ms
    AVBL_to_DB2_elec_syn_delay    250 ms
    AVBL_to_DB3_elec_syn_delay    500 ms
    ...
    AVBL_to_DB7_elec_syn_delay   1500 ms

with VB1–VB11 carrying the same structure and the AVB conductances a matching
ramp. A 250 ms per segment delay ladder from the command interneuron **is** a
travelling wave: the circuit is not computing the phase gradient, the parameter
file contains it. All 238 bioparameters there carry the authors' own
`certainty 0.1`.

Testing the *mechanism* rather than importing the gradient — one uniform delay,
swept 5 to 100 ms — produces no propagation (hypothesis 8). So c302's wave depends
on its delay being *graded*, which is the thing a model is supposed to explain.

**MetaWorm / BAAIWorm 2024** has multicompartment morphology, a closed
brain–body–environment loop, and reproduces zigzag chemotaxis. Its connection
weights and synaptic delays are "optimized … by iteratively adjusting model
parameters" to match electrophysiology.

The OpenWorm consortium's own assessment is direct: *"the level of detail that we
have incorporated to date is inadequate for biological research"*, with completion
of Hodgkin–Huxley channel parameter extraction named as a key remaining component.
That extraction is what §3.3 reports the result of.

So the field's position, stated from the other direction: **locomotion in a
*C. elegans* model is currently something you put in, not something you get out.**

## 7. Why publish a null

Three reasons.

**It is a strong null.** Nine hypotheses tested and rejected, each with a
mechanism and a measurement, with a positive control demonstrating the apparatus
detects the target behaviour when it is present — and with the one the others
pointed at, real intrinsic conductances, among the rejected.

**It localises the problem.** The failure is not the body, not the gap junctions,
not the missing inhibition, not the antagonist, not the receptive field, not the
dropped synapses, not conduction delay, and not the simplicity of the cells. It is
that a population receiving common input produces common output, and nothing
downstream of that can repair it.

**The obvious repairs are the published ones, and they are fitted.** Boyle et al.'s
model undulates, and each of the components responsible is a fitted modelling
choice rather than a measurement. Adopting them produces a gait that belongs to
the model, not to the connectome. That distinction is the whole point of building
from measured anatomy, and it is worth stating that the honest version of this
model does not yet crawl.

---

*Every measurement here is reproducible from the committed code. Commands,
parameter values and the full derivation — including every retraction along the
way — are in [`model_assumptions.md`](model_assumptions.md) §5C–§5AC. The
retractions are not incidental: five numbers in this document were first reported
with the opposite sign or magnitude, each from a measurement window shorter than
the transient it was measuring, and §5AC.2 records the general form of the
mistake.*
