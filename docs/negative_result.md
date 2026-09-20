# What a graded connectome model does not do

**A negative result.** An embodied *C. elegans* model built from the measured
whole-animal connectome, measured neuromuscular polarity, and the standard graded
leaky-integrator dynamics does not produce a travelling undulatory wave. This
document states that claim precisely, gives the measurements behind it, lists the
seven explanations tested and rejected, and is explicit about what it does *not*
establish.

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

### 5.3 The measurements themselves were wrong twice

Both are documented rather than quietly fixed, and both inflate confidence in the
wrong direction if not known:

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

Neither error changes the conclusion — re-measured correctly, the connectome
body's phase is −0.8° against the control's +23.0° — but every magnitude in
earlier analyses was on a compressed scale.

### 5.4 Missing biology that could matter

- **Extrasynaptic and neuropeptidergic signalling** (Ripoll-Sánchez et al. 2023)
  is entirely absent. It is not a small system in *C. elegans*.
- **Intrinsic conductances.** Every neuron here is graded and passive apart from
  its synaptic input. Real B-type cells may not be. The conductance-based models
  of Nicoletti et al. 2024 — now imported as data, not yet runnable (§5X) — are
  the direct test.
- **Sign coverage.** 29% of synaptic weight has no measured or predicted polarity.
  All four handling policies were tested, but none of them is the truth.
- **Muscle model.** Activation is read as contraction, and the neuromuscular
  transform is ours.

### 5.5 The negative result is conditional on the class of model

The claim is about a *graded, non-spiking, connectome-constrained* model with
proprioceptive feedback. It says nothing about whether some other dynamics over
the same connectome would work.

## 6. Why publish a null

Three reasons.

**It is a strong null.** Seven hypotheses tested and rejected, each with a
mechanism and a measurement, with a positive control demonstrating the apparatus
detects the target behaviour when it is present.

**It localises the problem.** The failure is not the body, not the gap junctions,
not the missing inhibition, not the antagonist, not the delay, and not the dropped
synapses. It is that a population receiving common input produces common output,
and nothing downstream of that can repair it.

**The obvious repairs are the published ones, and they are fitted.** Boyle et al.'s
model undulates, and each of the components responsible is a fitted modelling
choice rather than a measurement. Adopting them produces a gait that belongs to
the model, not to the connectome. That distinction is the whole point of building
from measured anatomy, and it is worth stating that the honest version of this
model does not yet crawl.

---

*Every measurement here is reproducible from the committed code. Commands,
parameter values and the full derivation — including retractions — are in
[`model_assumptions.md`](model_assumptions.md) §5C–§5X.*
