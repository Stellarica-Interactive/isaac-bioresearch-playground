# What is left

Status as of 2026-10-03. The project has a **working body and a non-working
nervous system**, which is a cleaner position than it has been in: when the worm
fails to crawl it is now possible to say the failure is neural without
qualification.

Measured on the scripted positive control, which is the same body, drag, solver
and torque the connectome gets:

| | model | real *C. elegans* |
|---|---|---|
| speed | 0.211 BL/s | 0.1–0.3 BL/s |
| slip | 0.283 | 0.1–0.3 |
| tail off its own track | 1.18 mm | about a body radius (1.09 mm) |
| body rotation per cycle | 12.8° | not established |
| heading held over 90 s | within 18.6° | — |
| self-intersection | none (`clear` ≥ +1.57 mm) | impossible |

Self-contact is implemented and off by default; with it on, a curl deep enough to
fold the body keeps **+0.79 mm** of clearance where it previously overlapped by
5.80 mm.

And on the connectome driving it: **0.717 BL in 60 s, inter-joint phase −6.1°,
`travel` +0.00, B-type adjacent correlation +1.00.** One signal, so no wave can
travel. That is the problem.

**Twelve hypotheses have now been tested and rejected.** The ten of
`negative_result.md` and §5Q–§5U, plus §5AK (a proprioceptive law that cannot
latch) and §5AL (a groove with memory). Both of the new ones are rejected with
their controls done, and both are documented with the measurement rather than the
intention.

---

## 1. The central problem

### 1.0 Start here: the loop is unstable on the wrong mode

The committed loop does not settle (`model_assumptions.md` §5AN.8). Its
equilibrium, a gentle 7° dorsal C, is **unstable**: a slow, tail-led, standing
mode at 0.065 Hz grows into a coiling cycle, a dorsal coil every ~22 s, each one
turning the body ~117° the same way. Over 200 s it turns three times and ends
6.6 mm from where it started. Isaac shows the same thing, a coil every 23–24 s.
The "curls again … tumbles through 115°" of §1.2 below is this cycle.

The forward-travelling mode — the gait — is damped at ζ 0.94 at that
equilibrium, and §5AN.9's pathway-by-pathway sensitivity says no reweighting of
the connectome within reason changes that: the smallest coordinated change
reaching zero damping is ~5 in log conductance, ×9.6 on one pathway. The coil,
by contrast, is set by the dorsal–ventral balance of the neuromuscular junctions
(DA/DB onto dorsal muscle stabilise it, VA/VB onto ventral feed it) and sits
within ×1.5 of its threshold.

What has been ruled out since, all at equilibria (§5AN.10–5AN.12): no measured
cell amplifies anywhere in 0.05–2 Hz (VB6 is a resistor across the band; RMD's
known bistability is the positive control); the threshold offset, which sets the
gain of every chemical synapse, never takes the gait mode below ζ 0.94; and
negative conductance in any group of motor neurons or muscles feeds the coil or
latches first. The gait mode lives in the body-wall muscle cells (66%); the
B-type neurons carry 3% of it.

**Why (§5AN.13):** with proprioception cut the loop is exactly reciprocal — it
cannot tell head from tail. Tonic proprioception is the only directional element,
and it is directional at 0 Hz (posterior/anterior 3.3, which builds the
tail-heavy coil) and almost reciprocal at the gait frequency (1.02).

**The first lever that works (§5AN.14): phasic proprioception.** Sensing the rate
of bending is exactly reciprocal at rest — the coil disappears — and strongly
tailward at 0.5 Hz (up to 8×). The gait mode's damping then falls monotonically
with gain, 0.907 → 0.573 from 20 to 320 mV/rad, its frequency crossing 0.5 Hz.

**And the loop crawls — backward (§5AN.15–5AN.17).** Under phasic
proprioception the connectome loop produces a sustained, regular, travelling
undulation that propels the body: in Isaac, 3.6 body lengths in 90 s, B-type
neighbours anti-correlated (−0.42). Every attractor found in a 24-run search
crawls tail first (−0.03 to −0.10 BL/s). The cause is phase, not wiring: the
rate sensor leads the bend by 90°, the five stages between sensing and muscle lag
by an amount growing with frequency, and the attractor runs at ~0.38 Hz, where
the lead still wins (headward, +6° per joint). At 0.5 Hz the loop already
propagates forward (−13° per joint, 8× stronger tailward) — it just does not
oscillate there.

**And it crawls forward (§5AN.18–5AN.20).** Weakening the gap junctions
shortens the attractor's wavelength and raises its frequency past the point
where propagation turns tailward. At phasic 400 mV/rad and `g_gap_ps` 5, Isaac
crawls forward at 1.9 s, `travel` +0.55, 6.36 body lengths in 120 s, heading
within 45° the whole run — predicted beforehand from the standalone loop. On
the way, a harness bug was found and fixed: `Loop` had 40 extra muscle cells
(§5AN.18); §5AN.19 re-ran everything on the runner's network and every earlier
conclusion stands.

1. **Controls in Isaac at 5 pS — DONE.** With `--no-proprioception` or
   `--lesion DB,VB` the body does not move at all (§5AN.20).
2. **The three assumptions it rests on** — rate-sensing proprioception, a
   400 mV/rad gain, gap junctions at 5 pS — are each a tenfold-or-more departure
   from committed values with nothing measured behind them. Which of them has
   biological support, and whether a smaller departure in all three does the
   same, is the question that decides whether this is a model of the worm or a
   tuned oscillator. Gap-junction strength per class (rather than one global
   value) is the obvious refinement: §5AN.9 found the muscle–muscle junctions
   matter most.
3. **What is still wrong with the gait**: wavelength about 0.9 body lengths
   against 0.65, speed a quarter of the positive control's, circling at 7–15 pS
   (a ventral bend bias of −2 to −3.5°), and the body overlapping itself by 1 mm
   with self-contact off.

### 1.0b The steady response is rank-one (§5AM)

`tools/analyse_modes.py` computes the linear transfer function of the whole
sensorimotor loop — the block of `(A - F)^-1` from B-type inputs to B-type
voltages, where `A` is the conductance matrix `neuron_models` already solves the
resting potentials with and `F` is the chemical feedback. Its answer (§5AM) is
more specific than anything twelve simulation experiments produced.

**The gap junctions make the loop rank-one.** Scaling `Ggap` alone:

| gap scale | first mode carries | condition number | sign changes in mode 1 |
|---|---|---|---|
| 1.00, committed | **99.38%** | 137 | **0** |
| 0.01 | 21.52% | 4.1 | 0 |
| 0.00 | **7.38%** | **1.3** | **12** |

`Ggap` is a graph Laplacian and a Laplacian's lowest mode is the constant vector,
so gap junctions make the network a diffusive medium — and diffusion smooths
rather than propagates. At 100 pS across 2883 junctions that one whole-body mode
carries 99.4% of the loop's response. `g_gap_ps` is **assumed**, not measured.

**And the thing to fix is not what it looked like.** A population `x_i = a_i s(t)`
has 100% of its variance in one component whatever the spatial weights are, so
"one signal" was never a statement about spatial uniformity. Removing the gap
junctions gives adjacent correlation −0.41 with a 12-sign-change mode and *still*
100% in one component: a **standing** wave. A travelling wave is
`sin(kx)cos(wt) - cos(kx)sin(wt)` — **two** temporal components in quadrature.

So the question is not how to decorrelate neighbouring cells. Removing gap
junctions already does that, and the body still does nothing (amplitude 0.05°).
**The question is what would give the loop a second temporal degree of freedom**,
which a first-order relaxation toward a single fixed point cannot have. That
points squarely at §5C.9's and §5C.5's remaining options — an intrinsic
oscillator, or a delayed branch — and away from the connectivity, the sign
policy, the proprioceptive law and the mechanics, which is where all twelve
rejected hypotheses were looking.

Three parameter explanations were proposed and refuted along the way, each
recorded before testing: homogeneous sigmoid gain, synaptic sign policy, and the
inhibitory reversal potential. `(A - F)^-1` has zero negative entries in 190,969,
so Perron–Frobenius makes the dominant mode a whole-body contraction as a
theorem, not a tuning outcome — which is why no parameter rescued it.

### 1.1 The symptom

The connectome's motor output is one signal. Every B-type motor neuron does the
same thing at the same time, so there is no phase gradient along the body and
nothing for the mechanics to turn into locomotion. Ten hypotheses have been
tested and rejected — gap junctions, B-type bistability, the D-class antagonist,
B-to-B inhibition, receptive-field offset, synaptic sign policy, drag anisotropy,
conduction delay, intrinsic conductances, and the body mechanics
(`negative_result.md`, `model_assumptions.md` §5Q–§5U, §5AF).

§5P establishes something worse than a failure to start: handed a scripted wave
and then released, **the loop destroys it**. So a symmetry-breaking mechanism
would not help; something in the model is actively synchronising.

### 1.2 The hand-over, re-run properly — DONE, and it is worse than §5P said

§5P's hand-over was measured on the inertial solver with a seed of `travel +0.12`,
and its gait statistic pooled the whole recording so the seeded phase kept
propping the number up afterwards. Both objections are now closed: the runner
gives the seeded phase its own report at hand-over and then clears the history
and the distance origin, so everything after is the loop on its own data.

Handed the §5AI.4 gait — 3.8 body lengths, `travel +0.92`, a clean head-to-tail
wave — the loop **destroyed it within one cycle**:

| | distance | amplitude | travel | period |
|---|---|---|---|---|
| the seed, at hand-over | 3.796 BL / 20 s | 16.17° | **+0.92** | 2.0 s |
| 2 s later | 0.025 BL | 4.55° | **+0.00** | *none* |
| 18 s later | 0.027 BL | **2.48°** | +0.00 | *none* |
| 60 s later, whole window | 0.250 BL | 7.19° | **+0.00** | *none* |

The wave dies in one cycle and never comes back. The body does not stay still —
over the full minute amplitude recovers to 7.19°, it curls again to an extent of
0.54 and tumbles through 115° — but `travel` is +0.00 and no period is resolvable
for the whole run. It resumes moving without ever propagating. Written up as
§5P.3.

**The sharpest open question in the project** comes out of it. While the scripted
wave drives the body the B-type population is *decorrelated* — adjacent
correlation **−0.40**, 57% shared variance — against **+1.00** and 100% under the
connectome. Proprioception does differentiate these cells when the body genuinely
undulates. So the circuit can represent a travelling wave and cannot hold one,
and the mechanism that collapses it the moment the external drive stops is not
identified. Finding it is the most direct route to the central problem, and it is
a different question from any of the ten hypotheses already rejected.

### 1.3 Missing biology that could matter

Listed in `negative_result.md` §5.4 and not yet tried:

* **Extrasynaptic and neuropeptidergic signalling** (Ripoll-Sánchez et al. 2023).
  The model has only wired chemical and electrical synapses. The neuropeptide
  network is a separate, denser graph.
* **Directional proprioception — partly answered, and the remaining part is
  narrow.** `Proprioception.asymmetric` already implements Boyle's directional
  stretch/compression gains and latches the model. §5AK then tested the other
  axis — responding to curvature *rate* rather than curvature, which is §5C.5's
  option 3 — and rejected it: the latch goes, the wave does not come, and at a
  compensating drive of 114 mV/rad the bend returns to 18.9° with `travel` still
  +0.00. What is left untried is **spatial** directionality: a receptor that
  distinguishes a bend arriving from the anterior from one arriving posteriorly,
  which is neither of the two things now tested.
* **AVB–B-class gap-junction coupling** as a travelling-wave mechanism rather
  than as ordinary resistive coupling.
* **Muscle-level dynamics.** Muscles are a linear map from activation to joint
  torque, with no calcium dynamics and no cell-level coupling.

---

## 2. Loose ends from the mechanics work

All measured, all open.

* **Self-collision — DONE.** `QuasiStaticBody(self_contact=True)` adds stiff drag
  along the line between two non-neighbouring segments that have closed to
  touching, on approach only, and `--self-collision` now reaches the quasi-static
  solver instead of being silently ignored. Worst clearance on a realistic curl
  goes from **−5.80 mm to +0.79 mm**. The band width is bounded by a measurement
  rather than chosen: the scripted gait's closest self-approach is 1.57 mm, so a
  wider band would engage during ordinary crawling — at 1 mm the gait's speed is
  identical to four decimals with contact on and off. Three tests cover it.
* **The 60° joint limit permits shapes no physical body can adopt.** Found while
  testing the above: on a hard curl the contact model barely helps, because at 60°
  per joint the body wraps **3.8 times** where a closed circle of it needs only
  **15.7°** per joint. Self-intersection is geometrically forced at that limit and
  no contact model can undo it. The limit was chosen to match the Isaac
  articulation (`worm/isaac_limits.py`), not from anatomy. An anatomical ceiling
  would make the contact model's job tractable everywhere rather than only at
  realistic amplitudes.
* **§5AF.1's Isaac timestep row** needs re-taking at default settings on the
  fixed render (240 Hz and 60 Hz). It was withdrawn under §5AJ and the standalone
  table already supports the claim, so this is confirmation rather than evidence.
* **§5Q–§5U distance columns are still inertial-solver upper bounds.** Never
  re-taken on the quasi-static solver. The phase columns, which is what the nine
  rejected hypotheses rest on, are unaffected.
* **The nonlinear drag law costs 0.65× real time** at 240 Hz — 26 linear solves
  per step after warm-starting, against one for the linear law. Acceptable as an
  opt-in flag; it would need optimising before becoming a default.
* **Why tapering the wave costs track-following is unexplained.** The obvious
  hypothesis — that it denies the head the excursion needed to cut a groove — was
  refuted: head excursion is preserved at a taper width of 0.12 and the track
  degrades anyway (§5AI.3). A gentler envelope over fewer segments has not been
  tried.
* **Decide whether to keep `--wave-taper` at all.** It is a documented negative
  result that only reduces rotation by reducing bend. Keeping it costs a flag and
  a test; removing it loses the record of having tried.

### 2.1 Parameters that are fitted rather than measured

Three, and they should be read together because they are the project's weakest
numbers:

| | value | status |
|---|---|---|
| drag exponent | 0.6 | ASSUMED — Rabets et al. 2014 established agar is nonlinear and groove-dominated, published no exponent |
| yield force | 2e-5 N | ASSUMED — no published value for a worm on agar |
| drag ratio | 40 | above the measured 3–10, but the conventional neuromechanical value (Niebur & Erdös 1991) |
| seed wavelength | 0.55 BL | **shorter than the animal's 0.65** — straightness bought outside the real gait |

The first two are two free parameters fitted against two targets, which is no
constraint at all (§5AH.3).

**A track-memory term was the proposed replacement, and the measurement says
there is nothing for it to do** (§5AL). In its final form -- a restoring force
along the body's own normal, acting only outside the channel -- it is correct and
never engages: the gait's lateral offset is **0.701 mm** against a channel
half-width of one body radius, so **the body already stays inside its own
groove**, with 0.391 mm to spare at the tightest point. The three parameters
below are already producing track-conforming motion. Four earlier formulations
immobilised the animal and each ruled out a different explanation.
The groove is implemented, remembers the head's path, and immobilises the animal
— 0.241 BL/s to between 0.008 and 0.066 — for a reason that is measured rather
than guessed: the track's tangent sits 25–34° off the body's own, so a drag
tensor built on the track normal applies 7 to 12 times the body's tangential
resistance against the direction it is travelling. Gripping harder *across the
body* costs nothing (the anisotropy sweep gives 0.237 to 0.239 BL/s from ratio 20
to 160); gripping across *the track* costs everything.

### 2.2 The constraint formulation — the top mechanical item

Two independent lines now reach the same architectural answer, which is why this
is worth doing properly rather than approximating a third time.

**Self-contact** (§5AF.7): two segments that have closed to touching must not
approach further. Implemented as a stiff drag penalty, which works at realistic
amplitudes and cannot work at the joint limit.

**The groove** (§5AL): a segment lying in its own channel may move along it and
not out of it. Implemented as a drag tensor, which cannot express "zero this way,
free that way" without resisting whatever else projects onto it.

Both are **constraints, not forces** — and the constraint formulation has now
been tried, for the groove, and it fails. §5AL: a hard `v·n = 0` at a dozen
segments stops the animal dead, because the condition demands the body lie
exactly along its remembered track and the 0.74 mm by which it never does makes
the rows conflict. A dead zone did not rescue it.

**That failure does not transfer to self-contact, and the difference is the
point.** The groove asks a body to conform to a *remembered* path, so any
accumulated mismatch is fatal. Self-contact is a condition between two *current*
positions — no memory, nothing to conform to, and the standard unilateral
complementarity problem. It remains untried and is still worth doing: the penalty
version works at realistic amplitudes and cannot work at the joint limit, where
§5AF.6's active set already handles exactly this shape of condition.

So the item is narrower than it was: **self-contact as a row in the active set**,
and the groove left behind `groove=False` as correct-but-idle (§5AL.2).

---

## 3. Neural model gaps

* **AVAL misses its published resting potential by 14 mV.** A strict xfail,
  undiagnosed. To be chased through the repository's other parameter sets rather
  than by tuning (§5AD).
* **`cadiff.mod` is adapted from a 2002 Purkinje cell model** (ModelDB 48332), so
  four of the seven Nicoletti 2024 cells carry non-worm calcium dynamics.
  Documented; not resolved.
* **Validation is against one scalar per cell.** `VB6_simulated_IV_SS_LEAK.txt`
  and the other published current–voltage curves would test the dynamics rather
  than a single resting potential.
* **29% of the animal is unsigned.** 1846 of 4879 chemical connections have no
  known sign, carrying 8239 of 28113 synaptic weight.
* **Decision 6.1 is not implemented.** Every experiment whose result could depend
  on polarity is supposed to be run both with and without the Fenyves 2020
  overlay, and the comparison reported. That has not been set up.

---

## 4. Behaviour, once there is a gait

Not blocked on anything but the gait itself.

* **Touch** works — proximity-based, the probe deliberately carries no collider
  (§5D). Scripted pokes with sham controls are in the runner.
* **Chemotaxis** has a module (`worm/body/chemotaxis.py`) and sensors wired, but
  no behavioural result: gradient climbing needs directed locomotion.
* **Escape response (W4)** not started. The classic harsh-touch reversal is the
  obvious first behaviour because it has a known circuit and an unambiguous
  readout.
* **Part II, *Drosophila*** — `fly/README.md` only, by design. The tree gets
  built when a FlyWire or BANC dataset is actually being parsed, not before.

---

## 5. Process debt

**The Isaac layer is structurally untestable here.** pytest runs in the project
venv, which has no Isaac Sim, so nothing in `worm/tests` can construct an
articulation. 603 tests passed on a kinematic drive that was writing the head's
pose onto segment 11 and would have gone on passing (§5AJ.6).

What stands in for tests:

* `tools/check_kinematic_drive.py` — compares the articulation's segment
  positions against the solver's, segment by segment. **Run it by hand after any
  change to the kinematic drive or to `build_scene`.** It found 59 mm of error
  and now reports 0.000 mm.
* `tools/check_root_pose.py` — confirms a written root orientation takes effect.
* `tools/diagnose_gait.py` — track following and the rotation decomposition.
* `tools/diagnose_slip.py`, `tools/diagnose_yaw.py` — slip, creep, yaw.

Six defects in one session were found by a person watching the viewport and none
by the suite. That is not an argument for fewer tests; it is an argument against
mistaking a green suite for a working simulation, and for the standing rule that
Isaac runs are watched rather than run headless.

### 5.1 A recurring failure mode worth naming

Nine statistics in this project have been built on signals that did not have the
structure the statistic assumed, and the pattern is always the same: **a
conclusion drawn from a table that did not contain the relevant column.** The
cases are recorded in §5R, §5AC, §5AE, §5AG.2, §5AI.1, §5AI.3, §5AJ.3 and
`negative_result.md` §5.3.

The cheap guards that actually caught things: measure the settling time before
choosing a window; sample faster than the period you are studying; check where a
number comes from before quoting or withdrawing it; and prefer two genuinely
independent measures over one computed twice.
