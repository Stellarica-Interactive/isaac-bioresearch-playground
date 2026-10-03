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

---

## 1. The central problem

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

### 1.1 The hand-over, re-run properly — DONE, and it is worse than §5P said

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

### 1.2 Missing biology that could matter

Listed in `negative_result.md` §5.4 and not yet tried:

* **Extrasynaptic and neuropeptidergic signalling** (Ripoll-Sánchez et al. 2023).
  The model has only wired chemical and electrical synapses. The neuropeptide
  network is a separate, denser graph.
* **Directional proprioception.** Real B-class stretch reception is thought to be
  directional along the body; ours reads local curvature symmetrically. §5T
  rejected a simple receptive-field *offset*, which is not the same as rejecting
  directionality.
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
constraint at all (§5AH.3). The honest next step is a measurement, not a better
fit: the groove is a *history-dependent* lateral constraint and nothing in the
model has any memory of where the body has been. A track-memory drag term is the
mechanism these three parameters are standing in for.

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
