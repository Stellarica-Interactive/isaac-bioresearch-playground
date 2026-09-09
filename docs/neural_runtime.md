# The neural runtime

What one timestep does, in equations and in code, and what it rests on.

Companion reading: [biology.md](biology.md) for why *C. elegans* neurons are
modelled as graded rather than spiking, and
[model_assumptions.md](model_assumptions.md) for the register of everything here
that is not measured.

---

## 1. The short version

Every cell holds **two numbers**:

| Symbol | Meaning | Unit |
|---|---|---|
| `V` | membrane potential | mV |
| `s` | synaptic activation — how hard this cell is driving its targets | dimensionless, 0 to 1 |

A timestep computes four currents into each cell, converts the total into a rate
of change of voltage, updates `s` toward a sigmoid of the cell's own voltage, and
advances both by `dt`.

There is no spike, no threshold crossing, no event. Everything is continuous —
which is the point: this is an *analog* nervous system.

## 2. The equations

```
C dV_i/dt = -G_leak (V_i - E_leak)                    leak
            + Σ_j Ggap[i,j] (V_j - V_i)               gap junctions
            + Σ_j Gsyn[i,j] · s_j · (E[i,j] - V_i)    chemical synapses
            + I_ext,i                                 injected (sensory) current

  ds_i/dt = a_r · φ(V_i) · (1 - s_i) - a_d · s_i
     φ(V) = 1 / (1 + exp(-β (V - Vth_i)))
```

Reading each term as a sentence:

- **Leak** pulls the cell toward `E_leak`. Alone, it relaxes exponentially with
  time constant `τ = C / G_leak` = 100 ms.
- **Gap junctions** pull the cell toward each electrically coupled neighbour,
  proportionally to the voltage difference. Ohm's law, nothing more.
- **Chemical synapses** pull the cell toward each incoming synapse's reversal
  potential `E[i,j]`, scaled by how active the presynaptic cell is (`s_j`). If
  `E[i,j]` is above the cell's voltage the synapse depolarizes it; below, it
  hyperpolarizes.
- **`I_ext`** is where sensation enters. A "smell" is a current injected into a
  chemosensory neuron. Nothing else in the model knows about the world.
- **`ds/dt`** says activation chases a sigmoid of the cell's own voltage, rising
  at rate `a_r` and decaying at `a_d`. This is the graded analogue of "firing".

### Two details that are not cosmetic

**`E[i,j]`, not `E_j`.** Reversal potential is stored per *connection*, not per
presynaptic neuron. The receptor that decides whether a synapse excites or
inhibits sits on the **postsynaptic** cell, so one neuron can excite one target
and inhibit another. Most published implementations of this model collapse sign
to one value per presynaptic neuron and cannot express that; ours keeps the full
matrix, because the polarity data is per connection.

**`Vth_i` is solved for, not chosen.** See §4.

## 3. Units

Chosen so that the numbers stay readable and the arithmetic stays consistent:

| Quantity | Unit |
|---|---|
| voltage | mV |
| time | ms |
| capacitance | pF |
| conductance | nS *(config files use pS, converted on load)* |
| current | pA |

The consistency check: `C dV/dt` in pF·mV/ms is pA, and `G·V` in nS·mV is also pA.

## 4. Where the resting state comes from

`Vth_i` — the midpoint of each neuron's release sigmoid — is not a free parameter.
It is set to that neuron's **resting potential**, so every neuron sits in the
responsive middle of its range at rest.

Finding those resting potentials is a single linear solve rather than something to
integrate towards. Pin every sigmoid at `φ = 0.5`; then synaptic activation settles
at a constant

```
s* = 0.5·a_r / (0.5·a_r + a_d) = 1/11 ≈ 0.0909
```

and with `s` fixed the voltage equation becomes linear in `V`:

```
A V = b
A = diag(G_leak + rowsum(Ggap) + s*·rowsum(Gsyn)) - Ggap
b = G_leak·E_leak + s*·rowsum(Gsyn ∘ E) + I_ext
```

One `numpy.linalg.solve` and every neuron has a resting potential.

> **This is an assumption, not a measurement.** It asserts that a resting worm's
> neurons sit centred in their dynamic range. Nobody has measured resting
> potentials across the *C. elegans* nervous system. It is also self-referential in
> a way worth noticing: the resting state is defined by the model, then used to
> parameterise the model.
>
> A side effect: a resting synapse is already passing about 9% of its maximum
> conductance, because `s* = 1/11`. Synapses here are never fully off.

## 5. One timestep, in code

From `common/neural/neuron_models.py`, unabridged:

```python
def currents(self, v, s, i_ext, network):
    p = self.params
    i_leak = -p.g_leak_ns * (v - p.e_leak_mv)
    i_gap  = network.g_gap @ v - v * network.gap_row_sum
    i_syn  = network.g_syn_e_rev @ s - v * (network.g_syn @ s)
    return {"leak": i_leak, "gap": i_gap, "chemical": i_syn, "external": i_ext}

def derivatives(self, state, i_ext, network):
    p = self.params
    v, s = state[0], state[1]
    terms = self.currents(v, s, i_ext, network)
    dv = (terms["leak"] + terms["gap"] + terms["chemical"] + terms["external"]) / p.c_m_pf
    phi = sigmoid(v, self._thresholds(network), p.beta_per_mv)
    ds = p.a_r_per_ms * phi * (1.0 - s) - p.a_d_per_ms * s
    return np.vstack([dv, ds])
```

Three matrix-vector products and some element-wise arithmetic. That is the entire
nervous system.

The two rewrites in `currents` are algebraically identical to the obvious form but
avoid allocating an `n×n` temporary on every call, of which a run needs millions:

```
Σ_j Gsyn[i,j]·s_j·E[i,j]  ==  (Gsyn ∘ E) @ s     ← precomputed as g_syn_e_rev
Σ_j Gsyn[i,j]·s_j         ==  Gsyn @ s
```

Skipping that is roughly a 30× slowdown.

## 6. Advancing the state

Three integrators, `common/neural/runtime.py`:

| | What it does | When to use it |
|---|---|---|
| `EULER` | `state += dt · derivatives(state)` | Simple enough to check by hand. Teaching and comparison only. |
| `RK4` | Classical fourth-order Runge–Kutta, four derivative evaluations | Reference results. Error falls as `dt⁴`. |
| `EXPONENTIAL` | Solves each cell's own linear leak-plus-conductance term **exactly** over the step, treating coupling explicitly | **Default.** Stable where the explicit methods are not. |

### Why the default is not RK4

This network is **stiff**. A neuron's time constant is `C / G_total`, and
`G_total` includes every gap junction it makes. An isolated neuron has τ = 100 ms;
the busiest cells in Cook 2019 come out near **0.008 ms**, more than ten thousand
times faster. The whole network then has to be integrated at the pace of its
fastest cell.

Measured on `cook_2019_herm`, 50 ms simulated, against an RK4 reference at
dt = 0.002 ms (`python tools/benchmark_runtime.py accuracy`):

| integrator | dt (ms) | max error (mV) | wall (s) |
|---|---:|---:|---:|
| rk4 | 0.005 | 0.000000 | 8.31 |
| rk4 | 0.020 | **diverged** | — |
| exponential | 0.020 | 0.001232 | 1.16 |
| exponential | 0.100 | 0.013426 | 0.20 |
| exponential | 0.500 | 0.211947 | 0.03 |

At dt = 0.1 ms the exponential integrator is **40× faster than RK4** for 0.013 mV
of error. That is the default. Anything that matters should be re-checked with
`rk4` at dt ≤ 0.005 ms.

`stability_report()` gives the limit for the current configuration, and `step()`
raises with an actionable message rather than quietly returning `NaN`.

## 7. Sensory in, motor out

```python
from worm.neural.config import build_runtime
from common.neural.synapses import UnknownSignPolicy

runtime, report = build_runtime("cook_2019_herm", unknown_sign=UnknownSignPolicy.EXCLUDE)
print(report.summary())          # what this run rests on

runtime.inject_many({"ASHL": 20.0, "ASHR": 20.0})   # pA into the nociceptor pair
runtime.run(duration_ms=500)

runtime.voltage("AVAL")                # membrane potential
runtime.activation("AVAL")             # drive onto its targets — the motor readout
runtime.current_breakdown("AVAL")      # why it moved: leak / gap / chemical / external
```

Sensation is a current into a sensory neuron. Motor output is the activation of
motor neurons. Nothing between them is scripted.

`current_breakdown` matters for the project's aims: it answers *"the worm moved
because this pathway carried current"* with numbers rather than assertion.

## 8. Checkpointing

`save_checkpoint` writes every state variable of every cell, the clock, the step
count and the standing input, plus a **fingerprint of the network** — so a
checkpoint cannot be restored into different wiring and quietly produce nonsense.
Restoring gives a bit-identical continuation, which the test suite asserts.

A lesion experiment deliberately wants to restore a neural state into an altered
network; that requires `allow_different_network=True`, so it cannot happen by
accident.

> This is continuity of a **numerical simulation**. It is not preservation of a
> mind. The connectome it runs on contains no memories to preserve.

## 9. Determinism

Same connectome, parameters, integrator, timestep and inputs → bit-identical
results, across runs and processes. No randomness is used anywhere. If noise is
added later it must come from a seeded generator stored in the checkpoint.
Asserted by `TestDeterminism`.

## 10. What this rests on — read before believing any output

**Nine of the ten biophysical parameters are assumptions.** Only whole-cell
capacitance has measured support, and even that is one figure applied to every
neuron. In this model's lineage the **connectivity is the only measured
quantity**.

So the obvious question is: how much does an answer depend on the connectome, and
how much on our guesses? `tools/benchmark_runtime.py sweep` measures exactly that.
Varying one assumption at a time, and asking how far AVAL moves when the ASH
nociceptor pair is driven:

```
baseline (exclude unsigned, linear scaling): +0.4137 mV

unknown-sign policy = excitatory       +0.1414    -65.8%
weight scaling = binary                +2.0154   +387.2%
weight scaling = sqrt                  +1.0620   +156.7%
a_d_per_ms x 2.0                       +0.7923    +91.5%
g_syn_ps  x 0.5                        +0.7615    +84.1%
a_r_per_ms x 0.5                       +0.7346    +77.6%
...
g_leak_ps x 2.0                        +0.4075     -1.5%
```

**The magnitude of the response spans a factor of fourteen across defensible
choices.** It is not a property of the worm; it is a property of our parameters.

What can be said more robustly is the *sign and ordering* of effects — which
pathways carry current, which cells move first, which are unaffected — because
those follow from connectivity, which is the part that was actually measured.
Structure the experiments accordingly, and report the sweep alongside the result.

Two further limits, both in [model_assumptions.md](model_assumptions.md):

- **Roughly half of chemical connections have no predicted sign**, and the
  unknown-sign policy is a required argument precisely because there is no
  defensible default.
- **No neuromuscular junction has a sign at all.** Fenyves et al. covers
  interneuronal connections only, so the synapses that actually drive muscle are
  unsigned. That has to be resolved before locomotion (W2).
