"""Conductance-based single neurons, for the few cells where one was measured.

    model = ConductanceModel.load("rmd_nicoletti2019")
    state = model.initial_state(n_cells=6)
    state = model.run(state, duration_ms=400.0, i_ext_pa=10.0)

Why this exists
---------------

The graded leaky integrator in :mod:`common.neural.neuron_models` has no
intrinsic currents. It has one first-order membrane equation and one first-order
synaptic activation per cell, and no voltage-gated conductances at all. As
``docs/model_assumptions.md`` §5C.9 establishes empirically, that means it has no
limit cycle **anywhere**: driven at any cell, at any amplitude tried, every
trajectory relaxes to a fixed point and stays there. Sustained variation after
the transient measured 4e-5 mV, which is floating-point residue.

So the connectome cannot supply a rhythm under that neuron model, and no amount
of rewiring changes it. The missing ingredient is intrinsic dynamics.

Scope: only where it was actually measured
------------------------------------------

*C. elegans* neurons are mostly non-spiking and graded, but not entirely, and
``docs/model_assumptions.md`` §5 already records the exceptions. **RMD** is the
one that matters for locomotion: it is a head motor neuron that innervates body
wall muscle directly (§5C.8), and it has plateau potentials -- a slow, bistable,
regenerative depolarisation, which is exactly the nonlinearity a relaxation model
lacks.

This module deliberately does **not** extend conductance-based dynamics to
neurons where nobody has measured them. Nicoletti et al. characterised two cells
in 2019 (AWC^on and RMD) and six more in 2024; those are the cells that may have
this model, and every other cell stays graded. Extrapolating a channel complement
across 300 neurons would be the single largest silent assumption in the project,
which is why §6.2 rejected Option C wholesale.

Where the numbers come from
---------------------------

Parameters are **not typed in here.** They are machine-extracted from the
authors' own XPPAUT source by ``tools/import_neuron_model.py`` into
``worm/neural/models/*.toml``, each carrying the line it came from. This module
implements the equations; the constants are never retyped. See that tool for why.

The equations below are transcribed from the same source rather than from the
paper, deliberately: Nicoletti et al. 2019 carries a published Correction
(doi:10.1371/journal.pone.0256930) reporting a missing plus sign in a Boltzmann
denominator and errors in twelve supplementary equations. The code the authors
ran does not contain those errors.

Units are mV, ms, pF, nS, pA -- already this repository's convention, and
mutually consistent, since nS x mV = pA.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Any

import numpy as np

MODELS_PACKAGE = "worm.neural.models"

#: Integration timestep, ms. **Check a trajectory, never an endpoint.**
#:
#: The authors ran this model with a stiff solver at dt = 0.01 ms. Exponential
#: gating lets us go considerably coarser, but measuring by how much has a trap in
#: it: rest and plateau are both *equilibria*, so their final voltages are
#: identical to four decimal places from dt = 0.01 all the way to dt = 1.0 ms.
#: Comparing final voltages therefore "shows" that every timestep is fine, while
#: saying nothing at all about the path taken to get there.
#:
#: Comparing the trajectory instead, peak error against dt = 0.01 over a 120 ms
#: current step:
#:
#:     dt 0.02 -> 0.07 mV     dt 0.50 -> 3.2 mV
#:     dt 0.05 -> 0.29 mV     dt 1.00 -> 5.5 mV
#:     dt 0.10 -> 0.65 mV     dt 2.00 -> 12.9 mV
#:
#: The plateau itself survives even a coarse step -- it is a robust feature, not a
#: fragile one -- so the cost of a large dt is accuracy in the transient rather
#: than loss of the behaviour. 0.05 ms keeps that under a third of a millivolt at
#: a price worth paying.
DEFAULT_DT_MS = 0.05

#: Order of the state array's first axis. Voltage first so that ``state[0]`` is
#: the membrane potential for every model, matching the graded runtime.
RMD_STATE_NAMES: tuple[str, ...] = (
    "v",
    "m_shal",
    "hf_shal",
    "hs_shal",
    "m_shak",
    "h_shak",
    "m1_egl36",
    "m2_egl36",
    "m3_egl36",
    "m_kir",
    "m_unc2",
    "h_unc2",
    "m_egl19",
    "hs_egl19",
    "m_cca1",
    "h_cca1",
    "mbk",
    "mbk2",
    "mslo1",
    "mslo2",
    "m_sk",
    "ca_intra1",
)


@lru_cache(maxsize=4)
def load_parameters(model_id: str) -> dict[str, float]:
    """Numeric constants for a model, as extracted from its published source."""
    text = (files(MODELS_PACKAGE) / f"{model_id}.toml").read_text(encoding="utf-8")
    data = tomllib.loads(text)
    return {k: float(v) for k, v in data["parameters"].items()}


@lru_cache(maxsize=4)
def model_provenance(model_id: str) -> dict[str, Any]:
    text = (files(MODELS_PACKAGE) / f"{model_id}.toml").read_text(encoding="utf-8")
    return dict(tomllib.loads(text)["source"])


def _boltzmann(v: np.ndarray, midpoint: float, slope: float) -> np.ndarray:
    """``1 / (1 + exp(-(v - midpoint) / slope))``.

    The plus sign in the denominator is the one the published Correction restores;
    the authors' source always had it.
    """
    return 1.0 / (1.0 + np.exp(-(v - midpoint) / slope))


@dataclass(frozen=True)
class ConductanceModel:
    """A Hodgkin-Huxley style neuron, vectorised over however many cells use it.

    Stateless: all state lives in the array passed to :meth:`derivatives` and
    :meth:`step`, so one instance serves every RMD cell in the network at once.
    """

    model_id: str
    p: dict[str, float]
    state_names: tuple[str, ...] = RMD_STATE_NAMES

    @classmethod
    def load(cls, model_id: str = "rmd_nicoletti2019") -> ConductanceModel:
        return cls(model_id=model_id, p=load_parameters(model_id))

    @property
    def n_state(self) -> int:
        return len(self.state_names)

    def index(self, name: str) -> int:
        return self.state_names.index(name)

    # -- initial conditions ------------------------------------------------

    def initial_state(self, n_cells: int = 1) -> np.ndarray:
        """The published initial values, one column per cell.

        Taken from the source's ``init`` statements, not chosen by us. The resting
        state that matters is the one the model relaxes to, which
        :meth:`settle` finds.
        """
        text = (files(MODELS_PACKAGE) / f"{self.model_id}.toml").read_text(encoding="utf-8")
        initial = tomllib.loads(text)["initial"]
        state = np.zeros((self.n_state, n_cells), dtype=np.float64)
        for i, name in enumerate(self.state_names):
            state[i] = float(initial.get(name, 0.0))
        return state

    def settle(self, state: np.ndarray, *, duration_ms: float = 2000.0) -> np.ndarray:
        """Relax to rest with no input. The model's own resting state, not ours."""
        return self.run(state, duration_ms=duration_ms, i_ext_pa=0.0)

    # -- the model ---------------------------------------------------------

    def currents(self, state: np.ndarray) -> dict[str, np.ndarray]:
        """Every ionic current, in pA. Positive is outward, as in the source."""
        p = self.p
        s = {name: state[i] for i, name in enumerate(self.state_names)}
        v = s["v"]

        i_shal = (
            p["gshal"]
            * s["m_shal"] ** 3
            * (0.7 * s["hf_shal"] + 0.3 * s["hs_shal"])
            * (v - p["ek"])
        )
        i_shak = p["gshak"] * s["m_shak"] * s["h_shak"] * (v - p["ek"])
        i_egl36 = (
            p["gegl36"]
            * (p["a1"] * s["m1_egl36"] + p["a2"] * s["m2_egl36"] + p["a3"] * s["m3_egl36"])
            * (v - p["ek"])
        )
        i_kir = p["gkir"] * s["m_kir"] * (v - p["ek"])
        i_unc2 = p["gunc2"] * s["m_unc2"] * s["h_unc2"] * (v - p["eca"])
        i_egl19 = p["gegl19"] * s["m_egl19"] * s["hs_egl19"] * (v - p["eca"])
        i_cca1 = p["gcca1"] * s["m_cca1"] ** 2 * s["h_cca1"] * (v - p["eca"])
        i_bk = p["gbk"] * s["mbk"] * s["h_unc2"] * (v - p["ek"])
        i_bk2 = p["gbk2"] * s["mbk2"] * s["hs_egl19"] * (v - p["ek"])
        i_slo1 = p["gslo1"] * s["mslo1"] * s["hs_egl19"] * (v - p["ek"])
        i_slo2 = p["gslo2"] * s["mslo2"] * s["h_unc2"] * (v - p["ek"])
        i_sk = p["gca"] * s["m_sk"] * (v - p["ek"])
        i_leak = p["gleak"] * (v - p["eleak"])
        i_nca = p["gnca"] * (v - p["ena"])

        return {
            "shal": i_shal,
            "shak": i_shak,
            "egl36": i_egl36,
            "kir": i_kir,
            "unc2": i_unc2,
            "egl19": i_egl19,
            "cca1": i_cca1,
            "bk": i_bk,
            "bk2": i_bk2,
            "slo1": i_slo1,
            "slo2": i_slo2,
            "sk": i_sk,
            "leak": i_leak,
            "nca": i_nca,
        }

    def _gates(self, state: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        """Steady-state values and time constants for every gating variable."""
        p = self.p
        s = {name: state[i] for i, name in enumerate(self.state_names)}
        v = s["v"]

        inf: dict[str, np.ndarray] = {}
        tau: dict[str, np.ndarray] = {}

        # --- SHL-1 (Kv4) --------------------------------------------------
        inf["m_shal"] = _boltzmann(v, p["vashal"] - p["shalsfhit"], p["kashal"])
        tau["m_shal"] = (
            p["ptmshal1"]
            / (
                np.exp(-(v - p["ptmshal2"]) / p["ptmshal3"])
                + np.exp((v - p["ptmshal4"]) / p["ptmshal5"])
            )
            + p["ptmshal6"]
        ) * p["cashal"]
        h_shal = 1.0 / (1.0 + np.exp((v - p["vishal"] + p["shalsfhit"]) / p["kishal"]))
        inf["hf_shal"] = h_shal
        inf["hs_shal"] = h_shal
        tau["hf_shal"] = (
            p["pthfshal1"] / (1.0 + np.exp((v - p["pthfshal2"]) / p["pthfshal3"])) + p["pthfshal4"]
        ) * p["cshal"]
        tau["hs_shal"] = (
            p["pthsshal1"] / (1.0 + np.exp((v - p["pthsshal2"]) / p["pthsshal3"])) + p["pthsshal4"]
        ) * p["cshal"]

        # --- SHK-1 (Kv1) --------------------------------------------------
        shift = p["shiftV05"]
        inf["m_shak"] = _boltzmann(v, p["vashak"] - shift, p["kashak"])
        tau["m_shak"] = (
            p["ptmshak1"]
            / (
                np.exp(-(v - (p["ptmshak2"] + shift)) / p["ptmshak4"])
                + np.exp((v - (p["ptmshak2"] + shift)) / p["ptmshak3"])
            )
            + p["ptmshak5"]
        )
        inf["h_shak"] = 1.0 / (1.0 + np.exp((v - p["vishak"] + shift) / p["kishak"]))
        tau["h_shak"] = np.full_like(v, p["pthshak"])

        # --- EGL-36 (Kv3), three kinetic components sharing one activation --
        m_egl36 = _boltzmann(v, p["va_egl36"], p["ka_egl36"])
        for i, key in enumerate(("m1_egl36", "m2_egl36", "m3_egl36"), start=1):
            inf[key] = m_egl36
            tau[key] = np.full_like(v, p[f"t{i}_egl36"])

        # --- IRK (inward rectifier) ---------------------------------------
        # The source writes `(v - va_kir + 30)`, an inactivation-style Boltzmann
        # with a hard-coded +30 mV shift; reproduced as written.
        inf["m_kir"] = 1.0 / (1.0 + np.exp((v - p["va_kir"] + 30.0) / p["ka_kir"]))
        tau["m_kir"] = (
            p["p1tmkir"]
            / (
                np.exp(-(v - p["p2tmkir"]) / p["p3tmkir"])
                + np.exp((v - p["p4tmkir"]) / p["p5tmkir"])
            )
            + p["p6tmkir"]
        )

        # --- UNC-2 (CaV2) -------------------------------------------------
        inf["m_unc2"] = _boltzmann(v, p["va_unc2"] - p["stm2"], p["ka_unc2"])
        tau["m_unc2"] = (
            p["p1tmunc2"]
            / (
                np.exp(-(v - p["p2tmunc2"] + p["shiftmunc2"]) / p["p3tmunc2"])
                + np.exp((v - p["p2tmunc2"] + p["shiftmunc2"]) / p["p4tmunc2"])
            )
            + p["p5tmunc2"]
        ) * p["constmunc2"]
        inf["h_unc2"] = 1.0 / (1.0 + np.exp((v - p["vi_unc2"] + p["sth2"]) / p["ki_unc2"]))
        tau["h_unc2"] = (
            p["p1thunc2"] / (1.0 + np.exp((v - p["p2thunc2"] + p["shifthunc2"]) / p["p3thunc2"]))
            + p["p4thunc2"] / (1.0 + np.exp(-(v - p["p5thunc2"] + p["shifthunc2"]) / p["p6thunc2"]))
        ) * p["consthunc2"]

        # --- EGL-19 (CaV1) ------------------------------------------------
        inf["m_egl19"] = _boltzmann(v, p["va_egl19"] - p["stm19"], p["ka_egl19"])
        tau["m_egl19"] = (
            p["pdg1"]
            + p["pdg2"] * np.exp(-((v - p["pdg3"] + p["stau19"]) ** 2) / p["pdg4"] ** 2)
            + p["pdg5"] * np.exp(-((v - p["pdg6"] + p["stau19"]) ** 2) / p["pdg7"] ** 2)
        )
        inf["hs_egl19"] = (
            p["p1hegl19"] / (1.0 + np.exp(-(v - p["p2hegl19"] + p["sth19"]) / p["p3hegl19"]))
            + p["p4hegl19"]
        ) * (
            p["p5hegl19"] / (1.0 + np.exp((v - p["p6hegl19"] + p["sth19"]) / p["p7hegl19"]))
            + p["p8hegl19"]
        )
        tau["hs_egl19"] = p["pds1"] * (
            (p["pds2"] * p["pds3"]) / (1.0 + np.exp((v - p["pds4"] + p["shiftdps"]) / p["pds5"]))
            + p["pds6"]
            + (p["pds7"] * p["pds8"]) / (1.0 + np.exp((v - p["pds9"] + p["shiftdps"]) / p["pds10"]))
            + p["pds11"]
        )

        # --- CCA-1 (CaV3, T-type) -----------------------------------------
        inf["m_cca1"] = _boltzmann(v, p["va_cca1"] - p["sscca1"], p["ka_cca1"] * p["fcca"])
        tau["m_cca1"] = (
            p["p1tmcca1"]
            / (1.0 + np.exp(-(v - p["p2tmcca1"] + p["stmcca1"]) / (p["p3tmcca1"] * p["f3ca"])))
            + p["p4tmcca1"]
        ) * p["constmcca1"]
        inf["h_cca1"] = 1.0 / (
            1.0 + np.exp((v - p["vi_cca1"] + p["sshcca1"]) / (p["ki_cca1"] * p["f2cca1"]))
        )
        tau["h_cca1"] = (
            p["p1thcca1"]
            / (1.0 + np.exp((v - p["p2thcca1"] + p["sthcca1"]) / (p["p3thcca1"] * p["f4ca"])))
            + p["p4thcca1"]
        ) * p["consthcca1"]

        # --- calcium nanodomains, seen only by BK/SLO ----------------------
        # A local calcium concentration at the mouth of a channel, far higher than
        # the bulk cytosolic value and decaying over nanometres. It is what gates
        # the BK-type channels here.
        cao_nano = (
            np.abs(p["gsc"] * (v - p["eca"]) * 1e-3)
            / (8.0 * np.pi * p["r"] * p["d"] * p["F"])
            * np.exp(-p["r"] / np.sqrt(p["d"] / (p["kb"] * p["b"])))
        ) * 1e6 * 1e-3 + p["backgr"]
        cac_nano = p["backgr"]

        def rates(
            w_om: float,
            w_yx: float,
            k_yx: float,
            n_yx: float,
            w_op: float,
            w_xy: float,
            k_xy: float,
            n_xy: float,
        ) -> tuple[Any, Any, Any]:
            kcm = w_om * np.exp(-w_yx * v) / (1.0 + (cac_nano / k_yx) ** n_yx)
            kom = w_om * np.exp(-w_yx * v) / (1.0 + (cao_nano / k_yx) ** n_yx)
            kop = w_op * np.exp(-w_xy * v) / (1.0 + (k_xy / cao_nano) ** n_xy)
            return kcm, kom, kop

        # UNC-2-coupled set (BK and SLO-2) and EGL-19-coupled set (SLO-1, BK2).
        kcm_a, kom_a, kop_a = rates(
            p["wom"], p["wyx"], p["kyx"], p["nyx"], p["wop"], p["wxy"], p["kxy"], p["nxy"]
        )
        kcm_b, kom_b, kop_b = rates(
            p["wom1"],
            p["wyx1"],
            p["kyx1"],
            p["nyx1"],
            p["wop1"],
            p["wxy1"],
            p["kxy1"],
            p["nxy1"],
        )

        alpha = inf["m_unc2"] / tau["m_unc2"]
        beta = 1.0 / tau["m_unc2"] - alpha
        alpha1 = inf["m_egl19"] / tau["m_egl19"]
        beta1 = 1.0 / tau["m_egl19"] - alpha1

        def bk(
            carrier: np.ndarray, a: np.ndarray, b_: np.ndarray, kcm: Any, kom: Any, kop: Any
        ) -> tuple[Any, Any]:
            denom = (kop + kom) * (kcm + a) + b_ * kcm
            return carrier * kop * (a + b_ + kcm) / denom, (a + b_ + kcm) / denom

        inf["mbk"], tau["mbk"] = bk(s["m_unc2"], alpha, beta, kcm_a, kom_a, kop_a)
        inf["mbk2"], tau["mbk2"] = bk(s["m_egl19"], alpha1, beta1, kcm_b, kom_b, kop_b)
        # NOTE: the source computes minf_slo1 with `kop` from the UNC-2 rate set
        # while using `kop2`/`kom2`/`kcm2` -- numerically identical duplicates of
        # the UNC-2 set -- in the denominator. Reproduced as written rather than
        # "corrected", because the published behaviour is what this must match.
        inf["mslo1"], tau["mslo1"] = bk(s["m_egl19"], alpha1, beta1, kcm_a, kom_a, kop_a)
        inf["mslo2"], tau["mslo2"] = bk(s["m_unc2"], alpha, beta, kcm_b, kom_b, kop_b)

        # --- KCNL (SK), gated by bulk calcium ------------------------------
        ca = state[self.index("ca_intra1")]
        inf["m_sk"] = ca / (p["k_sk2"] + ca)
        tau["m_sk"] = np.full_like(v, p["t_sk"])

        return inf, tau

    def step(
        self, state: np.ndarray, dt_ms: float, i_ext_pa: np.ndarray | float = 0.0
    ) -> np.ndarray:
        """Advance one timestep.

        Gating variables use an exact exponential update, which is possible
        because every one of them obeys ``dx/dt = (x_inf - x)/tau`` with the
        coefficients held fixed across the step. That removes most of the
        stiffness this model would otherwise impose -- the authors ran it with a
        stiff solver at dt = 0.01 ms -- and leaves only the membrane and calcium
        equations, which are integrated explicitly.
        """
        p = self.p
        inf, tau = self._gates(state)
        currents = self.currents(state)

        nxt = state.copy()
        for name, target in inf.items():
            i = self.index(name)
            nxt[i] = target + (state[i] - target) * np.exp(-dt_ms / tau[name])

        # Intracellular calcium: influx from the three calcium currents, with a
        # first-order return to background. The source only accumulates on influx.
        i_ca = currents["unc2"] + currents["egl19"] + currents["cca1"]
        alpha_ca = 1.0 / (2.0 * p["vol"] * p["fd"])
        backgr2 = p["backgr2"]
        ca = state[self.index("ca_intra1")]
        d_ca = np.where(
            i_ca < 0.0,
            -p["fca"] * alpha_ca * i_ca - (ca - backgr2) / p["t_ca"],
            (backgr2 - ca) / p["t_ca"],
        )
        nxt[self.index("ca_intra1")] = np.maximum(ca + dt_ms * d_ca, 0.0)

        i_tot = sum(currents.values())
        nxt[0] = state[0] + dt_ms * (i_ext_pa - i_tot) / p["c"]
        return nxt

    def run(
        self,
        state: np.ndarray,
        *,
        duration_ms: float,
        i_ext_pa: np.ndarray | float = 0.0,
        dt_ms: float = DEFAULT_DT_MS,
    ) -> np.ndarray:
        """Integrate for ``duration_ms`` and return the final state.

        See :data:`DEFAULT_DT_MS` before raising ``dt_ms`` to go faster.
        """
        for _ in range(int(round(duration_ms / dt_ms))):
            state = self.step(state, dt_ms, i_ext_pa)
        return state

    def run_recorded(
        self,
        state: np.ndarray,
        *,
        duration_ms: float,
        i_ext_pa: np.ndarray | float = 0.0,
        dt_ms: float = DEFAULT_DT_MS,
        every: int = 1,
    ) -> tuple[np.ndarray, np.ndarray]:
        """As :meth:`run`, also returning the membrane potential trace.

        Separate from :meth:`run` rather than switched by a flag: a function whose
        return *type* depends on an argument cannot be checked, and the trap this
        model actually has -- see :data:`DEFAULT_DT_MS` -- is only visible in a
        trajectory, so recording one needs to be easy and unambiguous.
        """
        trace: list[np.ndarray] = []
        for k in range(int(round(duration_ms / dt_ms))):
            state = self.step(state, dt_ms, i_ext_pa)
            if k % every == 0:
                trace.append(state[0].copy())
        return state, np.array(trace)
