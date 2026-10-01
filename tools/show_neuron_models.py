"""Print every conductance-based neuron model this project can run.

    python tools/show_neuron_models.py

Each cell is settled from its own initial condition and reported with the sum of
its currents at that state. That sum is the useful column: nothing in the code
forces it to zero, so a value near zero is evidence the channel kinetics and the
per-cell conductances were assembled consistently, and a model with a channel at
the wrong density will show a residual.

It is NOT evidence of agreement with the published figures -- see
docs/model_assumptions.md 5AB.6, which says what has and has not been checked.
"""

from __future__ import annotations

import numpy as np

from common.neural.conductance import ConductanceModel, channels_for, model_provenance

MODELS = (
    "rmd_nicoletti2019",
    "awc_nicoletti2019",
    "aval_nicoletti2024",
    "avar_nicoletti2024",
    "aiy_nicoletti2024",
    "rim_nicoletti2024",
    "va5_nicoletti2024",
    "vd5_nicoletti2024",
    "vb6_nicoletti2024",
)


def main() -> int:
    header = f"{'model':24s} {'ch':>3} {'states':>6} {'resting':>11} {'sum I':>11}"
    print(header)
    print("-" * len(header))
    for model_id in MODELS:
        model = ConductanceModel.load(model_id)
        rest = model.resting_state()
        total = float(sum(np.atleast_1d(v)[0] for v in model.currents(rest).values()))
        print(
            f"{model_id:24s} {len(channels_for(model_id)):3d} "
            f"{len(model.state_names):6d} {float(rest[0][0]):+10.2f} mV "
            f"{total:+11.5f} pA"
        )
    print("\nchannels per cell:")
    for model_id in MODELS:
        print(f"  {model_id:24s} {', '.join(channels_for(model_id))}")
    print("\nprovenance:")
    seen = set()
    for model_id in MODELS:
        source = model_provenance(model_id)
        key = source.get("source_id", "?")
        if key in seen:
            continue
        seen.add(key)
        print(f"  {key}: {source.get('citation', '')[:96]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
