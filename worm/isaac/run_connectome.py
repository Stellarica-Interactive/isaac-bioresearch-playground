"""Drive the worm body with the connectome. **No scripted wave.**

    C:\\isaacsim\\python.bat worm\\isaac\\run_connectome.py --headless --seconds 20

This is the experiment the rest of the project exists to make possible. Nothing
here computes a gait. The loop is:

    body curvature
        -> proprioceptive current into B-type motor neurons   (Wen et al. 2012)
        -> 397 coupled ODEs over the measured connectome
        -> muscle cell activation, via measured neuromuscular junctions
        -> joint torque
        -> body moves, curvature changes
        -> back to the top

There is no oscillator, no pattern generator and no phase variable anywhere in
this file. If the body undulates, the undulation came out of the wiring and the
mechanics together.

Read the result carefully
-------------------------

Sustained bending is **not** by itself evidence that the connectome produces a
gait. Any feedback loop with enough gain and delay will oscillate, and this one
has a proprioceptive gain we chose. The controls that matter:

* ``--lesion`` removes named cells. If the wave survives removing the B-type
  motor neurons, it was not coming from the circuit we think.
* ``--no-proprioception`` cuts the feedback. The wave should stop, because the
  animal has no central pattern generator for forward crawling.
* ``--shuffle-sign`` randomises which synapses excite and which inhibit, keeping
  the anatomy. If behaviour survives that, it never depended on the biology.

Compare against ``run_body.py``, which drives the same body from a scripted
wave. That is the mechanical control: it shows the body *can* crawl, so a
failure here is attributable to the neural model rather than the mechanics.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace

# --- 1. Argument parsing and app launch, before any Isaac import --------------

parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
parser.add_argument("--headless", action="store_true")
parser.add_argument("--seconds", type=float, default=20.0)
parser.add_argument("--physics-hz", type=float, default=240.0)
parser.add_argument("--neural-dt-ms", type=float, default=1.0)
parser.add_argument(
    "--torque-scale",
    default="5e-4",
    help="N m per unit antagonist activation difference. Comma-separate to sweep "
    "several in one Isaac session. Purely an engineering calibration constant: "
    "see docs/model_assumptions.md.",
)
parser.add_argument("--stiffness", type=float, default=1e-4)
parser.add_argument("--damping", type=float, default=2e-5)
parser.add_argument("--armature", type=float, default=2.0e-8)
parser.add_argument("--drag-ratio", type=float, default=None)
parser.add_argument(
    "--proprioceptive-mv",
    type=float,
    default=20.0,
    help="Depolarisation of a B-type motor neuron per radian of sensed "
    "curvature, mV. Still the parameter that decides whether the loop oscillates "
    "at all, and still assumed -- but now in a unit where the assumption can be "
    "judged. The previous 400 pA/rad implied 124 mV per radian at DB1.",
)
parser.add_argument("--sensing-offset", type=float, default=None)
parser.add_argument(
    "--command",
    default="AVBL,AVBR",
    help="Forward-locomotion command interneurons to hold depolarised. This is "
    "biology, not a gait: AVB gates forward crawling and is gap-junction coupled "
    "to the B-type motor neurons. It carries no rhythm and no spatial pattern -- "
    "it is one constant number applied to two cells. Pass '' to remove it.",
)
parser.add_argument(
    "--command-mv",
    type=float,
    default=20.0,
    help="How far the command interneurons are depolarised, mV. Expressed as a "
    "voltage rather than a current because a picoamp means different things to "
    "different cells -- total conductance varies fourteenfold across this network "
    "-- and because a millivolt can be judged against the -80..+30 mV a neuron "
    "actually occupies. The current is solved for at startup.",
)
parser.add_argument(
    "--no-proprioception",
    action="store_true",
    help="control: cut the sensory feedback. The wave should stop.",
)
parser.add_argument(
    "--lesion",
    default="",
    help="control: comma-separated cells or classes to silence, e.g. DB,VB",
)
parser.add_argument(
    "--shuffle-sign",
    type=int,
    default=0,
    help="control: seed for randomising synaptic signs, keeping anatomy fixed. "
    "Behaviour that survives this never depended on the biology.",
)
parser.add_argument(
    "--noise-pa",
    type=float,
    default=0.0,
    help="Std of a fluctuating current injected into every cell, pA, resampled "
    "each neural step. A real nervous system is not noiseless, and this model's "
    "perfectly symmetric resting state is partly an artefact of determinism. "
    "The magnitude is ASSUMED, not measured.",
)
parser.add_argument("--seed", type=int, default=0, help="seed for --noise-pa")
parser.add_argument(
    "--probe",
    action="store_true",
    help="Add a draggable sphere. Move it onto the worm in the viewport (select "
    "it, press W, drag) and the body wall touch receptors fire. Touch is a "
    "transient, so unlike locomotion it is not blocked by the model's inability "
    "to oscillate -- but the escape it should trigger needs working locomotion, "
    "so expect the circuit to respond and the animal not to get away.",
)
parser.add_argument(
    "--touch-mv",
    type=float,
    default=20.0,
    help="Receptor potential a full-strength touch produces, mV. Tens of mV is "
    "the order seen in real mechanoreceptor recordings (O'Hagan, Chalfie & "
    "Goodman 2005). Solved per receptor: ALM and PLM differ twofold in "
    "conductance, so one current underdrives one of them by half.",
)
parser.add_argument(
    "--sham",
    default="",
    help="START:STOP -- measure the readout over a window with NO touch at all. "
    "The control for the touch readout: whatever this reports is what the body's "
    "own ongoing motion contributes, and a real touch has to beat it.",
)
parser.add_argument(
    "--probe-pushes",
    action="store_true",
    help="give the probe a collider so it physically shoves the worm. Off by "
    "default: a rigid probe against 0.1 g segments is an impact, not a poke, and "
    "it throws the animal across the scene.",
)
parser.add_argument(
    "--poke",
    default="",
    help="Scripted touch: START:STOP:FRACTION, seconds and position along the body "
    "(0 head, 1 tail), e.g. 3:4:0.15 for a head touch in the fourth second. Drives "
    "the same probe the mouse does, so an interactive result can be reproduced "
    "headless and put in a test.",
)
parser.add_argument(
    "--tint-change",
    action="store_true",
    help="colour the body by how the muscle drive is CHANGING rather than what it "
    "is. A touch shifts the drive by about 2%% of its resting value, which is "
    "invisible against a large static bend; this rescales to the change so it can "
    "be seen. The percentage is printed so the amplification is not mistaken for a "
    "large effect.",
)
parser.add_argument(
    "--no-tint",
    action="store_true",
    help="stop colouring the body by muscle activation. Display only -- the tint "
    "changes no number a run reports.",
)
parser.add_argument(
    "--food",
    default="",
    help="Place a food source at X,Y in body lengths from the start, e.g. 3,0. "
    "The amphid chemosensors then report the CHANGE in concentration as the "
    "animal moves -- not the concentration -- because that is the quantity a "
    "biased random walk runs on. Producing chemotaxis needs working locomotion, "
    "which 5C says we do not have; this is the sensory half.",
)
parser.add_argument(
    "--food-mv",
    type=float,
    default=20.0,
    help="Depolarisation a full-scale concentration change produces, mV.",
)
parser.add_argument(
    "--no-grid",
    action="store_true",
    help="hide the checkerboard. It is scenery with no collider, so this changes "
    "nothing a run measures.",
)
parser.add_argument(
    "--self-collision",
    action="store_true",
    help="Stop the body passing through itself when it folds. Off by default "
    "because it adds contact forces beside the drag model that stands in for "
    "the substrate, and every W2 number was measured without it.",
)
parser.add_argument("--report-every", type=float, default=2.0)
args, _ = parser.parse_known_args()

from isaacsim import SimulationApp  # noqa: E402

simulation_app = SimulationApp({"headless": args.headless})

# --- 2. Everything else -------------------------------------------------------

import numpy as np  # noqa: E402
import omni.timeline  # noqa: E402
from isaacsim.core.experimental.prims import Articulation, RigidPrim  # noqa: E402
from isaacsim.core.simulation_manager import SimulationManager  # noqa: E402

from common.data.schemas import CellCategory, Connectome, Sign  # noqa: E402
from common.neural.runtime import NeuralRuntime  # noqa: E402
from common.neural.stimulus import (  # noqa: E402
    scale_for_depolarisation,
    solve_for_depolarisation,
)
from common.neural.synapses import UnknownSignPolicy  # noqa: E402
from worm.body.chemotaxis import (  # noqa: E402
    CHEMOSENSORS,
    DEFAULT_DECAY_BODY_LENGTHS,
    Chemosensation,
    FoodSource,
)
from worm.body.drag import DragParameters, GroundDrag  # noqa: E402
from worm.body.geometry import BodyPlan  # noqa: E402
from worm.body.muscles import MuscleModel, MuscleParameters  # noqa: E402
from worm.body.neural_bridge import MuscleDrive, Proprioception  # noqa: E402
from worm.body.touch import TOUCH_RECEPTORS, TouchField  # noqa: E402
from worm.importers.naming import body_wall_muscle_ids  # noqa: E402
from worm.isaac.stage import (  # noqa: E402
    DEFAULT_PROBE_RADIUS_SCALE,
    ActivityTint,
    add_camera,
    add_ground_grid,
    add_probe,
    build_scene,
    dof_order,
)  # noqa: E402
from worm.loader import load  # noqa: E402
from worm.neural.config import RUNTIME_OVERLAYS, build_runtime  # noqa: E402


def main() -> int:
    torque_scales = [float(s) for s in str(args.torque_scale).split(",")]
    plan = BodyPlan()

    connectome, _ = load("cook_2019_herm", annotations=RUNTIME_OVERLAYS)
    if args.shuffle_sign:
        connectome = _shuffle_signs(connectome, args.shuffle_sign)
        print(f"CONTROL: synaptic signs shuffled with seed {args.shuffle_sign}")

    muscles_wanted = body_wall_muscle_ids()
    cells = tuple(
        c.id
        for c in connectome.cells
        if c.category is CellCategory.NEURON or c.id in muscles_wanted
    )
    runtime, report = build_runtime(
        "cook_2019_herm",
        unknown_sign=UnknownSignPolicy.EXCLUDE,
        cells=cells,
        connectome=connectome,
        dt_ms=args.neural_dt_ms,
    )
    print(report.summary())

    lesioned = _resolve_lesion(connectome, args.lesion)
    lesion_idx = np.array(
        [runtime.network.index(c) for c in sorted(lesioned) if c in runtime.network.cell_ids],
        dtype=np.int64,
    )
    if lesioned:
        print(f"\nCONTROL: silencing {len(lesion_idx)} cells: {sorted(lesioned)[:8]} ...")

    bridge = MuscleDrive.build(plan, runtime.network.cell_ids)
    proprio = Proprioception.build(
        connectome,
        plan,
        **({"offset": args.sensing_offset} if args.sensing_offset else {}),
    )
    b_type = [c for c in proprio.targets if c in runtime.network.cell_ids]
    proprio = replace(
        proprio,
        gain_pa_per_rad=scale_for_depolarisation(runtime, b_type, args.proprioceptive_mv),
    )
    # Every input is a target depolarisation, converted to a current against this
    # network's own conductances. Done once here rather than per step: each solve
    # settles the network several times. See common/neural/stimulus.py for why the
    # targets are voltages and not picoamps.
    print("\n  solving inputs for their target depolarisations ...")
    runtime.run(2000.0)

    command_cells = [
        c.strip()
        for c in args.command.split(",")
        if c.strip() and c.strip() in runtime.network.cell_ids
    ]
    command = (
        solve_for_depolarisation(runtime, command_cells, args.command_mv) if command_cells else {}
    )

    food = None
    chemo = None
    if args.food:
        fx, fy = (float(v) for v in args.food.split(","))
        food = FoodSource(
            x_m=fx * plan.total_length_m,
            y_m=fy * plan.total_length_m,
            decay_m=DEFAULT_DECAY_BODY_LENGTHS * plan.total_length_m,
        )
        sensor_cells = [s.cell for s in CHEMOSENSORS if s.cell in runtime.network.cell_ids]
        chemo = Chemosensation.build(
            cells=runtime.network.cell_ids,
            per_cell_pa=solve_for_depolarisation(runtime, sensor_cells, args.food_mv),
        )

    receptor_cells = [r.cell for r in TOUCH_RECEPTORS if r.cell in runtime.network.cell_ids]
    per_cell = solve_for_depolarisation(runtime, receptor_cells, args.touch_mv)
    touch = TouchField.build(plan, cells=runtime.network.cell_ids, per_cell_pa=per_cell)
    print(
        f"\n  muscle slots driven by the connectome: {bridge.covered}/{4 * plan.n_segments}"
        f"\n  proprioceptive targets: {len(proprio.targets)} B-type neurons"
        f"{' (DISABLED)' if args.no_proprioception else ''}"
        f"\n  proprioceptive gain: {args.proprioceptive_mv:g} mV/rad"
        f" = {proprio.gain_pa_per_rad:.1f} pA/rad"
        f"\n  noise: {args.noise_pa:g} pA std (seed {args.seed})"
        f"\n  command drive: {', '.join(command) or 'none'} at "
        f"{args.command_mv:g} mV = "
        f"{np.mean(list(command.values())) if command else 0.0:.1f} pA "
        f"(constant, no rhythm)"
        + (
            f"\n  touch receptors: {', '.join(touch.cells)} at "
            f"{args.touch_mv:g} mV "
            f"({min(per_cell.values()):.0f}-{max(per_cell.values()):.0f} pA)"
            f" -- drag the red sphere onto the body"
            if args.probe
            else ""
        )
        + (
            f"\n  food at ({args.food}) BL: {', '.join(chemo.cells)} at "
            f"{args.food_mv:g} mV per unit change"
            if chemo is not None
            else ""
        )
    )

    dt = 1.0 / args.physics_hz
    neural_substeps = max(1, int(round(dt * 1000.0 / args.neural_dt_ms)))
    print(f"  {neural_substeps} neural steps per physics step")

    # The neural state at rest, so every condition in a sweep starts from the same
    # nervous system rather than inheriting the previous one's.
    rest_state = runtime.state.copy()

    for scale in torque_scales:
        if len(torque_scales) > 1:
            print(f"\n--- peak_torque_scale = {scale:g} " + "-" * 40)
        runtime.state = rest_state.copy()
        runtime.t_ms = 0.0
        _run_condition(
            plan=plan,
            runtime=runtime,
            touch=touch,
            food=food,
            chemo=chemo,
            rng=np.random.default_rng(args.seed),
            noise_pa=args.noise_pa,
            bridge=bridge,
            proprio=proprio,
            command=command,
            lesion_idx=lesion_idx,
            torque_scale=scale,
            dt=dt,
            neural_substeps=neural_substeps,
        )

    print(
        "\n  amp is how much the body is bending *over time*, so a body frozen in a\n"
        "  bent shape reads zero. travel is head-to-tail propagation: +1 a clean\n"
        "  wave, 0 a standing oscillation, and it subtracts the reverse direction\n"
        "  so a static bend cannot fake it. pinned is the fraction of joints against\n"
        "  the 60 deg limit -- above zero the mechanics are setting the body shape.\n"
        "  Compare against --no-proprioception and --lesion DB,VB before believing\n"
        "  any of it."
    )
    simulation_app.close()
    return 0


def _run_condition(  # noqa: PLR0913 - one experimental condition, all of it explicit
    *,
    plan: BodyPlan,
    runtime: NeuralRuntime,
    touch: TouchField,
    food: FoodSource | None,
    chemo: Chemosensation | None,
    rng: np.random.Generator,
    noise_pa: float,
    bridge: MuscleDrive,
    proprio: Proprioception,
    command: dict[str, float],
    lesion_idx: np.ndarray,
    torque_scale: float,
    dt: float,
    neural_substeps: int,
) -> None:
    """One closed-loop run on a freshly built body."""
    muscle_model = MuscleModel(
        plan,
        MuscleParameters(
            peak_torque_scale=torque_scale,
            joint_stiffness=args.stiffness,
            joint_damping=args.damping,
        ),
    )
    drag = GroundDrag(
        plan, DragParameters(**({"ratio": args.drag_ratio} if args.drag_ratio else {}))
    )

    root_path, joint_paths = build_scene(plan, self_collision=args.self_collision)
    camera = add_camera(plan)
    if not args.no_grid:
        add_ground_grid(plan)
    poke = [float(x) for x in args.poke.split(":")] if args.poke else None
    sham = [float(x) for x in args.sham.split(":")] if args.sham else None
    probe_path = add_probe(plan, collider=args.probe_pushes) if (args.probe or poke) else None
    probe = RigidPrim(probe_path) if probe_path else None
    tint = ActivityTint.build(plan.n_segments, probe_path=probe_path) if not args.no_tint else None
    SimulationManager.set_physics_dt(dt)
    masses = plan.masses_kg()
    articulation = Articulation(root_path)
    links = RigidPrim([f"{root_path}/segment_{i:02d}" for i in range(plan.n_segments)])

    omni.timeline.get_timeline_interface().play()
    simulation_app.update()

    dofs = dof_order(articulation, joint_paths)
    articulation.set_dof_gains(
        stiffnesses=muscle_model.params.joint_stiffness,
        dampings=muscle_model.params.joint_damping,
    )
    articulation.set_dof_position_targets(np.zeros((1, plan.n_joints)), dof_indices=dofs)
    if args.armature:
        articulation.set_dof_armatures(args.armature)

    start = _centroid(links)
    touching = False
    # Rolling history, so a touch is always reported against what the network was
    # doing anyway. Without this the readout cannot tell a response from noise --
    # and with --noise-pa on it reported pure noise as a reversal.
    command_history: list[tuple[float, dict[str, float]]] = []
    contact_started_at = 0.0
    drive_reference: np.ndarray | None = None
    drive_at_contact: np.ndarray | None = None
    peak: dict[str, float] = {}
    baseline: dict[str, float] = {}
    history: list[np.ndarray] = []
    reported = 0.0
    t = 0.0
    angles = np.zeros(plan.n_joints)

    for _ in range(int(args.seconds * args.physics_hz)):
        angles = np.asarray(
            articulation.get_dof_positions(dof_indices=dofs), dtype=np.float64
        ).reshape(-1)

        # --- body -> nervous system -------------------------------------
        runtime.clear_inputs()
        runtime.inject_many(command)
        if chemo is not None and food is not None:
            # Sampled at the nose, because that is where the amphid openings
            # are, and because a worm swinging its head samples a gradient the
            # body centre never sees.
            nose = np.asarray(links.get_world_poses()[0])[0, :2]
            chemo.step(food.concentration(nose), dt_ms=dt * 1000.0)
            runtime.inject_many(chemo.currents())
        if not args.no_proprioception:
            runtime.inject_many(proprio.currents(angles))
        if sham is not None:
            # Same statistic, no stimulus. Anything it reports is confound.
            if sham[0] <= t <= sham[1]:
                if not touching:
                    baseline = _command_state(runtime)
                    drive_at_contact = (
                        muscle_model.dorsal_activation() - muscle_model.ventral_activation()
                    )
                    print(f"    t={t:5.1f}s  SHAM window opens (no touch)")
                touching = True
                peak = {
                    k: max(peak.get(k, 0.0), v - baseline[k], key=abs)
                    for k, v in _command_state(runtime).items()
                }
            elif touching:
                now = muscle_model.dorsal_activation() - muscle_model.ventral_activation()
                assert drive_at_contact is not None
                shift = float(np.abs(now - drive_at_contact).max())
                rest = max(float(np.abs(drive_at_contact).max()), 1e-12)
                print(
                    f"    t={t:5.1f}s  SHAM closes  -> "
                    + ", ".join(f"{k} {v:+.2f}mV" for k, v in sorted(peak.items()))
                    + f"; muscle drive moved {100 * shift / rest:.1f}%"
                )
                touching, peak = False, {}
        elif probe is not None:
            if poke is not None:
                _drive_poke(probe, links, plan, t, poke)
            contact, ventral = _probe_contact(probe, links, plan)
            if contact.any():
                currents = touch.currents_from_segments(contact, ventral=ventral)
                runtime.inject_many(currents)
                if not touching:
                    hit = np.flatnonzero(contact)
                    baseline = _command_state(runtime)
                    contact_started_at = t
                    drive_at_contact = (
                        muscle_model.dorsal_activation() - muscle_model.ventral_activation()
                    )
                    print(
                        f"    t={t:5.1f}s  TOUCH segments {hit.min()}-{hit.max()} "
                        f"({touch.segment_fraction(int(hit.mean())):.2f} along body) "
                        f"-> {', '.join(f'{c} {v:.1f}pA' for c, v in sorted(currents.items()))}"
                    )
                touching = True
                peak = {
                    k: max(peak.get(k, 0.0), v - baseline[k], key=abs)
                    for k, v in _command_state(runtime).items()
                }
            elif touching:
                # The result that matters. Anterior touch should raise the backward
                # command (AVA) and lower the forward one (AVB); posterior touch the
                # reverse. Reported as a change from the moment of contact, because
                # the absolute values are set by the standing command drive.
                muscle_note = ""
                if drive_at_contact is not None:
                    now = muscle_model.dorsal_activation() - muscle_model.ventral_activation()
                    shift = float(np.abs(now - drive_at_contact).max())
                    rest = max(float(np.abs(drive_at_contact).max()), 1e-12)
                    muscle_note = f"; muscle drive moved {100 * shift / rest:.1f}%"
                background = _background_swing(
                    command_history, contact_started_at, t - contact_started_at
                )
                verdict = _verdict(peak, background)
                print(
                    f"    t={t:5.1f}s  released -> "
                    + ", ".join(f"{k} {v:+.2f}mV" for k, v in sorted(peak.items()))
                    + muscle_note
                )
                print(
                    "                 background over an equal window with no touch: "
                    + ", ".join(f"{k} {v:.2f}mV" for k, v in sorted(background.items()))
                    + f"   -> {verdict}"
                )
                touching, peak = False, {}

        # --- nervous system ---------------------------------------------
        for _ in range(neural_substeps):
            if noise_pa:
                # Injected as a current rather than added inside the integrator, so
                # the runtime stays bit-for-bit deterministic and this is visibly an
                # input to the model rather than a change to it. Resampled per
                # neural step; the magnitude is ours, not a measurement.
                runtime.i_ext_pa += rng.normal(0.0, noise_pa, size=runtime.network.n)
            runtime.step()
            if lesion_idx.size:
                # Hold ablated cells at rest so they transmit nothing. Applied
                # after the step rather than by deleting them, so the anatomy is
                # unchanged and only the cell's output is removed -- which is
                # what an ablation does.
                runtime.state[1, lesion_idx] = 0.0

        # --- nervous system -> body -------------------------------------
        muscle_model.step(bridge.drive(runtime.state[1]), dt_ms=dt * 1000.0)
        command_history.append((t, _command_state(runtime)))
        if len(command_history) > 2400:  # ten seconds at 240 Hz
            command_history.pop(0)
        drive_now = muscle_model.dorsal_activation() - muscle_model.ventral_activation()
        if drive_reference is None and t > 1.0:
            drive_reference = drive_now.copy()
        if tint is not None and len(history) % 6 == 0:
            tint.update(
                drive_now,
                touching=touching,
                reference=drive_reference if args.tint_change else None,
            )
        articulation.set_dof_efforts(muscle_model.joint_torques().reshape(1, -1), dof_indices=dofs)
        _apply_drag(links, drag, dt, masses)

        simulation_app.update()
        t += dt
        history.append(angles.copy())
        # Keep the animal in frame; it travels several body lengths.
        if len(history) % 4 == 0:
            camera.follow(_centroid(links))

        if t - reported >= args.report_every:
            reported = t
            _report(t, start, links, plan, angles, history)

    _report(args.seconds, start, links, plan, angles, history, final=True)


def _background_swing(
    history: list[tuple[float, dict[str, float]]], contact_at: float, duration: float
) -> dict[str, float]:
    """How much each command group moved on its own, just before the touch.

    Every touch needs its own control, because "how much did AVB change while I was
    pressing" is not a measurement of the touch -- it is a measurement of the touch
    plus whatever the network was already doing. With ``--noise-pa`` on, the second
    term dominates completely: a sham window with no contact at all reports the same
    -21 mV swing as a real touch does.

    So the same statistic is computed over an equal-length window ending at the
    moment of contact, and printed beside the response. A response that does not
    exceed its own background is not a response.
    """
    window = [s for time, s in history if contact_at - duration <= time < contact_at]
    if len(window) < 2:
        return {}
    keys = window[0].keys()
    return {k: float(max(abs(s[k] - window[0][k]) for s in window)) for k in keys}


def _verdict(peak: dict[str, float], background: dict[str, float]) -> str:
    """Whether the response is distinguishable from the network's own fluctuation."""
    if not background:
        return "no background window yet"
    clears = {k: abs(v) for k, v in peak.items() if abs(v) > 2.0 * background.get(k, 0.0)}
    if not clears:
        return "INDISTINGUISHABLE FROM BACKGROUND -- this is not a touch response"
    reversal = peak.get("AVA", 0.0) > 0 > peak.get("AVB", 0.0)
    tag = "; AVA up + AVB down = reversal" if reversal else ""
    return f"clears background: {', '.join(sorted(clears))}{tag}"


def _command_state(runtime: NeuralRuntime) -> dict[str, float]:
    """Mean membrane potential of the forward and backward command interneurons.

    AVA/AVD/AVE drive backward locomotion and AVB/PVC forward; the touch circuit
    reaches the body only through them, so they are where a response has to appear
    if the pathway works at all.
    """
    groups = {
        "AVA": ("AVAL", "AVAR"),
        "AVB": ("AVBL", "AVBR"),
        "AVD": ("AVDL", "AVDR"),
        "PVC": ("PVCL", "PVCR"),
    }
    out: dict[str, float] = {}
    for name, cells in groups.items():
        present = [c for c in cells if c in runtime.network.cell_ids]
        if present:
            out[name] = float(np.mean([runtime.voltage(c) for c in present]))
    return out


def _drive_poke(
    probe: RigidPrim, links: RigidPrim, plan: BodyPlan, t: float, poke: list[float]
) -> None:
    """Move the probe onto the body between two times, then take it away.

    The scripted equivalent of dragging with the mouse. It exists so a touch
    experiment is reproducible and can be asserted in a test: an interactive demo
    that cannot be re-run identically is not evidence of anything.
    """
    start, stop, fraction = poke[0], poke[1], poke[2]
    positions = np.asarray(links.get_world_poses()[0])[:, :3]
    segment = int(np.clip(round(fraction * (plan.n_segments - 1)), 0, plan.n_segments - 1))

    if start <= t <= stop:
        # Placed against the *surface*, not the centre.
        #
        # Teleporting a collider to positions[segment] puts the sphere entirely
        # inside the animal. PhysX resolves that interpenetration explosively, the
        # body flails, and the failure then surfaces somewhere else entirely: joint
        # angles go wild, proprioception carries them into the network, and the
        # neural runtime reports a divergence that looks like a neuroscience
        # problem. It is not; it is a spawn position.
        axis = np.gradient(positions[:, :2], axis=0)[segment]
        norm = np.linalg.norm(axis)
        normal = np.array([-axis[1], axis[0]]) / norm if norm > 1e-12 else np.array([0.0, 1.0])
        # Just far enough in to register as a firm touch, not far enough to shove.
        reach = plan.segments()[segment].radius_m + DEFAULT_PROBE_RADIUS_SCALE * plan.max_radius_m
        target = positions[segment].copy()
        # Resting against the surface. With no collider there is nothing to stop
        # the probe going deeper, but there is no reason to: the proximity margin
        # already reads a body-contacting probe as a full-strength touch.
        target[:2] += normal * reach
    else:
        # Parked well clear of the animal, where it contacts nothing.
        target = positions.mean(axis=0) + np.array([0.0, 0.5 * plan.total_length_m, 0.0])
    if not np.all(np.isfinite(target)):
        return  # the body has already diverged; do not feed PhysX a NaN pose
    # Orientation is passed explicitly: PhysX rejects a pose whose quaternion it
    # cannot read, and an omitted one is not necessarily identity.
    probe.set_world_poses(
        positions=np.asarray(target, dtype=np.float32).reshape(1, 3),
        orientations=np.array([[1.0, 0.0, 0.0, 0.0]], dtype=np.float32),
    )


def _probe_contact(
    probe: RigidPrim, links: RigidPrim, plan: BodyPlan
) -> tuple[np.ndarray, np.ndarray]:
    """Which segments the probe is pressing on, and from which side.

    Contact is measured geometrically here rather than read from PhysX contact
    reports. The probe is kinematic and the worm's shape is known, so sphere-to-
    capsule overlap is exact, cheap, and -- unlike a contact report -- gives a
    graded depth we can use as stimulus strength. It also means the sensory signal
    is computed the same way whether or not the probe is physically pushing.

    The model is planar, so 'ventral' is a side of the body axis rather than a
    direction in space: the sign of the cross product between the local body axis
    and the direction to the probe. That distinction only matters for AVM and PVM,
    the two genuinely ventral receptors.
    """
    centre = np.asarray(probe.get_world_poses()[0]).reshape(-1)[:3]
    positions = np.asarray(links.get_world_poses()[0])[:, :3]
    seg_radii = np.array([s.radius_m for s in plan.segments()])
    probe_radius = DEFAULT_PROBE_RADIUS_SCALE * plan.max_radius_m

    delta = positions - centre
    distance = np.linalg.norm(delta, axis=1)
    # A margin standing in for the compliance a rigid body does not have.
    #
    # Without it this sensor is blind to the mouse. The probe is a rigid collider,
    # so PhysX holds it *outside* the body: the centre-to-centre distance settles at
    # exactly seg_radius + probe_radius and never goes below it, which makes any
    # test based on interpenetration read zero no matter how hard you press. The
    # scripted poke only appeared to work because it teleports the probe inside for
    # a single frame before the solver ejects it.
    #
    # A real worm is soft and a real touch indents it; our capsules cannot deform,
    # so the indentation that would have happened is represented by this margin.
    # Resting against the body counts as a firm touch, and the response falls to
    # zero as the probe is lifted a margin clear.
    margin = 0.5 * seg_radii
    overlap = (seg_radii + probe_radius + margin) - distance
    contact = np.clip(overlap / margin, 0.0, 1.0)
    contact[overlap <= 0.0] = 0.0

    # Body axis at each segment, from its neighbours.
    axis = np.gradient(positions[:, :2], axis=0)
    to_probe = -delta[:, :2]
    cross = axis[:, 0] * to_probe[:, 1] - axis[:, 1] * to_probe[:, 0]
    return contact, cross < 0.0


# --- controls -----------------------------------------------------------------


def _resolve_lesion(connectome: Connectome, spec: str) -> set[str]:
    """Expand a comma-separated list of cell ids or class names into cell ids."""
    wanted = {s.strip() for s in spec.split(",") if s.strip()}
    if not wanted:
        return set()
    return {
        c.id
        for c in connectome.cells
        if c.id in wanted or (c.class_name and c.class_name in wanted)
    }


def _shuffle_signs(connectome: Connectome, seed: int) -> Connectome:
    """Randomly reassign excitatory/inhibitory, keeping every synapse in place."""
    from dataclasses import replace

    rng = np.random.default_rng(seed)
    signed = [Sign.EXCITATORY, Sign.INHIBITORY]
    return replace(
        connectome,
        connections=tuple(
            replace(e, sign=signed[int(rng.integers(2))])
            if e.sign in (Sign.EXCITATORY, Sign.INHIBITORY)
            else e
            for e in connectome.connections
        ),
    )


# --- reporting ----------------------------------------------------------------


def _centroid(links: RigidPrim) -> np.ndarray:
    return np.asarray(links.get_world_poses()[0])[:, :2].mean(axis=0)


def _apply_drag(links: RigidPrim, drag: GroundDrag, dt_s: float, masses: np.ndarray) -> None:
    positions, _ = links.get_world_poses()
    xy = np.asarray(positions)[:, :2]
    vxy = np.asarray(links.get_velocities()[0])[:, :2]
    forces = np.zeros((len(xy), 3), dtype=np.float32)
    forces[:, :2] = drag.forces(xy, vxy, dt_s=dt_s, masses_kg=masses)
    if not np.all(np.isfinite(forces)):
        raise FloatingPointError("body diverged")
    links.apply_forces(forces)


def _wave_metrics(history: list[np.ndarray]) -> tuple[float, float]:
    """How much the body is bending, and whether the bend is *travelling*.

    Amplitude is the temporal standard deviation of joint angle, so a body frozen
    in a bent shape reads zero however bent it is.

    Travel is the harder measurement and the one worth getting right. A posterior
    joint doing what its anterior neighbour did a moment ago is a wave moving head
    to tail. A body flexing in place does the same thing at positive and negative
    lag, so correlating at one lag alone cannot tell the two apart -- and reports a
    large number for a standing oscillation, which is exactly the mistake this
    function made in its first version. The answer is the *difference* between the
    two directions, which is zero for anything symmetric in time:

        travel = corr(anterior(t), posterior(t+lag)) - corr(anterior(t), posterior(t-lag))

    Each joint's own temporal mean is removed first, so a static bend contributes
    nothing. +1 is a clean head-to-tail wave, -1 tail-to-head, 0 no propagation.
    """
    if len(history) < 120:
        return 0.0, 0.0
    recent = np.array(history[-240:])
    fluct = recent - recent.mean(axis=0, keepdims=True)
    amplitude = float(np.degrees(fluct.std()))

    lag = min(12, len(fluct) // 4)
    anterior = fluct[lag:-lag, :-1]

    def corr(shift: int) -> float:
        posterior = fluct[lag + shift : len(fluct) - lag + shift, 1:]
        a, b = anterior.ravel(), posterior.ravel()
        if a.std() < 1e-12 or b.std() < 1e-12:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    return amplitude, corr(lag) - corr(-lag)


def _report(  # noqa: PLR0913
    t: float,
    start: np.ndarray,
    links: RigidPrim,
    plan: BodyPlan,
    angles: np.ndarray,
    history: list[np.ndarray],
    final: bool = False,
) -> None:
    xyz = np.asarray(links.get_world_poses()[0])
    if not np.all(np.isfinite(xyz)):
        print(f"  t={t:5.1f}s   DIVERGED")
        return
    delta = np.asarray(links.get_world_poses()[0])[:, :2].mean(axis=0) - start
    extent = float(np.linalg.norm(xyz[-1] - xyz[0])) / plan.total_length_m
    amplitude, travel = _wave_metrics(history)
    label = "FINAL" if final else f"t={t:5.1f}s"
    # Heading tells crawling apart from being stirred. A body held in a fixed bend
    # can still travel a long way by rotating against anisotropic drag, and the
    # distance moved alone cannot distinguish that from locomotion.
    axis = xyz[-1, :2] - xyz[0, :2]
    heading = float(np.degrees(np.arctan2(axis[1], axis[0])))
    # Fraction of joints pinned against the 60 deg limit. Anything above zero
    # means the mechanics, not the nervous system, are setting the body shape.
    saturated = float(np.mean(np.abs(np.degrees(angles)) > 59.0))
    print(
        f"  {label}  moved {np.linalg.norm(delta) * 1e3:7.2f} mm "
        f"({np.linalg.norm(delta) / plan.total_length_m:5.3f} BL)  "
        f"bend {np.degrees(np.abs(angles).max()):5.1f}deg  "
        f"amp {amplitude:5.2f}deg  travel {travel:+5.2f}  "
        f"extent {extent:4.2f}  pinned {saturated:4.0%}  head {heading:+7.1f}deg"
    )
    if final:
        print(
            "\n  amp is how much the body is bending; travel is whether the bend "
            "moves along it\n  (0 = standing oscillation, positive = wave "
            "propagating head to tail).\n  Compare against --no-proprioception "
            "and --lesion DB,VB before believing any of it."
        )


if __name__ == "__main__":
    sys.exit(main())
