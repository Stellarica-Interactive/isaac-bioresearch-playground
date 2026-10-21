"""Presentation: a skin over the segments, an agar ground, and the connectome drawn live.

Display only. Everything here writes to the stage's appearance and reads the
simulation; nothing feeds back, so a run with any of it measures identically to
one without. The geometry is computed in :mod:`worm.body.visuals`, which is
tested without Isaac Sim; this module only turns it into USD.

* :class:`WormSkin` -- one smooth tube over the 24 capsules, a flesh tone shifted
  towards the orange or blue :class:`~worm.isaac.stage.ActivityTint` paints for
  the muscle drive. The capsules stay, invisible, as the physics bodies they are.
* :func:`add_agar_ground` -- a textured plane in place of the checkerboard. The
  texture is generated at run time, so nothing binary enters the repository.
* :func:`presentation_lighting` -- a softer key light and a fill, which the skin
  and agar colours are chosen under.
* :class:`ConnectomeOverlay` -- the neurons as spheres and the strongest synapses
  as curves that run along the body, coloured by synaptic activation: ``inside``
  the body's outline (drawn over the top of the skin, which from the camera above
  reads as inside, with no transparency to depend on), ``above`` it and beside
  it, or on ``screen`` as a panel.
* :class:`ViewportCapture` -- screenshots at intervals, for presentations and for
  checking a change without watching the window.

Colours given to ``displayColor`` go through :func:`~worm.body.visuals.linear`
first: RTX reads them as linear light, and a colour picked by eye renders washed
out otherwise.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from worm.body.geometry import BodyPlan
from worm.body.visuals import (
    Beside,
    ConnectomeDrawing,
    ConnectomePanel,
    activity_colours,
    agar_texture,
    encode_png,
    linear,
    ring_colours,
    skin_geometry,
    skin_rings,
    skin_tint,
    skin_topology,
    sphere,
)

CONNECTOME_MODES = ("off", "screen", "above", "inside")

#: How far beside the worm the ``above`` drawing floats, in body lengths: clear
#: of the body's sideways sweep and of the drawing's own width.
ABOVE_OFFSET_BODY_LENGTHS = 0.3


def _straight_centres(plan: BodyPlan, z: float) -> np.ndarray:
    return np.array([[s.centre_x_m, 0.0, z] for s in plan.segments()])


@dataclass
class WormSkin:
    """A smooth body over the segments, re-swept along the spine every update."""

    mesh: Any
    plan: BodyPlan
    per_segment: int = 4
    n_around: int = 16

    @classmethod
    def build(cls, plan: BodyPlan, *, root_path: str, path: str = "/World/WormSkin") -> WormSkin:
        """Coloured by ``displayColor`` with no material bound: RTX draws a
        prim's display colour directly, and ignores it once a UsdPreviewSurface
        reading it through a primvar reader is bound (found by screenshot)."""
        import isaacsim.core.experimental.utils.stage as stage_utils
        from pxr import UsdGeom, Vt

        stage = stage_utils.get_current_stage()
        skin = cls(mesh=UsdGeom.Mesh.Define(stage, path), plan=plan)
        n_rings = skin_rings(plan.n_segments, skin.per_segment)
        counts, indices = skin_topology(n_rings, skin.n_around)
        skin.mesh.CreateFaceVertexCountsAttr().Set(Vt.IntArray(counts.tolist()))
        skin.mesh.CreateFaceVertexIndicesAttr().Set(Vt.IntArray(indices.tolist()))
        skin.mesh.CreateSubdivisionSchemeAttr().Set(UsdGeom.Tokens.none)
        skin.mesh.CreateDoubleSidedAttr().Set(False)
        skin.mesh.SetNormalsInterpolation(UsdGeom.Tokens.vertex)
        # Authored as per-vertex explicitly. The schema supplies a displayColor
        # whether or not one is authored, constant by default, and with that RTX
        # paints the whole body the first colour -- the tint never showed.
        skin.mesh.CreateDisplayColorPrimvar(UsdGeom.Tokens.vertex)
        skin.update(_straight_centres(plan, plan.max_radius_m))
        # The capsules remain the physics bodies; only their drawing is replaced.
        for i in range(plan.n_segments):
            prim = stage.GetPrimAtPath(f"{root_path}/segment_{i:02d}")
            if prim.IsValid():
                UsdGeom.Imageable(prim).MakeInvisible()
        return skin

    def update(self, centres: np.ndarray, painted: np.ndarray | None = None) -> None:
        """``painted`` is what :class:`ActivityTint` gave each segment, or None."""
        from pxr import Vt

        points, normals = skin_geometry(
            centres, self.plan, per_segment=self.per_segment, n_around=self.n_around
        )
        self.mesh.GetPointsAttr().Set(Vt.Vec3fArray.FromNumpy(points.astype(np.float32)))
        self.mesh.GetNormalsAttr().Set(Vt.Vec3fArray.FromNumpy(normals.astype(np.float32)))
        lo, hi = points.min(axis=0), points.max(axis=0)
        self.mesh.GetExtentAttr().Set(
            Vt.Vec3fArray.FromNumpy(np.vstack([lo, hi]).astype(np.float32))
        )
        n_rings = skin_rings(self.plan.n_segments, self.per_segment)
        colours = ring_colours(skin_tint(painted, self.plan.n_segments), n_rings, self.n_around)
        self.mesh.GetDisplayColorPrimvar().Set(
            Vt.Vec3fArray.FromNumpy(linear(colours).astype(np.float32))
        )


def add_agar_ground(
    plan: BodyPlan,
    *,
    path: str = "/World/Agar",
    extent_body_lengths: float = 24.0,
    tile_body_lengths: float = 1.5,
    texture_size: int = 1024,
    seed: int = 0,
) -> str:
    """A textured plane under the worm, like an agar plate. Scenery, as the grid is.

    No collider and no rigid body, for the reason :func:`add_ground_grid` gives:
    the substrate is the drag model, and a contact surface would be a second one.
    """
    import isaacsim.core.experimental.utils.stage as stage_utils
    from pxr import Gf, Sdf, UsdGeom, UsdShade, Vt

    stage = stage_utils.get_current_stage()
    folder = Path(tempfile.gettempdir()) / "isaac_worm_textures"
    folder.mkdir(parents=True, exist_ok=True)
    texture = folder / f"agar_{texture_size}_{seed}.png"
    # Rewritten every run: it takes well under a second, and a cached file would
    # outlive any change to how the texture is made.
    texture.write_bytes(encode_png(agar_texture(texture_size, seed)))

    half = 0.5 * extent_body_lengths * plan.total_length_m
    z = -1.05 * plan.max_radius_m
    tiles = extent_body_lengths / tile_body_lengths
    mesh = UsdGeom.Mesh.Define(stage, path)
    mesh.CreatePointsAttr().Set(
        Vt.Vec3fArray(
            [
                Gf.Vec3f(-half, -half, z),
                Gf.Vec3f(half, -half, z),
                Gf.Vec3f(half, half, z),
                Gf.Vec3f(-half, half, z),
            ]
        )
    )
    mesh.CreateFaceVertexCountsAttr().Set(Vt.IntArray([4]))
    mesh.CreateFaceVertexIndicesAttr().Set(Vt.IntArray([0, 1, 2, 3]))
    mesh.CreateExtentAttr().Set(Vt.Vec3fArray([Gf.Vec3f(-half, -half, z), Gf.Vec3f(half, half, z)]))
    st = UsdGeom.PrimvarsAPI(mesh).CreatePrimvar(
        "st", Sdf.ValueTypeNames.TexCoord2fArray, UsdGeom.Tokens.vertex
    )
    st.Set(
        Vt.Vec2fArray(
            [Gf.Vec2f(0, 0), Gf.Vec2f(tiles, 0), Gf.Vec2f(tiles, tiles), Gf.Vec2f(0, tiles)]
        )
    )

    material = UsdShade.Material.Define(stage, "/World/Looks/Agar")
    surface = UsdShade.Shader.Define(stage, "/World/Looks/Agar/Surface")
    surface.CreateIdAttr("UsdPreviewSurface")
    reader = UsdShade.Shader.Define(stage, "/World/Looks/Agar/St")
    reader.CreateIdAttr("UsdPrimvarReader_float2")
    reader.CreateInput("varname", Sdf.ValueTypeNames.Token).Set("st")
    reader.CreateOutput("result", Sdf.ValueTypeNames.Float2)
    image = UsdShade.Shader.Define(stage, "/World/Looks/Agar/Texture")
    image.CreateIdAttr("UsdUVTexture")
    image.CreateInput("file", Sdf.ValueTypeNames.Asset).Set(str(texture))
    image.CreateInput("st", Sdf.ValueTypeNames.Float2).ConnectToSource(
        reader.ConnectableAPI(), "result"
    )
    image.CreateInput("wrapS", Sdf.ValueTypeNames.Token).Set("repeat")
    image.CreateInput("wrapT", Sdf.ValueTypeNames.Token).Set("repeat")
    image.CreateOutput("rgb", Sdf.ValueTypeNames.Float3)
    surface.CreateInput("diffuseColor", Sdf.ValueTypeNames.Color3f).ConnectToSource(
        image.ConnectableAPI(), "rgb"
    )
    surface.CreateInput("roughness", Sdf.ValueTypeNames.Float).Set(0.6)
    surface.CreateOutput("surface", Sdf.ValueTypeNames.Token)
    material.CreateSurfaceOutput().ConnectToSource(surface.ConnectableAPI(), "surface")
    UsdShade.MaterialBindingAPI.Apply(mesh.GetPrim()).Bind(material)
    return path


def presentation_lighting(
    *,
    key_path: str = "/World/DistantLight",
    key: float = 700.0,
    fill: float = 400.0,
    softness_deg: float = 10.0,
) -> None:
    """Soft light for the skin and the agar: the scene's key light dimmed and
    widened, so its shadow has a soft edge, and a white dome filling that shadow.

    Picked by eye from screenshots of six settings. The scene's own light, 2000
    on its own, overexposes -- colours wash towards white -- and leaves a black
    shadow along every bend.
    """
    import isaacsim.core.experimental.utils.stage as stage_utils
    from pxr import Gf, UsdLux

    stage = stage_utils.get_current_stage()
    prim = stage.GetPrimAtPath(key_path)
    if prim.IsValid():
        light = UsdLux.DistantLight(prim)
        light.GetIntensityAttr().Set(key)
        light.CreateAngleAttr().Set(softness_deg)
    dome = UsdLux.DomeLight.Define(stage, "/World/FillLight")
    dome.CreateIntensityAttr().Set(fill)
    dome.CreateColorAttr().Set(Gf.Vec3f(1.0, 1.0, 1.0))


@dataclass
class ConnectomeOverlay:
    """The connectome drawn live: spheres for neurons, curves for synapses.

    ``inside`` places it in the body, so it bends with the animal; ``above`` draws
    the same, larger, floating over the agar beside the worm -- directly over it,
    the camera overhead would see the drawing and not the animal; ``screen`` is a
    flat panel in its own window. Neurons are coloured by synaptic activation,
    blue at rest through rose to gold, and each synapse by its presynaptic neuron.
    """

    mode: str
    drawing: ConnectomeDrawing
    plan: BodyPlan
    #: The neurons: one mesh holding a small sphere per cell. A ``UsdGeom.Points``
    #: drew the same, but RTX rebuilt it on every change -- updating it every
    #: frame cost 20 ms a frame, against about 1 ms for a mesh (measured).
    dots: Any = None
    #: Radius of each neuron's sphere, metres.
    dot_radius: float = 0.0
    curves: Any = None
    provider: Any = None
    window: Any = None
    panel: ConnectomePanel | None = None
    beside: Beside | None = None
    #: The panel as last drawn (``screen`` only), for :class:`ViewportCapture`.
    image: np.ndarray | None = None
    #: Half the ``above`` drawing's offset: a camera aimed this far from the worm
    #: frames both.
    view_shift: np.ndarray = field(default_factory=lambda: np.zeros(2))
    _pre: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    _unit_sphere: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    @classmethod
    def build(cls, mode: str, drawing: ConnectomeDrawing, plan: BodyPlan) -> ConnectomeOverlay:
        if mode not in CONNECTOME_MODES or mode == "off":
            raise ValueError(f"mode must be one of {CONNECTOME_MODES[1:]}, not {mode!r}")
        overlay = cls(mode=mode, drawing=drawing, plan=plan)
        order = {n: i for i, n in enumerate(drawing.names)}
        overlay._pre = np.array([order[a] for a, _, _ in drawing.edges], dtype=np.int64)
        if mode == "screen":
            import omni.ui as ui

            overlay.panel = ConnectomePanel.build(drawing.layout, drawing.names, drawing.edges)
            overlay.window = ui.Window("Connectome", width=660, height=240)
            overlay.provider = ui.ByteImageProvider()
            with overlay.window.frame:
                ui.ImageWithProvider(overlay.provider)
            return overlay

        import isaacsim.core.experimental.utils.stage as stage_utils
        from pxr import UsdGeom, Vt

        stage = stage_utils.get_current_stage()
        n_cells, n_edges = len(drawing.names), len(drawing.edges)
        inside = mode == "inside"
        if not inside:
            overlay.beside = Beside(distance_m=ABOVE_OFFSET_BODY_LENGTHS * plan.total_length_m)
        overlay.dot_radius = (0.16 if inside else 0.275) * plan.max_radius_m
        unit, counts, indices = sphere()
        overlay._unit_sphere = unit
        overlay.dots = UsdGeom.Mesh.Define(stage, "/World/Connectome/Neurons")
        overlay.dots.CreateFaceVertexCountsAttr().Set(
            Vt.IntArray(np.tile(counts, n_cells).tolist())
        )
        overlay.dots.CreateFaceVertexIndicesAttr().Set(
            Vt.IntArray(
                (indices[None, :] + len(unit) * np.arange(n_cells)[:, None]).ravel().tolist()
            )
        )
        overlay.dots.CreateSubdivisionSchemeAttr().Set(UsdGeom.Tokens.none)
        # A sphere's normals are its own unit vertices, wherever it is moved to.
        overlay.dots.CreateNormalsAttr().Set(
            Vt.Vec3fArray.FromNumpy(np.tile(unit, (n_cells, 1)).astype(np.float32))
        )
        overlay.dots.SetNormalsInterpolation(UsdGeom.Tokens.vertex)
        overlay.dots.CreateDisplayColorPrimvar(UsdGeom.Tokens.vertex)
        overlay.curves = UsdGeom.BasisCurves.Define(stage, "/World/Connectome/Synapses")
        overlay.curves.CreateTypeAttr().Set(UsdGeom.Tokens.cubic)
        overlay.curves.CreateBasisAttr().Set(UsdGeom.Tokens.catmullRom)
        overlay.curves.CreateWrapAttr().Set(UsdGeom.Tokens.nonperiodic)
        overlay.curves.CreateCurveVertexCountsAttr().Set(
            Vt.IntArray([drawing.points_per_edge] * n_edges)
        )
        overlay.curves.CreateWidthsAttr().Set(
            Vt.FloatArray([(0.05 if inside else 0.1) * plan.max_radius_m])
        )
        overlay.curves.SetWidthsInterpolation(UsdGeom.Tokens.constant)
        overlay.curves.CreateDisplayColorPrimvar(UsdGeom.Tokens.uniform)
        # Complete -- positions, extents, colours -- before the first frame draws it.
        at_rest = np.full(int(drawing.indices.max()) + 1, drawing.rest)
        overlay.update(_straight_centres(plan, 0.0), at_rest, 0.0)
        return overlay

    def update(self, centres: np.ndarray, activation: np.ndarray, t: float = 0.0) -> None:
        """``activation`` is the network's synaptic activation for every cell, and
        ``t`` the simulated time in seconds, which ``above`` smooths its offset over."""
        from pxr import Vt

        level = self.drawing.level(activation)
        if self.mode == "screen" and self.panel is not None:
            self.image = self.panel.draw(level)
            self.provider.set_data_array(self.image, [self.image.shape[1], self.image.shape[0]])
            return
        inside = self.mode == "inside"
        # Within the outline, over the skin; or larger, and lifted off the agar.
        scale, lift_m = (0.8, 0.0) if inside else (2.2, 5.0 * self.plan.max_radius_m)
        xyz, paths = self.drawing.positions(
            centres, self.plan, scale=scale, lift_m=lift_m, on_surface=inside
        )
        if inside:
            xyz[:, 2] += 0.02 * self.plan.max_radius_m  # dots over their own curves
        elif self.beside is not None:
            shift = self.beside.offset(centres, t)
            xyz += shift
            paths += shift
            self.view_shift = 0.5 * shift[:2]
        colours = activity_colours(level)
        unit = self._unit_sphere
        dots = (xyz[:, None, :] + self.dot_radius * unit[None, :, :]).reshape(-1, 3)
        self.dots.GetPointsAttr().Set(Vt.Vec3fArray.FromNumpy(dots.astype(np.float32)))
        self.dots.GetDisplayColorPrimvar().Set(
            Vt.Vec3fArray.FromNumpy(
                np.repeat(linear(colours), len(unit), axis=0).astype(np.float32)
            )
        )
        self.curves.GetPointsAttr().Set(Vt.Vec3fArray.FromNumpy(paths.astype(np.float32)))
        # Each synapse in its presynaptic cell's colour, a shade darker than the cells.
        self.curves.GetDisplayColorPrimvar().Set(
            Vt.Vec3fArray.FromNumpy(linear(0.8 * colours[self._pre]).astype(np.float32))
        )
        both = np.vstack([dots, paths]) if len(paths) else dots
        extent = Vt.Vec3fArray.FromNumpy(
            np.vstack([both.min(axis=0), both.max(axis=0)]).astype(np.float32)
        )
        self.dots.GetExtentAttr().Set(extent)
        self.curves.GetExtentAttr().Set(extent)


@dataclass
class ViewportCapture:
    """Saves the viewport to ``folder`` every ``every_s`` seconds of simulated time.

    The first frame comes one interval in: at t = 0 the stage has not finished
    loading and the capture is of a grey void. A viewport capture holds the
    viewport only, so the on-screen connectome panel, when given, is saved beside
    each frame as ``panel_*.png``.
    """

    folder: Path
    every_s: float
    _last: float = 0.0
    count: int = 0

    def maybe(self, t: float, panel: np.ndarray | None = None) -> None:
        if self.every_s <= 0.0 or t - self._last < self.every_s:
            return
        from omni.kit.viewport.utility import capture_viewport_to_file, get_active_viewport

        viewport = get_active_viewport()
        if viewport is None:
            return
        self.folder.mkdir(parents=True, exist_ok=True)
        capture_viewport_to_file(viewport, file_path=str(self.folder / f"frame_{t:07.2f}s.png"))
        if panel is not None:
            (self.folder / f"panel_{t:07.2f}s.png").write_bytes(encode_png(panel))
        self._last = t
        self.count += 1
