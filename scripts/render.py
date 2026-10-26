"""Render a region in Blender from the viewer's data (viewer/data/<region>/), headless.

blender -b -P scripts/render.py -- <region>    # figures/render-<region>.jpg, about 30 s for Carina on an M4
    [--elev 15] [--azim 60] [--focus KPC] [--dist KPC] [--aspect 2.4] [--width 2400] [--samples 64] [--out x.jpg]

Every point is a star in the export: field stars at their Bailer-Jones distance, the members of each
recovered cluster at the cluster's distance, coloured by Hunt & Reffert's age as the viewer's age bands.
Missed clusters are red rings at their catalogue position. Points and rings are sized to be seen, not to
scale. By default the camera looks along the wedge at the recovered clusters' median distance, from as
far away, in perspective; `--dist 0` frames the whole wedge in an orthographic view. One Blender unit is 100 pc; the Sun is at the origin, x towards the Galactic centre, z to the
north Galactic pole.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import bpy
import numpy as np
from mathutils import Matrix

UNIT = 100.0  # pc per Blender unit
PAPER, INK, RULE, MISSED = "#f5f5f2", "#222222", "#c4c5bf", "#c8553d"
AGES = [(7, "#86b6ef"), (8, "#3987e5"), (9, "#1c5cab"), (math.inf, "#0d366b")]  # the viewer's bands, by log10(age / yr)

args = argparse.ArgumentParser()
args.add_argument("region")
args.add_argument("--elev", type=float, default=15, help="degrees above the Galactic plane; 90 is from above")
args.add_argument("--azim", type=float, default=60, help="degrees from straight across towards behind the Sun")
args.add_argument("--aspect", type=float, default=2.4)
args.add_argument("--dist", type=float, help="kpc from the camera to its focus (default: as far as the focus "
                  "from the Sun); 0 for an orthographic view of the whole wedge")
args.add_argument("--focus", type=float, help="kpc from the Sun along the wedge to look at (default: the recovered "
                  "clusters' median distance)")
args.add_argument("--lens", type=float, default=50)
args.add_argument("--fstop", type=float, default=0.5, help="in Blender units: smaller blurs more")
args.add_argument("--width", type=int, default=2400)
args.add_argument("--samples", type=int, default=64)
args.add_argument("--out")
a = args.parse_args(sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else [])
root = Path(__file__).resolve().parent.parent
data = root / "viewer" / "data" / a.region
meta = json.loads((data / "meta.json").read_text())
out = Path(a.out) if a.out else root / "figures" / f"render-{a.region}.jpg"


def linear(hex_):
    """sRGB hex to the linear RGBA Blender's colours take."""
    c = [int(hex_[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    return tuple(v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in c) + (1.0,)


def toward(l, b, d):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)]) * d


# Stars
lo, hi = np.array(meta["quant"]["lo"]), np.array(meta["quant"]["hi"])
field = np.vstack([np.frombuffer((data / lv["file"]).read_bytes(), "<u2").reshape(-1, 3) for lv in meta["field"]["levels"]])
field = (lo + field * (hi - lo) / 65535) / UNIT
m = meta["members"]["count"]
raw = (data / meta["members"]["file"]).read_bytes()
members = np.frombuffer(raw, "<f4", 3 * m).reshape(-1, 3) / UNIT
owner = np.frombuffer(raw, "<u2", m, offset=12 * m)
band = np.array([next(hex_ for top, hex_ in AGES if c["age"][1] < top) for c in meta["clusters"]])
missed = np.array([c["xyz"] for c in meta["clusters"] if not c["recovered"]]).reshape(-1, 3) / UNIT

# Scene
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.render.engine = "CYCLES"
cycles = bpy.context.preferences.addons["cycles"].preferences
cycles.compute_device_type = "METAL"
cycles.get_devices()
for device in cycles.devices:
    device.use = True
scene.cycles.device = "GPU"
scene.cycles.samples = a.samples
scene.cycles.use_denoising = False  # the denoiser smears single stars
scene.cycles.max_bounces = 0
scene.cycles.transparent_max_bounces = 1024  # the field near the Sun stacks hundreds of stars deep
scene.view_settings.view_transform = "Standard"  # hex in, hex out, as the viewer
world = bpy.data.worlds.new("paper")
world.use_nodes = True
world.node_tree.nodes["Background"].inputs["Color"].default_value = linear(PAPER)
scene.world = world


def flat(name, colour, alpha=1.0, attribute=None):
    """A material that paints its colour over what is behind it at `alpha`, unlit: alpha compositing, so many
    faint field stars add up to a density, as the viewer's log-scaled grey. With `attribute`, the colour is
    the point's own."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    emit, clear, mix, outp = (nodes.new(t) for t in ("ShaderNodeEmission", "ShaderNodeBsdfTransparent",
                                                      "ShaderNodeMixShader", "ShaderNodeOutputMaterial"))
    emit.inputs["Color"].default_value = linear(colour)
    if attribute:
        attr = nodes.new("ShaderNodeAttribute")
        attr.attribute_name = attribute
        links.new(attr.outputs["Color"], emit.inputs["Color"])
    mix.inputs["Fac"].default_value = alpha
    links.new(clear.outputs[0], mix.inputs[1])
    links.new(emit.outputs[0], mix.inputs[2])
    links.new(mix.outputs[0], outp.inputs["Surface"])
    return mat


def cloud(name, xyz, radius, material, colours=None):
    """Points, each a small sphere to Cycles."""
    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(len(xyz))
    mesh.vertices.foreach_set("co", np.ascontiguousarray(xyz, np.float32).ravel())
    if colours is not None:
        attr = mesh.attributes.new("colour", "FLOAT_COLOR", "POINT")
        attr.data.foreach_set("color", np.ascontiguousarray(colours, np.float32).ravel())
    obj = bpy.data.objects.new(name, mesh)
    scene.collection.objects.link(obj)
    tree = bpy.data.node_groups.new(name, "GeometryNodeTree")
    tree.interface.new_socket("Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
    tree.interface.new_socket("Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
    nodes, links = tree.nodes, tree.links
    gin, gout = nodes.new("NodeGroupInput"), nodes.new("NodeGroupOutput")
    points, setmat = nodes.new("GeometryNodeMeshToPoints"), nodes.new("GeometryNodeSetMaterial")
    points.inputs["Radius"].default_value = radius
    setmat.inputs["Material"].default_value = material
    links.new(gin.outputs[0], points.inputs["Mesh"])
    links.new(points.outputs[0], setmat.inputs["Geometry"])
    links.new(setmat.outputs[0], gout.inputs[0])
    obj.modifiers.new(name, "NODES").node_group = tree
    return obj


def rings(name, xyz, radius, thickness, material):
    """A ring lying in the Galactic plane at each point."""
    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(len(xyz))
    mesh.vertices.foreach_set("co", np.ascontiguousarray(xyz, np.float32).ravel())
    obj = bpy.data.objects.new(name, mesh)
    scene.collection.objects.link(obj)
    tree = bpy.data.node_groups.new(name, "GeometryNodeTree")
    tree.interface.new_socket("Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
    tree.interface.new_socket("Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
    nodes, links = tree.nodes, tree.links
    gin, gout = nodes.new("NodeGroupInput"), nodes.new("NodeGroupOutput")
    path, profile = nodes.new("GeometryNodeCurvePrimitiveCircle"), nodes.new("GeometryNodeCurvePrimitiveCircle")
    path.inputs["Radius"].default_value, profile.inputs["Radius"].default_value = radius, thickness
    profile.inputs["Resolution"].default_value = 6
    path.inputs["Resolution"].default_value = 48
    tube, place = nodes.new("GeometryNodeCurveToMesh"), nodes.new("GeometryNodeInstanceOnPoints")
    setmat = nodes.new("GeometryNodeSetMaterial")
    setmat.inputs["Material"].default_value = material
    links.new(path.outputs["Curve"], tube.inputs["Curve"])
    links.new(profile.outputs["Curve"], tube.inputs["Profile Curve"])
    links.new(gin.outputs[0], place.inputs["Points"])
    links.new(tube.outputs[0], place.inputs["Instance"])
    links.new(place.outputs[0], setmat.inputs["Geometry"])
    links.new(setmat.outputs[0], gout.inputs[0])
    obj.modifiers.new(name, "NODES").node_group = tree
    return obj


def guides(far, thickness, material):
    """The Galactic plane under the wedge: its two edges from the Sun and an arc every kiloparsec."""
    (l0, l1) = meta["box"]["l"]
    curve = bpy.data.curves.new("guides", "CURVE")
    curve.dimensions, curve.bevel_depth, curve.bevel_resolution = "3D", thickness, 2
    polylines = [[toward(l, 0, 0), toward(l, 0, far)] for l in (l0, l1)]
    polylines += [[toward(l0 + (l1 - l0) * s / 40, 0, d) for s in range(41)] for d in np.arange(10, far + 1e-6, 10)]
    for line in polylines:
        spline = curve.splines.new("POLY")
        spline.points.add(len(line) - 1)
        for p, co in zip(spline.points, line):
            p.co = (*co, 1.0)
    obj = bpy.data.objects.new("guides", curve)
    obj.data.materials.append(material)
    scene.collection.objects.link(obj)


far = math.ceil(max(meta["reach"], max(np.linalg.norm(c["xyz"]) for c in meta["clusters"])) / 1000) * 10
cloud("field", field, 0.012, flat("field", INK, 0.05))
cloud("members", members, 0.03, flat("members", INK, attribute="colour"),
      colours=np.array([linear(h) for h in band])[owner])
rings("missed", missed, 0.25, 0.018, flat("missed", MISSED))
guides(far, 0.006, flat("rule", RULE))

# Camera: it looks at the wedge from the side the Sun is on the left of, `elev` degrees above the Galactic
# plane and turned `azim` degrees from straight across towards the Sun's back. In perspective it looks at
# the wedge's axis `focus` kpc out; orthographic (`--dist 0`), it frames the whole wedge with a margin
if a.focus is None:
    a.focus = float(np.median([c["dist"] for c in meta["clusters"] if c["recovered"]])) / 1000
if a.dist is None:
    a.dist = a.focus
lc = sum(meta["box"]["l"]) / 2
axis, side, north = toward(lc, 0, 1), toward(lc + 90, 0, 1), np.array([0, 0, 1.0])
elev, azim = np.radians(a.elev), np.radians(a.azim)
back = np.cos(elev) * (-np.cos(azim) * side - np.sin(azim) * axis) + np.sin(elev) * north  # from the wedge to the camera
right = np.cross(-back, side if a.elev > 80 else north)
right /= np.linalg.norm(right)
up = np.cross(right, -back)
(l0, l1), (b0, b1) = meta["box"]["l"], meta["box"]["b"]
hull = np.array([toward(l, b, far) for l in np.linspace(l0, l1, 11) for b in (b0, b1)] + [[0, 0, 0]])
x, y = hull @ right, hull @ up
aspect = a.aspect
centre = right * (x.min() + x.max()) / 2 + up * (y.min() + y.max()) / 2
cam = bpy.data.objects.new("camera", bpy.data.cameras.new("camera"))
scene.collection.objects.link(cam)
scene.camera = cam
if a.dist:  # in perspective, from `dist` kpc away, looking at the axis `focus` kpc from the Sun, sharp there
    centre = axis * a.focus * 10
    location = centre + back * a.dist * 10
    cam.data.lens = a.lens
    cam.data.dof.use_dof, cam.data.dof.focus_distance, cam.data.dof.aperture_fstop = True, a.dist * 10, a.fstop
else:
    cam.data.type = "ORTHO"
    cam.data.ortho_scale = 1.08 * max(np.ptp(x), np.ptp(y) * aspect)
    location = centre + back * far * 3
cam.data.clip_end = far * 10
cam.matrix_world = Matrix([[*right, 0], [*up, 0], [*back, 0], [*location, 1]]).transposed()

scene.render.resolution_x = a.width
scene.render.resolution_y = round(a.width / aspect)
scene.render.filepath = str(out)
scene.render.image_settings.file_format = "JPEG" if out.suffix.lower() in (".jpg", ".jpeg") else "PNG"
scene.render.image_settings.quality = 92
out.parent.mkdir(exist_ok=True)
bpy.ops.render.render(write_still=True)
print(f"{a.region}: {len(field):,} field stars, {m:,} members, {len(missed)} missed; wrote {out}")
