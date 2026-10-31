"""A one-minute film of a region in Blender, from the viewer's data, night style; headless.

blender -b -P scripts/film.py -- <region>                 # frames into data/film/<region>/, 1920x1080, 30 fps
    [--preview] [--start F] [--end F] [--samples 32]     # --preview: 640x360, 8 samples, every third frame
uv run scripts/film_encode.py <region>                   # the frames, with titles, into figures/film-<region>.mp4

The same stars as scripts/render.py, on a night sky: the field as faint dust, recovered clusters glowing in
Hunt & Reffert's age bands, missed clusters as red rings. Stars drift as in the viewer's time slider: along
the sky at their proper motion, less the field's typical motion at their distance (meta.reference), never in
depth. Members get a modelled depth, as deep as their cluster is wide on the sky (the export has none). The
shots are in SHOTS; the titles, in film_encode.py, follow the same seconds. Frames already
rendered are kept, so a stopped run resumes. One Blender unit is 100 pc; the Sun is at the origin.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import bpy
import numpy as np

UNIT = 100.0  # pc per Blender unit
FPS = 30
MAS_MYR = 4.8481e-3  # radians crossed in a million years at 1 mas/yr, as the viewer
MYR = 3.0  # how far the time runs, in Myr
SKY, DUST, MISSED = "#03050a", "#7d8cb3", "#ff6a45"
AGES = [(7, "#d8ecff"), (8, "#9fd0ff"), (9, "#6fb2ff"), (math.inf, "#4f8fff")]  # the viewer's bands, lightened

args = argparse.ArgumentParser()
args.add_argument("region")
args.add_argument("--preview", action="store_true")
args.add_argument("--start", type=int, default=0)
args.add_argument("--end", type=int)
args.add_argument("--samples", type=int, default=32)
a = args.parse_args(sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else [])
root = Path(__file__).resolve().parent.parent
data = root / "viewer" / "data" / a.region
meta = json.loads((data / "meta.json").read_text())
frames = root / "data" / "film" / (a.region + ("-preview" if a.preview else ""))


def linear(hex_):
    """sRGB hex to the linear RGBA Blender's colours take."""
    c = [int(hex_[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    return tuple(v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in c) + (1.0,)


def toward(l, b, d):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)]) * d


def velocity(xyz_pc, pm):
    """Blender units per Myr: the viewer's drift(), the field's typical motion at each distance taken off."""
    d = np.linalg.norm(xyz_pc, axis=1)
    ref, step = np.array(meta["reference"]["pm"]), meta["reference"]["step"]
    x = np.clip(d / step - 0.5, 0, len(ref) - 1)
    i, f = x.astype(int), (x % 1)[:, None]
    typical = ref[i] * (1 - f) + ref[np.minimum(i + 1, len(ref) - 1)] * f
    r = xyz_pc / np.maximum(d, 1e-3)[:, None]
    along = np.column_stack([-r[:, 1], r[:, 0], np.zeros(len(r))])
    along /= np.linalg.norm(along, axis=1)[:, None]
    across = np.cross(r, along)
    k = (d * MAS_MYR / UNIT)[:, None]
    rel = pm - typical
    return k * (along * rel[:, :1] + across * rel[:, 1:])


# Stars
lo, hi = np.array(meta["quant"]["lo"]), np.array(meta["quant"]["hi"])
levels = meta["field"]["levels"]
field = np.vstack([np.frombuffer((data / lv["file"]).read_bytes(), "<u2").reshape(-1, 3) for lv in levels])
field = lo + field * (hi - lo) / 65535
field_pm = np.vstack([np.frombuffer((data / lv["pm"]).read_bytes(), "<i2").reshape(-1, 2) for lv in levels])
m = meta["members"]["count"]
raw = (data / meta["members"]["file"]).read_bytes()
members = np.frombuffer(raw, "<f4", 3 * m).reshape(-1, 3).astype(float)
owner = np.frombuffer(raw, "<u2", m, offset=12 * m)
# The export puts every member at its cluster's distance, so a cluster is a disc facing the Sun and, from the
# side, a line. The film gives each one a depth as deep as it is wide on the sky (clusters are roughly round):
# a modelled depth, not a measured one, for a camera that leaves the Sun's line of sight
rng = np.random.default_rng(0)
for k in np.unique(owner):
    mine = owner == k
    centre = members[mine].mean(axis=0)
    r = centre / np.linalg.norm(centre)
    offsets = members[mine] - centre
    offsets -= np.outer(offsets @ r, r)
    sigma = np.sqrt((offsets ** 2).sum(axis=1).mean() / 2)
    members[mine] += np.outer(rng.normal(0, sigma, mine.sum()), r)
members_pm = np.frombuffer((data / meta["members"]["pm"]).read_bytes(), "<i2").reshape(-1, 2)
step = meta["pm_step"]
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
scene.cycles.samples = 8 if a.preview else a.samples
scene.cycles.use_denoising = False  # the denoiser smears single stars
scene.cycles.max_bounces = 0
scene.cycles.transparent_max_bounces = 1024
scene.view_settings.view_transform = "Standard"
world = bpy.data.worlds.new("sky")
world.use_nodes = True
world.node_tree.nodes["Background"].inputs["Color"].default_value = linear(SKY)
scene.world = world


def key(socket, frame, value):
    socket.default_value = value
    socket.keyframe_insert("default_value", frame=frame)


def glowing(name, colour, attribute=None):
    """Unlit emission over what is behind it. Returns the material and its mix (alpha), colour mix and strength
    sockets, to animate: with `attribute`, `shift` goes from `colour` (0) to the point's own colour (1)."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    emit, clear, mix, outp = (nodes.new(t) for t in ("ShaderNodeEmission", "ShaderNodeBsdfTransparent",
                                                      "ShaderNodeMixShader", "ShaderNodeOutputMaterial"))
    shift = None
    if attribute:
        attr, blend = nodes.new("ShaderNodeAttribute"), nodes.new("ShaderNodeMix")
        attr.attribute_name = attribute
        blend.data_type = "RGBA"
        blend.inputs["A"].default_value = linear(colour)
        links.new(attr.outputs["Color"], blend.inputs["B"])
        links.new(blend.outputs["Result"], emit.inputs["Color"])
        shift = blend.inputs["Factor"]
    else:
        emit.inputs["Color"].default_value = linear(colour)
    links.new(clear.outputs[0], mix.inputs[1])
    links.new(emit.outputs[0], mix.inputs[2])
    links.new(mix.outputs[0], outp.inputs["Surface"])
    return mat, mix.inputs["Fac"], shift, emit.inputs["Strength"]


clocks = []  # every cloud's time, in Myr: keyed together


def cloud(name, xyz, vel, radius, material, colours=None):
    """Points, each a small sphere to Cycles, moved by `vel` (Blender units per Myr) times the clock."""
    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(len(xyz))
    mesh.vertices.foreach_set("co", np.ascontiguousarray(xyz, np.float32).ravel())
    v = mesh.attributes.new("vel", "FLOAT_VECTOR", "POINT")
    v.data.foreach_set("vector", np.ascontiguousarray(vel, np.float32).ravel())
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
    named, scale, clock = nodes.new("GeometryNodeInputNamedAttribute"), nodes.new("ShaderNodeVectorMath"), nodes.new("ShaderNodeValue")
    named.data_type = "FLOAT_VECTOR"
    named.inputs["Name"].default_value = "vel"
    scale.operation = "SCALE"
    move, points, setmat = (nodes.new(t) for t in ("GeometryNodeSetPosition", "GeometryNodeMeshToPoints",
                                                   "GeometryNodeSetMaterial"))
    points.inputs["Radius"].default_value = radius
    setmat.inputs["Material"].default_value = material
    links.new(named.outputs["Attribute"], scale.inputs[0])
    links.new(clock.outputs[0], scale.inputs["Scale"])
    links.new(scale.outputs["Vector"], move.inputs["Offset"])
    links.new(gin.outputs[0], move.inputs["Geometry"])
    links.new(move.outputs[0], points.inputs["Mesh"])
    links.new(points.outputs[0], setmat.inputs["Geometry"])
    links.new(setmat.outputs[0], gout.inputs[0])
    obj.modifiers.new(name, "NODES").node_group = tree
    clocks.append(clock.outputs[0])
    return points.inputs["Radius"]


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


DUST_ALPHA = 0.022
dust, dust_alpha, _, _ = glowing("field", DUST)
dust_radius = cloud("field", field / UNIT, velocity(field, field_pm * step), 0.012, dust)
lit, lit_alpha, lit_shift, lit_strength = glowing("members", DUST, attribute="colour")
lit_radius = cloud("members", members / UNIT, velocity(members, members_pm * step), 0.03, lit,
      colours=np.array([linear(h) for h in band])[owner])
red, red_alpha, _, red_strength = glowing("missed", MISSED)
red_strength.default_value = 3.0
rings("missed", missed, 0.25, 0.018, red)

# Bloom on what is brighter than the dust: the lit clusters and the rings
scene.use_nodes = True
tree = getattr(scene, "node_tree", None)
if tree is None:  # Blender 5: the compositor is a node group of its own
    tree = bpy.data.node_groups.new("compositor", "CompositorNodeTree")
    scene.compositing_node_group = tree
    tree.interface.new_socket("Image", in_out="OUTPUT", socket_type="NodeSocketColor")
    layers, out = tree.nodes.new("CompositorNodeRLayers"), tree.nodes.new("NodeGroupOutput")
else:
    layers, out = tree.nodes["Render Layers"], tree.nodes["Composite"]
glare = tree.nodes.new("CompositorNodeGlare")
for name, value in (("Type", "Bloom"), ("Threshold", 0.85), ("Strength", 1.0), ("Size", 0.7)):
    if name in glare.inputs:
        glare.inputs[name].default_value = value
tree.links.new(layers.outputs["Image"], glare.inputs["Image"])
tree.links.new(glare.outputs[0], out.inputs[0])

# Shots, in seconds. Points the camera passes through or looks at, along the wedge's axis (`axis`), across it
# (`side`, the Sun on the left) and up (z), in Blender units (100 pc).
lc = sum(meta["box"]["l"]) / 2
axis, side, up = toward(lc, 0, 1), toward(lc - 90, 0, 1), np.array([0, 0, 1.0])
focus = float(np.median([c["dist"] for c in meta["clusters"] if c["recovered"]])) / UNIT
big = max((c for c in meta["clusters"] if c["recovered"] and c["dist"] < 1000), key=lambda c: c["n"])
near = np.array(big["xyz"]) / UNIT
at = lambda along, across=0.0, height=0.0: axis * along + side * across + up * height
SHOTS = [  # (second, camera, target)
    (0, at(-1.5, 0, 0.5), at(10)),
    (6, at(-0.5, 0.3, 0.6), at(12)),
    (14, at(8, 1.5, 1.0), at(focus)),
    (20, at(focus - 9, 7, 5), at(focus)),
    (36, at(focus - 4, 9, 6), at(focus + 1)),
    (44, at(focus - 6, 24, 20), at(focus)),
    (47, at(focus - 6, 24, 20), at(focus)),
    (54, near + side * 0.9 + up * 0.35 - axis * 0.6, near),
    (60, near + side * 0.7 + up * 0.25 - axis * 0.4, near),
]
cam = bpy.data.objects.new("camera", bpy.data.cameras.new("camera"))
target = bpy.data.objects.new("target", None)
for obj in (cam, target):
    scene.collection.objects.link(obj)
scene.camera = cam
cam.data.lens, cam.data.clip_start, cam.data.clip_end = 35, 0.01, 1000
track = cam.constraints.new("TRACK_TO")
track.target, track.track_axis, track.up_axis = target, "TRACK_NEGATIVE_Z", "UP_Y"
for second, where, look in SHOTS:
    cam.location, target.location = where, look
    cam.keyframe_insert("location", frame=second * FPS)
    target.keyframe_insert("location", frame=second * FPS)

# The dust fades in; at 21 s the clock starts and the clusters light up; at 36 s it runs back to now and the
# missed clusters show.
s = lambda second: round(second * FPS)
for socket, steps in ((dust_alpha, [(0, 0.0), (4, DUST_ALPHA)]),
                      (lit_alpha, [(0, 0.0), (4, DUST_ALPHA), (23, DUST_ALPHA), (29, 1.0)]),
                      (lit_shift, [(0, 0.0), (23, 0.0), (29, 1.0)]),
                      (lit_strength, [(0, 1.0), (23, 1.0), (29, 3.0)]),
                      (red_alpha, [(0, 0.0), (40, 0.0), (42, 1.0), (47, 1.0), (49, 0.0)]),
                      # sized to be seen from kiloparsecs; in the close-up, a star is a star
                      (dust_radius, [(0, 0.012), (48, 0.012), (54, 0.0015)]),
                      (lit_radius, [(0, 0.03), (48, 0.03), (54, 0.004)])):
    for second, value in steps:
        key(socket, s(second), value)
for clock in clocks:
    for second, value in ((0, 0.0), (21, 0.0), (33, MYR), (36, MYR), (40, 0.0)):
        key(clock, s(second), value)

scene.render.fps = FPS
scene.frame_start, scene.frame_end = a.start, a.end if a.end is not None else 60 * FPS - 1
scene.frame_step = 3 if a.preview else 1
scene.render.resolution_x, scene.render.resolution_y = (640, 360) if a.preview else (1920, 1080)
scene.render.image_settings.file_format = "PNG"
scene.render.use_overwrite = False  # resume a stopped run
frames.mkdir(parents=True, exist_ok=True)
scene.render.filepath = str(frames / "####")
print(f"{a.region}: {len(field):,} field stars, {m:,} members, {len(missed)} missed; close-up on {big['name']}")
bpy.ops.render.render(animation=True)
