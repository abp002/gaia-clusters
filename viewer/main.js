// The gaia-clusters viewer: one scored region in three dimensions.
// The files come from scripts/export_viewer.py; their layout is in src/gaia_clusters/viewer.py.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

THREE.ColorManagement.enabled = false;  // hex in, hex out: the colours of the README figures

const HEX = { paper: 0xf5f5f2, ink: 0x222222, recovered: 0x1f3b73, missed: 0xc8553d, field: 0x9a9a9a,
              rule: 0xc4c5bf, unverified: 0x6b6b67 };
const coarse = matchMedia('(pointer: coarse)').matches;
const BUDGET = coarse ? 600_000 : Infinity;  // field stars loaded before asking
const MOVING = coarse ? 250_000 : 700_000;   // field stars drawn while the view moves
const DARKEST = 400;                         // field stars per pixel drawn in the darkest grey
const MIN_RADIUS = 4;                        // CSS pixels: no cluster is drawn smaller than this
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;
const RAD = Math.PI / 180;
const num = new Intl.NumberFormat('en-GB').format;
const pct = (x) => `${Math.round(100 * x)} %`;
const kpc = (pc) => `${(pc / 1000).toFixed(pc < 1000 ? 2 : 1)} kpc`;
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const $ = (s) => document.querySelector(s);
const toward = (l, b, d) => new THREE.Vector3(Math.cos(b * RAD) * Math.cos(l * RAD), Math.cos(b * RAD) * Math.sin(l * RAD),
                                              Math.sin(b * RAD)).multiplyScalar(d);

const box = $('#view-box'), canvas = $('#view'), labelsEl = $('#labels'), card = $('#card');
const narrow = matchMedia('(max-width: 900px)');  // the CSS breakpoint where the view gets tall
let renderer;
try {
  renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
} catch {
  box.innerHTML = '<p class="hint mono" style="position:static;padding:24px">This page needs WebGL, which this browser has turned off.</p>';
  throw new Error('no WebGL');
}
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
renderer.setClearColor(HEX.paper);
const dpr = renderer.getPixelRatio();
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(35, 1, 0.5, 100000);
camera.up.set(0, 0, 1);  // north Galactic pole up
const controls = new OrbitControls(camera, canvas);
Object.assign(controls, { enableDamping: true, dampingFactor: 0.12, zoomToCursor: true, screenSpacePanning: true,
                          minDistance: 2, maxDistance: 50000 });

// The field is drawn as a density, like the README's hexbin maps: each star adds one to its pixel in an
// off-screen buffer, and the count is shown in grey on a log scale. Field stars arrive as 16-bit steps
// between the region's corners (quant.lo, quant.hi).
const quant = { lo: { value: new THREE.Vector3() }, step: { value: new THREE.Vector3() } };
const fieldScene = new THREE.Scene();
// Half floats count far past the 255 of a byte; a browser without them gets the byte, which clips
const halfFloat = ['EXT_color_buffer_float', 'EXT_color_buffer_half_float'].some((e) => renderer.extensions.has(e));
const unit = halfFloat ? 1 : 1 / 255;
const density = new THREE.WebGLRenderTarget(1, 1, { depthBuffer: false, type: halfFloat ? THREE.HalfFloatType : THREE.UnsignedByteType,
                                                    minFilter: THREE.NearestFilter, magFilter: THREE.NearestFilter });
const fieldMaterial = new THREE.ShaderMaterial({
  uniforms: { ...quant, size: { value: dpr } },  // one CSS pixel, so the counts do not depend on the screen
  vertexShader: `uniform vec3 lo; uniform vec3 step; uniform float size;
    void main() { gl_Position = projectionMatrix * modelViewMatrix * vec4(lo + position * step, 1.0); gl_PointSize = size; }`,
  fragmentShader: `void main() { gl_FragColor = vec4(${unit.toFixed(8)}); }`,
  blending: THREE.CustomBlending, blendEquation: THREE.AddEquation, blendSrc: THREE.OneFactor, blendDst: THREE.OneFactor,
  transparent: true, depthTest: false, depthWrite: false,
});
const tone = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), new THREE.ShaderMaterial({
  uniforms: { density: { value: density.texture }, size: { value: new THREE.Vector2(1, 1) }, gain: { value: 1 },
              paper: { value: new THREE.Color(HEX.paper) }, dark: { value: new THREE.Color(HEX.field) } },
  vertexShader: 'void main() { gl_Position = vec4(position.xy, 0.0, 1.0); }',
  fragmentShader: `uniform sampler2D density; uniform vec2 size; uniform float gain; uniform vec3 paper; uniform vec3 dark;
    void main() { float n = texture2D(density, gl_FragCoord.xy / size).r * ${(1 / unit).toFixed(1)} * gain;
      gl_FragColor = vec4(mix(paper, dark, clamp(log(1.0 + n) / log(1.0 + ${DARKEST.toFixed(1)}), 0.0, 1.0)), 1.0); }`,
  depthTest: false, depthWrite: false,
}));
tone.frustumCulled = false;
tone.renderOrder = -1;
scene.add(tone);
// Drawing-buffer pixels per parsec at one parsec from the camera: for things drawn at their size in parsecs
const proj = { value: 1 };
// `hi` is the highlighted cluster or group, drawn in ink. A recovered cluster shows its members only then,
// in place of its sphere
const memberMaterial = new THREE.ShaderMaterial({
  uniforms: { size: { value: 3.2 * dpr }, ink: { value: new THREE.Color(HEX.ink) }, hi: { value: -1 } },
  vertexShader: `uniform float size; uniform float hi; attribute float cluster;
    void main() { gl_Position = abs(cluster - hi) < 0.5 ? projectionMatrix * modelViewMatrix * vec4(position, 1.0) : vec4(2.0, 2.0, 2.0, 1.0);
      gl_PointSize = size; }`,
  fragmentShader: `uniform vec3 ink;
    void main() { vec2 c = gl_PointCoord - 0.5; if (dot(c, c) > 0.25) discard; gl_FragColor = vec4(ink, 1.0); }`,
  transparent: true, depthTest: false,
});
const sphereMaterial = new THREE.MeshLambertMaterial({ color: HEX.recovered });
const ringMaterial = new THREE.ShaderMaterial({
  uniforms: { color: { value: new THREE.Color(HEX.missed) }, ink: { value: new THREE.Color(HEX.ink) }, hi: { value: -1 },
              line: { value: 1.2 * dpr }, scale: { value: dpr }, proj, least: { value: 2 * MIN_RADIUS * dpr } },
  vertexShader: `uniform float hi; uniform float scale; uniform float proj; uniform float least; attribute float radius; attribute float idx;
    varying float vHi; varying float vSize;
    void main() { vHi = abs(idx - hi) < 0.5 ? 1.0 : 0.0; vec4 p = modelViewMatrix * vec4(position, 1.0);
      vSize = max(least, 2.0 * radius * proj / -p.z) + 2.0 * scale * vHi;
      gl_Position = projectionMatrix * p; gl_PointSize = vSize; }`,
  fragmentShader: `uniform vec3 color; uniform vec3 ink; uniform float line; varying float vHi; varying float vSize;
    void main() { float r = length(gl_PointCoord - 0.5) * vSize, w = line * (1.0 + vHi);
      if (r > 0.5 * vSize || r < 0.5 * vSize - w) discard; gl_FragColor = vec4(mix(color, ink, vHi), 1.0); }`,
  transparent: true,
});
const squareMaterial = new THREE.ShaderMaterial({
  uniforms: { color: { value: new THREE.Color(HEX.unverified) }, ink: { value: new THREE.Color(HEX.ink) }, hi: { value: -1 },
              size: { value: 7 * dpr }, line: { value: 1.3 * dpr } },
  vertexShader: `uniform float size; uniform float hi; attribute float idx; varying float vHi;
    void main() { vHi = abs(idx - hi) < 0.5 ? 1.0 : 0.0;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); gl_PointSize = size * (1.0 + 0.4 * vHi); }`,
  fragmentShader: `uniform vec3 color; uniform vec3 ink; uniform float size; uniform float line; varying float vHi;
    void main() { vec2 c = abs(gl_PointCoord - 0.5); float m = max(c.x, c.y);
      if (m > 0.5 || m < 0.5 - line * (1.0 + vHi) / (size * (1.0 + 0.4 * vHi))) discard;
      gl_FragColor = vec4(mix(color, ink, vHi), 1.0); }`,
  transparent: true,
});
const lineMaterial = new THREE.LineBasicMaterial({ color: HEX.rule, transparent: true, depthTest: false });

let R = null;             // the region on screen
let loading = 0;          // bumps on every region change, so a slow load cannot land in the wrong region
let pinned = null, hovered = null;  // {kind: 'cluster' | 'group', i}
let needs = true, lastMove = 0, flight = null;

function points(geometry, material, order) {
  const p = new THREE.Points(geometry, material);
  p.frustumCulled = false;  // field positions are encoded, so three.js cannot bound them
  p.renderOrder = order;
  return p;
}

async function show(name) {
  const token = ++loading;  // also cancels a load in flight for another region
  for (const b of document.querySelectorAll('[data-region]')) b.setAttribute('aria-pressed', String(b.dataset.region === name));
  history.replaceState(null, '', `#${name}`);
  if (R?.name === name) return fieldStatus();
  $('#field-status').textContent = 'Loading…';
  const get = async (file) => {
    const r = await fetch(`data/${name}/${file}`);
    if (!r.ok) throw new Error(`${file}: HTTP ${r.status}`);
    return r;
  };
  const meta = await (await get('meta.json')).json();
  const [members, cmd] = await Promise.all([meta.members.file, meta.cmd.file].map(async (f) => (await get(f)).arrayBuffer()));
  if (token !== loading) return;
  clear();
  R = build(name, meta, members, cmd, get);
  showCard(null);
  for (let k = 0; k < meta.field.levels.length && token === loading; k++) {
    if (k > 0 && R.loaded + meta.field.levels[k].count > R.budget) break;
    await addLevel(R, k);
  }
}

async function addLevel(region, k) {
  const level = region.meta.field.levels[k];
  const buf = await (await region.get(level.file)).arrayBuffer();
  if (region !== R || region.levels.length !== k) return;
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(new Uint16Array(buf), 3));
  const p = points(g, fieldMaterial, 0);
  region.field.add(p);
  region.levels.push(p);
  region.loaded += level.count;
  fieldStatus();
  needs = true;
}

function build(name, meta, membersBuf, cmdBuf, get) {
  const group = new THREE.Group(), field = new THREE.Group();
  const region = { name, meta, get, group, field, levels: [], loaded: 0, budget: BUDGET, labels: [] };
  const lo = meta.quant.lo, step = meta.quant.hi.map((h, i) => (h - lo[i]) / 65535);
  quant.lo.value.fromArray(lo);
  quant.step.value.fromArray(step);

  // Members, in parsecs: for drawing, picking and flying to a cluster
  const M = meta.members.count;
  const pos = new Float32Array(membersBuf, 0, 3 * M), cluster = new Uint16Array(membersBuf, 12 * M, M);
  const mg = new THREE.BufferGeometry();
  mg.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  mg.setAttribute('cluster', new THREE.BufferAttribute(cluster, 1));
  group.add(points(mg, memberMaterial, 6));
  region.members = { pos, cluster };

  // Recovered clusters: spheres as wide as the cluster in Hunt & Reffert (its tidal radius), but never
  // smaller than a few pixels. drawSpheres() sizes them every frame
  region.found = meta.clusters.flatMap((c, i) => c.members ? [i] : []);
  region.spheres = new THREE.InstancedMesh(sphereGeometry, sphereMaterial, region.found.length);
  region.spheres.frustumCulled = false;
  group.add(region.spheres);

  // Missed clusters: rings at the catalogue position, sized the same way
  const missed = meta.clusters.map((c, i) => [c, i]).filter(([c]) => !c.recovered);
  const rg = new THREE.BufferGeometry();
  rg.setAttribute('position', new THREE.Float32BufferAttribute(missed.flatMap(([c]) => c.xyz), 3));
  rg.setAttribute('radius', new THREE.Float32BufferAttribute(missed.map(([c]) => c.radius), 1));
  rg.setAttribute('idx', new THREE.Float32BufferAttribute(missed.map(([, i]) => i), 1));
  group.add(points(rg, ringMaterial, 4));
  region.missed = missed.map(([, i]) => i);

  const ug = new THREE.BufferGeometry();
  ug.setAttribute('position', new THREE.Float32BufferAttribute(meta.unverified.xyz, 3));
  ug.setAttribute('idx', new THREE.Float32BufferAttribute(meta.unverified.n.map((_, i) => i), 1));
  region.unverified = points(ug, squareMaterial, 2);
  region.unverified.visible = $('#unverified').checked;
  group.add(region.unverified);

  region.cmd = { bp: new Int16Array(cmdBuf, 0, meta.cmd.count), g: new Uint16Array(cmdBuf, 2 * meta.cmd.count, meta.cmd.count),
                 cat: new Uint8Array(cmdBuf, 4 * meta.cmd.count, meta.cmd.count) };
  group.add(guides(region));
  for (const name of meta.named) {
    const i = meta.clusters.findIndex((c) => c.name === name);
    if (i >= 0) label(region, meta.clusters[i].label, centre(region, i), 'named');
  }
  scene.add(group);
  fieldScene.add(field);
  intro(meta);
  fillFind(meta);
  $('#unverified-label').textContent = `Show the ${num(meta.unverified.count)} unverified groups`;
  frame(region);
  return region;
}

const sphereGeometry = new THREE.SphereGeometry(1, 32, 16);
scene.add(new THREE.HemisphereLight(0xffffff, 0x8a8a87, 2.2));
const sun = new THREE.DirectionalLight(0xffffff, 1.4);
sun.position.set(-1, 1, 3);
scene.add(sun);

const inView = new THREE.Vector3();
const depth = (xyz) => -inView.fromArray(xyz).applyMatrix4(camera.matrixWorldInverse).z;
// A cluster's radius on screen, in CSS pixels
const screenRadius = (c) => Math.max(MIN_RADIUS, c.radius * proj.value / dpr / Math.max(depth(c.xyz), 1e-3));

const scaled = new THREE.Matrix4();
function drawSpheres() {
  const on = hovered ?? pinned;
  R.found.forEach((i, k) => {
    const c = R.meta.clusters[i], d = depth(c.xyz);  // behind the camera (d < 0), at its own size
    const r = on?.kind === 'cluster' && on.i === i ? 0 : Math.max(c.radius, MIN_RADIUS * dpr * d / proj.value);
    R.spheres.setMatrixAt(k, scaled.makeScale(r, r, r).setPosition(...c.xyz));
  });
  R.spheres.instanceMatrix.needsUpdate = true;
}

// The kiloparsec past the field's 99th percentile and the farthest cluster: some clusters lie beyond the field
const farEdge = (meta) => Math.ceil(Math.max(meta.reach, ...meta.clusters.map((c) => Math.hypot(...c.xyz))) / 1000) * 1000;

// The Galactic plane under the wedge: a ray every 2° of longitude from the Sun, an arc every half kiloparsec
function guides(region) {
  const { box: { l: [l0, l1] } } = region.meta;
  const far = farEdge(region.meta), steps = 40, pts = [];
  const arc = (b, d) => {
    for (let s = 0; s < steps; s++) pts.push(toward(l0 + (l1 - l0) * s / steps, b, d), toward(l0 + (l1 - l0) * (s + 1) / steps, b, d));
  };
  for (let l = l0; l <= l1 + 1e-9; l += 2) pts.push(new THREE.Vector3(), toward(l, 0, far));
  for (let d = 500; d <= far; d += 500) arc(0, d);
  for (let d = 1000; d <= far; d += 1000) {
    label(region, `${d / 1000} kpc`, toward(l1, 0, d), 'scale', d);  // frame() moves it to the edge in front
  }
  label(region, 'Sun', new THREE.Vector3(), 'scale sun');
  const lines = new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(pts), lineMaterial);
  lines.renderOrder = 1;
  return lines;
}

function label(region, text, at, kind, distance) {
  const el = document.createElement('span');
  el.className = kind;
  el.textContent = text;
  labelsEl.append(el);
  region.labels.push({ el, at, distance, rank: kind.includes('sun') ? 0 : kind === 'named' ? 1 : 2 });
  region.labels.sort((a, b) => a.rank - b.rank);  // placeLabels() hides a label that would cover one placed before it
}

const projected = new THREE.Vector3();
function placeLabels() {
  const w = box.clientWidth, h = box.clientHeight;
  const taken = [$('#hint'), $('.zoom')].map((el) => [el.offsetLeft, el.offsetTop, el.offsetLeft + el.offsetWidth, el.offsetTop + el.offsetHeight]);
  for (const lab of R.labels) {
    projected.copy(lab.at).project(camera);
    let off = projected.z > 1 || Math.abs(projected.x) > 1.05 || Math.abs(projected.y) > 1.05;
    if (!off) {
      const x = (projected.x + 1) / 2 * w + 7, y = (1 - projected.y) / 2 * h - 15;
      lab.w ||= lab.el.offsetWidth;
      lab.h ||= lab.el.offsetHeight;
      const r = [x - 3, y - 3, x + lab.w + 3, y + lab.h + 3];
      off = taken.some((t) => r[0] < t[2] && r[2] > t[0] && r[1] < t[3] && r[3] > t[1]);
      if (!off) {
        taken.push(r);
        lab.el.style.transform = `translate(${x.toFixed(1)}px, ${y.toFixed(1)}px)`;
      }
    }
    lab.el.style.visibility = off ? 'hidden' : '';
  }
}

function clear() {
  if (!R) return;
  scene.remove(R.group);
  fieldScene.remove(R.field);
  for (const g of [R.group, R.field]) g.traverse((o) => o.geometry?.dispose());
  labelsEl.replaceChildren();
  pinned = hovered = null;
  flight = null;
  R = null;
}

// The whole wedge in view from above, tilted a little so it reads as a solid, with the Sun at the left;
// on a tall screen, with the Sun at the bottom. The distance labels go on the edge nearest the viewer.
// On a wide screen the wedge is much longer than the view is tall, so the view shrinks to the wedge.
function frame(region, height) {
  const { box: { l: [l0, l1], b: [b0, b1] } } = region.meta;
  box.style.height = height ? `${height}px` : '';
  camera.aspect = box.clientWidth / box.clientHeight;
  camera.updateProjectionMatrix();
  const lc = (l0 + l1) / 2, far = farEdge(region.meta), tall = camera.aspect < 0.9;
  const axis = toward(lc, 0, 1), side = toward(lc + 90, 0, 1), north = new THREE.Vector3(0, 0, 1);
  const corners = [[l0, b0], [l0, b1], [l1, b0], [l1, b1]], edges = [[l0, 0], [l1, 0]];
  const hull = [new THREE.Vector3(), ...corners.map(([l, b]) => toward(l, b, far))];
  const target = axis.clone().multiplyScalar(0.45 * far), right = new THREE.Vector3(), up = new THREE.Vector3();
  const place = (dir, distance) => {
    camera.position.copy(target).addScaledVector(dir, distance);
    camera.lookAt(target);
    camera.updateMatrixWorld();
    right.setFromMatrixColumn(camera.matrixWorld, 0);
    up.setFromMatrixColumn(camera.matrixWorld, 1);
  };
  let dir;
  for (const sign of [1, -1]) {
    dir = (tall ? axis.clone().multiplyScalar(-0.3).addScaledVector(side, 0.15 * sign).addScaledVector(north, 0.94)
                : side.clone().multiplyScalar(0.12 * sign).addScaledVector(north, 0.99)).normalize();
    place(dir, far);
    if ((tall ? up : right).dot(axis) > 0) break;
  }
  // Fit: the distance at which the wedge fills 94 % of the view, then centre it, twice over
  const extent = () => {
    const ndc = hull.map((v) => v.clone().project(camera));
    const xs = ndc.map((v) => v.x), ys = ndc.map((v) => v.y);
    return { x: [Math.min(...xs), Math.max(...xs)], y: [Math.min(...ys), Math.max(...ys)] };
  };
  let distance = far;
  for (let pass = 0; pass < 2; pass++) {
    let lo = 0.1 * far, hi = 10 * far;
    for (let i = 0; i < 30; i++) {
      distance = (lo + hi) / 2;
      place(dir, distance);
      const e = extent();
      if (Math.max(-e.x[0], e.x[1], -e.y[0], e.y[1]) > 0.94) lo = distance; else hi = distance;
    }
    const e = extent(), half = distance * Math.tan(camera.fov * RAD / 2);
    target.addScaledVector(right, (e.x[0] + e.x[1]) / 2 * half * camera.aspect).addScaledVector(up, (e.y[0] + e.y[1]) / 2 * half);
  }
  if (!height && !narrow.matches) {  // the wedge fills 75 % of the height, leaving room for its labels
    const e = extent(), fit = Math.max(300, Math.ceil(box.clientHeight * (e.y[1] - e.y[0]) / 2 / 0.75));
    if (fit < box.clientHeight) return frame(region, fit);
  }
  place(dir, distance);
  controls.target.copy(target);
  controls.update();
  region.home = { dir, distance, target: target.clone() };
  const screen = ([l, b]) => toward(l, b, far).project(camera);
  const [l, b] = edges.sort((p, q) => tall ? screen(q).x - screen(p).x : screen(p).y - screen(q).y)[0];
  for (const lab of region.labels) if (lab.distance) lab.at.copy(toward(l, b, lab.distance));
  needs = true;
}

// Where a cluster is drawn, and how close to come to see it whole
function centre(region, i) {
  const m = region.meta.clusters[i].members;
  if (!m) return new THREE.Vector3().fromArray(region.meta.clusters[i].xyz);
  const c = new THREE.Vector3(), p = region.members.pos;
  for (let j = m[0]; j < m[0] + m[1]; j++) c.x += p[3 * j], c.y += p[3 * j + 1], c.z += p[3 * j + 2];
  return c.divideScalar(m[1]);
}

function viewDistance(region, i, at) {
  const m = region.meta.clusters[i].members;
  if (!m) return 40;
  const p = region.members.pos, d = [];
  for (let j = m[0]; j < m[0] + m[1]; j++) d.push(Math.hypot(p[3 * j] - at.x, p[3 * j + 1] - at.y, p[3 * j + 2] - at.z));
  d.sort((a, b) => a - b);
  return Math.max(12, 5.5 * d[Math.floor(0.9 * (d.length - 1))]);
}

function flyTo(at, distance, dir) {
  const from = { target: controls.target.clone(), distance: camera.position.distanceTo(controls.target),
                 dir: camera.position.clone().sub(controls.target).normalize() };
  flight = { t0: performance.now(), from, to: { target: at.clone(), distance, dir: dir ?? from.dir } };
  if (reduceMotion) flight.t0 -= 1e4;
}

function fly() {
  const t = Math.min(1, (performance.now() - flight.t0) / 1100);
  const e = t < 0.5 ? 4 * t ** 3 : 1 - (-2 * t + 2) ** 3 / 2;
  const { from, to } = flight, dir = from.dir.clone().lerp(to.dir, e).normalize();
  controls.target.lerpVectors(from.target, to.target, e);
  const d = Math.exp(Math.log(from.distance) + (Math.log(to.distance) - Math.log(from.distance)) * e);
  camera.position.copy(controls.target).addScaledVector(dir, d);
  lastMove = performance.now();
  if (t === 1) flight = null;
}

// What is under the pointer: a sphere, a ring, or (if shown) an unverified group
const viewProj = new THREE.Matrix4();
function pick(px, py) {
  if (!R) return null;
  const m = viewProj.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse).elements;
  const w = box.clientWidth, h = box.clientHeight;
  const dist2 = (x, y, z) => {  // squared pixels from the pointer; Infinity behind the camera
    const cw = m[3] * x + m[7] * y + m[11] * z + m[15];
    if (cw <= 0) return Infinity;
    const dx = ((m[0] * x + m[4] * y + m[8] * z + m[12]) / cw + 1) / 2 * w - px;
    const dy = (1 - (m[1] * x + m[5] * y + m[9] * z + m[13]) / cw) / 2 * h - py;
    return dx * dx + dy * dy;
  };
  let best = null, bestD = Infinity;
  const consider = (d, reach, hit) => { if (d < reach * reach && d < bestD) { bestD = d; best = hit; } };
  for (const i of [...R.found, ...R.missed]) {
    const c = R.meta.clusters[i];
    consider(dist2(...c.xyz), screenRadius(c) + 4, { kind: 'cluster', i });
  }
  if (R.unverified.visible) {
    const u = R.meta.unverified.xyz;
    for (let i = 0; i < u.length / 3; i++) consider(dist2(u[3 * i], u[3 * i + 1], u[3 * i + 2]), 8, { kind: 'group', i });
  }
  return best;
}

function highlight() {
  const on = hovered ?? pinned;
  const cluster = on?.kind === 'cluster' ? on.i : -1, group = on?.kind === 'group' ? on.i : -1;
  memberMaterial.uniforms.hi.value = ringMaterial.uniforms.hi.value = cluster;
  squareMaterial.uniforms.hi.value = group;
  needs = true;
}

const same = (a, b) => a?.kind === b?.kind && a?.i === b?.i;

function hover(hit) {
  if (same(hit, hovered)) return;
  hovered = hit;
  box.classList.toggle('pointing', !!hit);
  highlight();
  if (hit) showCard(hit);
  else if (pinned) showCard(pinned);
}

function pin(hit, go = false) {
  pinned = hit;
  highlight();
  showCard(hit);
  if (hit && go) goTo(hit);
}

// The first kiloparsec from above, where the clusters crowd; or the whole region again
function zoomTo(view) {
  if (!R) return;
  const { home } = R, { box: { l: [l0, l1] } } = R.meta;
  if (view === 'home') return flyTo(home.target, home.distance, home.dir);
  const across = camera.aspect < 0.9 ? 1 : camera.aspect;  // the wedge runs across a wide view, up a tall one
  flyTo(toward((l0 + l1) / 2, 0, 500), 600 / (Math.tan(camera.fov * RAD / 2) * across), home.dir);
}

function goTo(hit) {
  if (hit.kind === 'cluster') {
    const at = centre(R, hit.i);
    flyTo(at, viewDistance(R, hit.i, at));
  } else {
    flyTo(new THREE.Vector3().fromArray(R.meta.unverified.xyz, 3 * hit.i), 40);
  }
}

function showCard(hit) {
  if (!R) return;
  card.innerHTML = !hit ? emptyCard() : hit.kind === 'cluster' ? clusterCard(hit.i) : groupCard(hit.i);
}

function emptyCard() {
  const m = R.meta;
  return `<div class="info"><p class="kicker">${esc(m.title)}</p>
    <p>Point at a blue sphere or a red ring, or tap it, for its card. Every cluster is also in the list under the map.</p>
    <dl class="facts"><dt>Recovered</dt><dd>${m.recovered} of ${m.recoverable} (${pct(m.recovered / m.recoverable)})</dd>
    <dt>Highly reliable</dt><dd>${m.reliable_recovered} of ${m.reliable} (${pct(m.reliable_recovered / m.reliable)})</dd>
    <dt>Missed</dt><dd>${m.recoverable - m.recovered}</dd><dt>Field stars</dt><dd>${num(m.field.count)}</dd></dl></div>`;
}

function clusterCard(i) {
  const c = R.meta.clusters[i], plain = c.name.replaceAll('_', ' ');
  const status = c.recovered ? '<p class="kicker rec">Recovered</p>' : `<p class="kicker miss">Missed · ${esc(c.how.replace(/ \(.*\)$/, ''))}</p>`;
  const others = c.runs.filter((k) => k !== c.run);
  const group = !c.size ? `No group of the leaf runs holds any of its members.`
    : `${c.recovered ? 'Shown' : 'Closest'}: a group of ${num(c.size)} stars in the run with min_cluster_size ${c.run}.`
      + (others.length ? ` Also recovered at ${list(others)}.` : '');
  const figure = cmdFigure(c);
  return `<div class="info"><div class="head"><div>${status}<h3>${esc(c.label)}</h3>${c.label !== plain ? `<span class="alias">${esc(plain)} in Hunt &amp; Reffert</span>` : ''}</div>
    ${same(pinned, { kind: 'cluster', i }) ? '<button type="button" class="close" data-act="unpin" aria-label="Close">×</button>' : ''}</div>
    <p class="meta">${kpc(c.dist)} · ${num(c.n)} members in our sample</p>
    <div class="scores"><div><span>Recovery</span><b>${pct(c.recovery)}</b></div><div><span>Purity</span><b>${pct(c.purity)}</b></div></div>
    <p>${group} ${c.reliable ? 'Highly reliable in Hunt &amp; Reffert.' : 'Not rated highly reliable by Hunt &amp; Reffert.'}</p>
    <button type="button" class="go" data-act="go" data-i="${i}">Go to ${esc(c.label)}</button></div>
    ${figure ? `<div class="plot">${figure}</div>` : ''}`;
}

const list = (xs) => xs.length < 2 ? `${xs[0]}` : `${xs.slice(0, -1).join(', ')} and ${xs.at(-1)}`;

function groupCard(i) {
  const u = R.meta.unverified, at = new THREE.Vector3().fromArray(u.xyz, 3 * i);
  return `<div class="info"><div class="head"><div><p class="kicker">Unverified group</p><h3>${num(u.n[i])} stars</h3></div></div>
    <p class="meta">${kpc(at.length())}, the median of its stars · min_cluster_size ${R.meta.unverified_run}</p>
    <p>No Hunt &amp; Reffert cluster holds most of its stars. Hunt &amp; Reffert discard groups like this with a significance
    test and a colour-magnitude classifier, which this project does not have, so it is not counted as a cluster.</p></div>`;
}

// Colour-magnitude diagram: G against BP − RP, blue in both, red ours only, grey theirs only
function cmdFigure(c) {
  const [start, count] = c.cmd;
  if (!count) return '';
  const { bp, g, cat } = R.cmd, W = 300, H = 250, L = 34, B = 38, T = 6, Rm = 6;
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (let j = start; j < start + count; j++) {
    x0 = Math.min(x0, bp[j]); x1 = Math.max(x1, bp[j]); y0 = Math.min(y0, g[j]); y1 = Math.max(y1, g[j]);
  }
  [x0, x1, y0, y1] = [x0 / 1000 - 0.1, x1 / 1000 + 0.1, y0 / 1000 - 0.3, y1 / 1000 + 0.3];
  const sx = (x) => L + (x - x0) / (x1 - x0) * (W - L - Rm), sy = (y) => T + (y - y0) / (y1 - y0) * (H - T - B);
  const dots = ['', '', ''], r = count > 1500 ? 1.1 : count > 300 ? 1.5 : 2;
  for (let j = start; j < start + count; j++) {
    const x = sx(bp[j] / 1000), y = sy(g[j] / 1000);
    dots[cat[j]] += `M${(x - r).toFixed(1)} ${y.toFixed(1)}a${r} ${r} 0 1 0 ${2 * r} 0a${r} ${r} 0 1 0 ${-2 * r} 0`;
  }
  const tx = ticks(x0, x1, 4).map((v) => `<line class="axis" x1="${sx(v)}" x2="${sx(v)}" y1="${H - B}" y2="${H - B + 4}"/>
    <text x="${sx(v)}" y="${H - B + 15}" text-anchor="middle">${+v.toFixed(2)}</text>`).join('');
  const ty = ticks(y0, y1, 5).map((v) => `<line class="axis" x1="${L - 4}" x2="${L}" y1="${sy(v)}" y2="${sy(v)}"/>
    <text x="${L - 7}" y="${sy(v) + 3.5}" text-anchor="end">${+v.toFixed(1)}</text>`).join('');
  const n = [0, 1, 2].map((k) => { let s = 0; for (let j = start; j < start + count; j++) s += cat[j] === k; return s; });
  return `<svg class="cmd" viewBox="0 0 ${W} ${H}" role="img" aria-label="Colour-magnitude diagram of ${esc(c.label)}: ${n[0]} stars in both, ${n[1]} in our group only, ${n[2]} Hunt and Reffert members it misses">
      <path class="theirs" d="${dots[2]}"/><path class="both" d="${dots[0]}"/><path class="ours" d="${dots[1]}"/>
      <line class="axis" x1="${L}" x2="${L}" y1="${T}" y2="${H - B}"/><line class="axis" x1="${L}" x2="${W - Rm}" y1="${H - B}" y2="${H - B}"/>
      ${tx}${ty}<text x="${(L + W - Rm) / 2}" y="${H - 4}" text-anchor="middle">BP − RP</text><text x="${L + 6}" y="${T + 9}">G</text>
    </svg>
    <div class="cmd-key"><span><i style="background:var(--recovered)"></i>in both</span><span><i style="background:var(--missed)"></i>ours only</span><span><i style="background:var(--theirs)"></i>theirs only</span></div>`;
}

function ticks(a, b, n) {
  const raw = (b - a) / n, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((s) => s * mag).find((s) => s >= raw);
  const out = [];
  for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) out.push(v);
  return out;
}

function intro(meta) {
  $('#title').textContent = meta.title;
  const [l0, l1] = meta.box.l, [b0, b1] = meta.box.b;
  const share = `<strong>${meta.recovered} of the ${meta.recoverable}</strong> Hunt &amp; Reffert clusters with at least 10 members in our sample come back (${pct(meta.recovered / meta.recoverable)})`;
  const standing = meta.protocol
    ? `, with leaf selection. That is post hoc: the protocol used ${meta.protocol.selection} and recovered ${meta.protocol.recovered}.`
    : ', with leaf selection fixed in advance.';
  const deg = (v) => `${v < 0 ? '\u2212' : v > 0 ? '+' : ''}${Math.abs(v)}°`;
  $('#lede').innerHTML = `${share}${standing} ${num(meta.stars)} Gaia DR3 stars, Galactic longitude ${l0}–${l1}°, latitude ${deg(b0)} to ${deg(b1)}.`;
  $('#depth-sigma').textContent = `${meta.depth.sigma} pc`;
  $('#depth-r50').textContent = `${meta.depth.r50} pc`;
  $('#depth-region').textContent = `In ${meta.title}`;
  document.title = `${meta.title} · gaia-clusters`;
}

function fillFind(meta) {
  const order = (a, b) => meta.clusters[a].label.localeCompare(meta.clusters[b].label, 'en', { numeric: true });
  const all = meta.clusters.map((_, i) => i);
  const option = (i) => `<option value="${i}">${esc(meta.clusters[i].label)}</option>`;
  const rec = all.filter((i) => meta.clusters[i].recovered).sort(order), miss = all.filter((i) => !meta.clusters[i].recovered).sort(order);
  $('#find').innerHTML = `<option value="">${all.length} clusters</option>
    <optgroup label="Recovered (${rec.length})">${rec.map(option).join('')}</optgroup>
    <optgroup label="Missed (${miss.length})">${miss.map(option).join('')}</optgroup>`;
}

function fieldStatus() {
  const { loaded, meta, levels, budget } = R, next = meta.field.levels[levels.length];
  $('#field-status').textContent = `Field: ${num(loaded)} of ${num(meta.field.count)} stars`;
  $('#more').hidden = !next || loaded + next.count <= budget;  // shown only where the budget stopped the loading
}

// Events
if (coarse) $('#hint').textContent = 'Drag to turn · pinch to zoom · two fingers to pan · tap a cluster';
const failed = (err) => { $('#field-status').textContent = `Could not load the region: ${err.message}`; };
for (const b of document.querySelectorAll('[data-region]')) b.addEventListener('click', () => show(b.dataset.region).catch(failed));
$('#find').addEventListener('change', (e) => { if (e.target.value !== '') pin({ kind: 'cluster', i: +e.target.value }, true); });
$('#unverified').addEventListener('change', (e) => {
  if (!R) return;
  R.unverified.visible = e.target.checked;
  document.querySelector('.unverified-key').hidden = !e.target.checked;
  if (!e.target.checked && (pinned?.kind === 'group' || hovered?.kind === 'group')) { hovered = null; pin(null); }
  needs = true;
});
$('#more').addEventListener('click', async () => {
  const region = R;
  region.budget = Infinity;
  $('#more').hidden = true;
  for (let k = region.levels.length; k < region.meta.field.levels.length && region === R; k++) await addLevel(region, k);
});
for (const b of document.querySelectorAll('[data-zoom]')) b.addEventListener('click', () => zoomTo(b.dataset.zoom));
card.addEventListener('click', (e) => {
  const act = e.target.closest('[data-act]')?.dataset.act;
  if (act === 'unpin') { pin(null); $('#find').value = ''; }
  if (act === 'go') goTo({ kind: 'cluster', i: +e.target.closest('[data-act]').dataset.i });
});
addEventListener('keydown', (e) => { if (e.key === 'Escape' && pinned) { pin(null); $('#find').value = ''; } });

let down = null, queued = null;
canvas.addEventListener('pointerdown', (e) => { down = { x: e.clientX, y: e.clientY, t: performance.now() }; flight = null; });
canvas.addEventListener('pointerup', (e) => {
  if (!down || Math.hypot(e.clientX - down.x, e.clientY - down.y) > 5 || performance.now() - down.t > 600) return;
  const r = canvas.getBoundingClientRect(), hit = pick(e.clientX - r.left, e.clientY - r.top);
  pin(hit);
  $('#find').value = hit?.kind === 'cluster' ? hit.i : '';
});
canvas.addEventListener('pointermove', (e) => {
  if (e.buttons || e.pointerType === 'touch') return;
  const r = canvas.getBoundingClientRect();
  if (!queued) requestAnimationFrame(() => { hover(pick(queued.x, queued.y)); queued = null; });
  queued = { x: e.clientX - r.left, y: e.clientY - r.top };
});
canvas.addEventListener('pointerleave', () => hover(null));
controls.addEventListener('change', () => { lastMove = performance.now(); needs = true; });

// Frame again when the width changes; not on height alone, which a phone's address bar changes while scrolling
let framedWidth = innerWidth, reframe;
addEventListener('resize', () => {
  if (innerWidth === framedWidth) return;
  framedWidth = innerWidth;
  clearTimeout(reframe);
  reframe = setTimeout(() => R && frame(R), 150);
});

new ResizeObserver(() => {
  const w = box.clientWidth, h = box.clientHeight;
  renderer.setSize(w, h, false);
  const size = renderer.getDrawingBufferSize(tone.material.uniforms.size.value);
  proj.value = size.y / (2 * Math.tan(camera.fov * RAD / 2));
  density.setSize(size.x, size.y);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  needs = true;
}).observe(box);

// Fewer field stars while the view moves: any prefix of the field is a uniform sample
function draw() {
  requestAnimationFrame(draw);
  if (flight) fly();
  controls.update();
  const moving = performance.now() - lastMove < 250;
  if (!R || (!needs && !moving && performance.now() - lastMove > 400)) return;
  needs = false;
  let drawn = 0;
  for (const p of R.levels) {
    p.visible = !moving || drawn < MOVING;
    if (p.visible) drawn += p.geometry.attributes.position.count;
  }
  tone.material.uniforms.gain.value = R.meta.field.count / Math.max(drawn, 1);  // each star drawn stands for this many
  drawSpheres();
  renderer.setRenderTarget(density);
  renderer.setClearColor(0x000000, 0);
  renderer.render(fieldScene, camera);
  renderer.setRenderTarget(null);
  renderer.setClearColor(HEX.paper, 1);
  renderer.render(scene, camera);
  placeLabels();
}

const start = location.hash.slice(1);
show(['carina', 'auriga'].includes(start) ? start : 'carina').catch(failed);
draw();
