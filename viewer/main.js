// The gaia-clusters viewer: one scored region in three dimensions.
// The files come from scripts/export_viewer.py; their layout is in src/gaia_clusters/viewer.py.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TOURS } from './tour.js';

THREE.ColorManagement.enabled = false;  // hex in, hex out: the colours of the README figures

const HEX = { paper: 0xf5f5f2, ink: 0x222222, recovered: 0x1f3b73, missed: 0xc8553d, field: 0x9a9a9a,
              rule: 0xc4c5bf, unverified: 0x6b6b67 };
const coarse = matchMedia('(pointer: coarse)').matches;
const BUDGET = coarse ? 600_000 : Infinity;  // field stars loaded before asking
const MOVING = coarse ? 250_000 : 700_000;   // field stars drawn while the view moves
const DARKEST = 400;                         // field stars per pixel drawn in the darkest grey
const MIN_RADIUS = 4;                        // CSS pixels: no cluster is drawn smaller than this
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;
const MAS_MYR = 4.8481e-3;                     // radians a star crosses in a million years at 1 mas/yr
const SPAN = 5;                                // Myr either side of now that the time slider reaches
// Colour by age: Hunt & Reffert's median log10(age / yr) in four bands, an ordinal one-hue ramp, lighter = younger
const AGES = [{ below: 7, label: 'under 10 Myr', hex: '#86b6ef' }, { below: 8, label: '10–100 Myr', hex: '#3987e5' },
              { below: 9, label: '0.1–1 Gyr', hex: '#1c5cab' }, { below: Infinity, label: 'over 1 Gyr', hex: '#0d366b' }];
const ageBand = (c) => AGES.find((a) => c.age[1] < a.below);
const RAD = Math.PI / 180;
const num = new Intl.NumberFormat('en-GB').format;
const pct = (x) => `${Math.round(100 * x)} %`;
const kpc = (pc) => `${(pc / 1000).toFixed(pc < 1000 ? 2 : 1)} kpc`;
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const $ = (s) => document.querySelector(s);
$('#age-key').innerHTML = `age: ${AGES.map((a) => `<span><i style="background:${a.hex}"></i>${a.label}</span>`).join('')}`;
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
// Moving in time: every star drifts along the sky at its proper motion (pm, in steps of pmStep mas/yr). Unless
// seen from the Sun (relative = 0), the field's typical motion at the star's distance is taken off: the Sun's own
// motion and the Galaxy's rotation shift every star at a distance alike (meta.reference, in bins of refStep pc).
// Most stars have no radial velocity, so nothing moves in depth, and the paths are straight lines.
const REF_BINS = 64;
const motion = { myr: { value: 0 }, pmStep: { value: 0.01 }, relative: { value: 1 }, refStep: { value: 250 }, refCount: { value: 1 },
                 ref: { value: Array.from({ length: REF_BINS }, () => new THREE.Vector2()) } };
const drift = `uniform float myr; uniform float pmStep; uniform float relative; uniform float refStep; uniform float refCount;
  uniform vec2 ref[${REF_BINS}]; attribute vec2 pm;
  vec2 typical(float d) { float x = clamp(d / refStep - 0.5, 0.0, refCount - 1.0);
    return mix(ref[int(x)], ref[int(min(x + 1.0, refCount - 1.0))], fract(x)); }
  vec3 drift(vec3 p) { float d = length(p); if (myr == 0.0 || d < 1e-3) return p;
    vec3 r = p / d, el = normalize(vec3(-r.y, r.x, 0.0)), eb = cross(r, el);
    vec2 m = pm * pmStep - relative * typical(d);
    return p + d * ${MAS_MYR} * myr * (m.x * el + m.y * eb); }`;
const fieldScene = new THREE.Scene();
// Half floats count far past the 255 of a byte; a browser without them gets the byte, which clips
const halfFloat = ['EXT_color_buffer_float', 'EXT_color_buffer_half_float'].some((e) => renderer.extensions.has(e));
const unit = halfFloat ? 1 : 1 / 255;
const density = new THREE.WebGLRenderTarget(1, 1, { depthBuffer: false, type: halfFloat ? THREE.HalfFloatType : THREE.UnsignedByteType,
                                                    minFilter: THREE.NearestFilter, magFilter: THREE.NearestFilter });
const fieldMaterial = new THREE.ShaderMaterial({
  uniforms: { ...quant, ...motion, size: { value: dpr } },  // one CSS pixel, so the counts do not depend on the screen
  vertexShader: `uniform vec3 lo; uniform vec3 step; uniform float size; ${drift}
    void main() { gl_Position = projectionMatrix * modelViewMatrix * vec4(drift(lo + position * step), 1.0); gl_PointSize = size; }`,
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
  uniforms: { ...motion, size: { value: 3.2 * dpr }, ink: { value: new THREE.Color(HEX.ink) }, hi: { value: -1 } },
  vertexShader: `uniform float size; uniform float hi; attribute float cluster; ${drift}
    void main() { gl_Position = abs(cluster - hi) < 0.5 ? projectionMatrix * modelViewMatrix * vec4(drift(position), 1.0) : vec4(2.0, 2.0, 2.0, 1.0);
      gl_PointSize = size; }`,
  fragmentShader: `uniform vec3 ink;
    void main() { vec2 c = gl_PointCoord - 0.5; if (dot(c, c) > 0.25) discard; gl_FragColor = vec4(ink, 1.0); }`,
  transparent: true, depthTest: false,
});
const sphereMaterial = new THREE.MeshLambertMaterial({ color: 0xffffff });  // each sphere's colour is its instance colour (paint)
const ringMaterial = new THREE.ShaderMaterial({
  uniforms: { ...motion, ink: { value: new THREE.Color(HEX.ink) }, hi: { value: -1 },
              line: { value: 1.2 * dpr }, scale: { value: dpr }, proj, least: { value: 2 * MIN_RADIUS * dpr } },
  vertexShader: `uniform float hi; uniform float scale; uniform float proj; uniform float least; attribute float radius; attribute float idx;
    attribute vec3 tint; varying vec3 vTint; varying float vHi; varying float vSize; ${drift}
    void main() { vHi = abs(idx - hi) < 0.5 ? 1.0 : 0.0; vTint = tint; vec4 p = modelViewMatrix * vec4(drift(position), 1.0);
      vSize = max(least, 2.0 * radius * proj / -p.z) + 2.0 * scale * vHi;
      gl_Position = projectionMatrix * p; gl_PointSize = vSize; }`,
  fragmentShader: `uniform vec3 ink; uniform float line; varying vec3 vTint; varying float vHi; varying float vSize;
    void main() { float r = length(gl_PointCoord - 0.5) * vSize, w = line * (1.0 + vHi);
      if (r > 0.5 * vSize || r < 0.5 * vSize - w) discard; gl_FragColor = vec4(mix(vTint, ink, vHi), 1.0); }`,
  transparent: true,
});
const squareMaterial = new THREE.ShaderMaterial({
  uniforms: { ...motion, color: { value: new THREE.Color(HEX.unverified) }, ink: { value: new THREE.Color(HEX.ink) }, hi: { value: -1 },
              size: { value: 7 * dpr }, line: { value: 1.3 * dpr } },
  vertexShader: `uniform float size; uniform float hi; attribute float idx; varying float vHi; ${drift}
    void main() { vHi = abs(idx - hi) < 0.5 ? 1.0 : 0.0;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(drift(position), 1.0); gl_PointSize = size * (1.0 + 0.4 * vHi); }`,
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
  const [members, cmd, cmdPm, nearPm] = await Promise.all([meta.members.file, meta.cmd.file, meta.cmd.pm, meta.near.file]
    .map(async (f) => (await get(f)).arrayBuffer()));
  if (token !== loading) return;
  clear();
  R = build(name, meta, get);
  R.near = new Int16Array(nearPm);
  R.score = makeView(R, null, meta, members, cmd, cmdPm);
  setView(R, R.score);
  runStrip(R);
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
  if (region.motions) addMotion(region, k);
  region.loaded += level.count;
  fieldStatus();
  needs = true;
}

const pmAttribute = async (region, file) => new THREE.BufferAttribute(new Int16Array(await (await region.get(file)).arrayBuffer()), 2);
async function addMotion(region, k) {
  region.levels[k].geometry.setAttribute('pm', await pmAttribute(region, region.meta.field.levels[k].pm));
  needs = true;
}
// The members' and the field's proper motions, the first time the clock leaves now. Levels that arrive later
// bring their own (addLevel), and so does a view shown later (setView)
function loadMotions(region) {
  region.motions ??= Promise.all([viewMotion(region, region.view), ...region.levels.map((_, k) => addMotion(region, k))]);
  return region.motions;
}
function viewMotion(region, v) {
  v.motion ??= pmAttribute(region, v.data.members.pm).then((a) => { v.members.geometry.setAttribute('pm', a); needs = true; });
  return v.motion;
}

function build(name, meta, get) {
  const group = new THREE.Group(), field = new THREE.Group();
  const region = { name, meta, get, group, field, levels: [], loaded: 0, budget: BUDGET, labels: [], views: new Map() };
  const lo = meta.quant.lo, step = meta.quant.hi.map((h, i) => (h - lo[i]) / 65535);
  quant.lo.value.fromArray(lo);
  quant.step.value.fromArray(step);

  motion.pmStep.value = meta.pm_step;
  const ref = meta.reference.pm.slice(0, REF_BINS);
  ref.forEach(([l, b], k) => motion.ref.value[k].set(l, b));
  motion.refCount.value = ref.length;
  motion.refStep.value = meta.reference.step;
  region.typical = (d) => {
    const x = Math.max(0, Math.min(d / meta.reference.step - 0.5, ref.length - 1)), i = Math.floor(x), j = Math.min(i + 1, ref.length - 1);
    return ref[i].map((v, c) => v + (ref[j][c] - v) * (x - i));
  };

  // Recovered clusters: spheres as wide as the cluster in Hunt & Reffert (its tidal radius), but never
  // smaller than a few pixels. drawSpheres() sizes them every frame; a view uses the first view.found.length
  region.spheres = new THREE.InstancedMesh(sphereGeometry, sphereMaterial, meta.clusters.length);
  region.spheres.frustumCulled = false;
  group.add(region.spheres);
  group.add(guides(region));
  for (const name of meta.named) {  // at the catalogue position: the label must not jump between views
    const c = meta.clusters.find((c) => c.name === name);
    if (c) label(region, c.label, new THREE.Vector3().fromArray(c.xyz), 'named').pm = c.pm;
  }
  scene.add(group);
  fieldScene.add(field);
  intro(meta);
  galaxy(name);
  frame(region);
  return region;
}

// A view of the region: the score's (run null: each cluster from the run that recovers it best, as in meta.json)
// or one run's (run-<k>.json). It brings the recovered clusters' members, in parsecs, the missed clusters' rings,
// the unverified groups and the diagrams' stars; the clusters' fixed facts (place, motion, age) stay meta's
function makeView(region, run, data, membersBuf, cmdBuf, cmdPmBuf) {
  const { meta } = region;
  const clusters = run === null ? meta.clusters : meta.clusters.map((c, i) => ({ ...c, ...data.clusters[i], run, runs: [] }));
  const M = data.members.count, C = data.cmd.count;
  const pos = new Float32Array(membersBuf, 0, 3 * M), cluster = new Uint16Array(membersBuf, 12 * M, M);
  const mg = new THREE.BufferGeometry();
  mg.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  mg.setAttribute('cluster', new THREE.BufferAttribute(cluster, 1));

  // Missed clusters: rings at the catalogue position, sized as the spheres
  const missed = clusters.flatMap((c, i) => c.recovered ? [] : [i]);
  const rg = new THREE.BufferGeometry();
  rg.setAttribute('position', new THREE.Float32BufferAttribute(missed.flatMap((i) => clusters[i].xyz), 3));
  rg.setAttribute('radius', new THREE.Float32BufferAttribute(missed.map((i) => clusters[i].radius), 1));
  rg.setAttribute('idx', new THREE.Float32BufferAttribute(missed, 1));
  rg.setAttribute('pm', new THREE.Float32BufferAttribute(missed.flatMap((i) => clusters[i].pm.map((v) => v / meta.pm_step)), 2));
  rg.setAttribute('tint', new THREE.Float32BufferAttribute(new Float32Array(3 * missed.length), 3));

  const u = data.unverified, ug = new THREE.BufferGeometry();
  ug.setAttribute('position', new THREE.Float32BufferAttribute(u.xyz, 3));
  ug.setAttribute('idx', new THREE.Float32BufferAttribute(u.n.map((_, i) => i), 1));
  ug.setAttribute('pm', new THREE.Float32BufferAttribute(u.pm.map((v) => v / meta.pm_step), 2));

  return { run, data, clusters, missed, found: clusters.flatMap((c, i) => c.members ? [i] : []),
           unverifiedRun: run ?? meta.unverified_run,
           members: { pos, cluster, geometry: mg, points: points(mg, memberMaterial, 6) },
           rings: points(rg, ringMaterial, 4), unverified: points(ug, squareMaterial, 2),
           cmd: { bp: new Int16Array(cmdBuf, 0, C), g: new Uint16Array(cmdBuf, 2 * C, C), cat: new Uint8Array(cmdBuf, 4 * C, C),
                  pm: new Int16Array(cmdPmBuf) } };
}

function setView(region, v) {
  const old = region.view;
  if (old) region.group.remove(old.members.points, old.rings, old.unverified);
  region.view = v;
  region.group.add(v.members.points, v.rings, v.unverified);
  v.unverified.visible = $('#unverified').checked;
  region.spheres.count = v.found.length;
  if (region.motions) viewMotion(region, v);
  $('#unverified-label').textContent = `Show the ${num(v.data.unverified.count)} unverified groups`;
  if (old && old !== v) {  // each run has its own unverified groups
    if (pinned?.kind === 'group') pinned = null;
    if (hovered?.kind === 'group') hovered = null;
  }
  fillFind(v);
  paint(region);
  highlight();
  showCard(hovered ?? pinned);
}

// A single run's view, fetched the first time it is asked for
function loadView(region, run) {
  if (!region.views.has(run)) {
    region.views.set(run, (async () => {
      const data = await (await region.get(`run-${run}.json`)).json();
      const [m, c, cp] = await Promise.all([data.members.file, data.cmd.file, data.cmd.pm]
        .map(async (f) => (await region.get(f)).arrayBuffer()));
      return makeView(region, run, data, m, c, cp);
    })());
  }
  return region.views.get(run);
}

// The run slider: the score's view, or one run at a time from the sweep (scripts/sweep_region.py)
function runStrip(region) {
  const sweep = region.meta.runs?.sweep ?? [];
  $('#runs').hidden = !sweep.length;
  if (!sweep.length) return;
  $('#run').max = sweep.length - 1;
  chooseRun().catch(failed);
}
async function chooseRun() {
  const region = R;
  if (!region?.meta.runs?.sweep.length) return;
  const { sweep, protocol } = region.meta.runs, on = $('#one-run').checked, s = sweep[+$('#run').value];
  $('#run').disabled = !on;
  const n = region.meta.clusters.length;
  $('#run-out').innerHTML = on
    ? `min_cluster_size <b>${s.run}</b> · ${s.recovered} of ${n} recovered · ${num(s.groups)} groups, ${num(s.unverified)} unverified`
    : `the score: any of runs ${list(protocol)} · ${region.meta.recovered} of ${n} recovered`;
  region.want = on ? s.run : null;
  if (!on) return region.view === region.score || setView(region, region.score);
  const v = await loadView(region, s.run);
  if (region === R && region.want === s.run && region.view !== v) setView(region, v);
}

// Spheres blue and rings red, or both by age
const tint = new THREE.Color();
function paint(region) {
  const byAge = $('#by-age').checked, colour = (c, plain) => tint.set(byAge ? ageBand(c).hex : plain), v = region.view;
  v.found.forEach((i, k) => region.spheres.setColorAt(k, colour(v.clusters[i], HEX.recovered)));
  if (region.spheres.instanceColor) region.spheres.instanceColor.needsUpdate = true;
  const t = v.rings.geometry.attributes.tint;
  v.missed.forEach((i, k) => colour(v.clusters[i], HEX.missed).toArray(t.array, 3 * k));
  t.needsUpdate = true;
  $('#age-key').hidden = !byAge;
  needs = true;
}

const sphereGeometry = new THREE.SphereGeometry(1, 32, 16);
scene.add(new THREE.HemisphereLight(0xffffff, 0x8a8a87, 2.2));
const sun = new THREE.DirectionalLight(0xffffff, 1.4);
sun.position.set(-1, 1, 3);
scene.add(sun);

// Where a point is at the current time, as the shaders' drift() puts it; pm in mas/yr
const along = new THREE.Vector3(), across = new THREE.Vector3(), radial = new THREE.Vector3();
function drifted(p, pm, out = new THREE.Vector3()) {
  const t = motion.myr.value, d = p.length();
  out.copy(p);
  if (!t || !pm || d < 1e-3) return out;
  radial.copy(p).divideScalar(d);
  along.set(-radial.y, radial.x, 0).normalize();
  across.crossVectors(radial, along);
  const k = d * MAS_MYR * t, [l, b] = motion.relative.value ? R.typical(d) : [0, 0];
  return out.addScaledVector(along, k * (pm[0] - l)).addScaledVector(across, k * (pm[1] - b));
}
const clusterAt = (c, out = new THREE.Vector3()) => drifted(out.fromArray(c.xyz), c.pm, out);

const inView = new THREE.Vector3(), at = new THREE.Vector3();
const depth = (v) => -inView.copy(v).applyMatrix4(camera.matrixWorldInverse).z;
// A cluster's radius on screen, in CSS pixels
const screenRadius = (c) => Math.max(MIN_RADIUS, c.radius * proj.value / dpr / Math.max(depth(clusterAt(c, at)), 1e-3));

const scaled = new THREE.Matrix4();
function drawSpheres() {
  const on = hovered ?? pinned;
  R.view.found.forEach((i, k) => {
    const c = R.meta.clusters[i], p = clusterAt(c, at), d = depth(p);  // behind the camera (d < 0), at its own size
    const r = on?.kind === 'cluster' && on.i === i ? 0 : Math.max(c.radius, MIN_RADIUS * dpr * d / proj.value);
    R.spheres.setMatrixAt(k, scaled.makeScale(r, r, r).setPosition(p));
  });
  R.spheres.instanceMatrix.needsUpdate = true;
}

// Recovery by bin, every region side by side: horizontal bars, one colour per region in the buttons' order. Drawn
// W units wide for a figure W pixels wide, so its text keeps its size on a phone instead of scaling down with it
const REGION_HEX = ['#2a78d6', '#eda100'];
function rateChart(ms, by, title, W) {
  const bins = ms[0].meta.rates[by].map((r) => r.bin), bar = 12, gap = 2, band = ms.length * (bar + gap) + 14;
  const L = 74, Rm = 128, T = 8, H = T + bins.length * band + 22, sx = (v) => L + v * (W - L - Rm);
  const grid = [0, 0.5, 1].map((v) => `<line class="grid" x1="${sx(v)}" x2="${sx(v)}" y1="${T}" y2="${H - 22}"/>
    <text x="${sx(v)}" y="${H - 6}" text-anchor="middle">${v * 100} %</text>`).join('');
  const rows = bins.map((b, i) => {
    const y = T + i * band;
    return `<text x="${L - 10}" y="${y + ms.length * (bar + gap) / 2 + 4}" text-anchor="end">${esc(b)}</text>` + ms.map(({ region, meta }, k) => {
      const { recovered, total } = meta.rates[by][i], share = total ? recovered / total : 0, yb = y + k * (bar + gap);
      const w = Math.max(0, sx(share) - L), d = w > 4 ? `M${L} ${yb}h${w - 4}a4 4 0 0 1 4 4v${bar - 8}a4 4 0 0 1 -4 4h${-(w - 4)}z` : '';
      const text = total ? `${pct(share)} · ${recovered} of ${total}` : 'none';
      return `<g><title>${esc(meta.title)}, ${esc(b)}: ${recovered} of ${total} recovered</title>
        <path d="${d}" fill="${REGION_HEX[k]}"/><rect x="${L}" y="${yb}" width="${W - L - Rm}" height="${bar}" fill="transparent"/>
        <text class="value" x="${sx(share) + 6}" y="${yb + bar - 2}">${text}</text></g>`;
    }).join('');
  }).join('');
  return `<figure><figcaption>${title}</figcaption><svg class="rates-chart" viewBox="0 0 ${W} ${H}" role="img"
    aria-label="${title}: ${bins.map((b, i) => `${b}, ${ms.map(({ meta }) => `${meta.title} ${meta.rates[by][i].recovered} of ${meta.rates[by][i].total}`).join(', ')}`).join('; ')}">${grid}${rows}</svg></figure>`;
}

// One line per region over the sweep's runs, min_cluster_size on a log axis: `value` of each run, on a 0–100 %
// axis or, with `log`, a log one. Drawn W units wide for a figure W pixels wide, as rateChart
function sweepChart(ms, title, W, { value, text, log }) {
  const runs = ms[0].meta.runs.sweep.map((s) => s.run), L = 46, Rm = 64, T = 10, B = 34, H = 230;
  const lk = (k) => Math.log(k), sx = (k) => L + (lk(k) - lk(runs[0])) / (lk(runs.at(-1)) - lk(runs[0])) * (W - L - Rm);
  const all = ms.flatMap(({ meta }) => meta.runs.sweep.map((s) => value(s, meta)));
  const [y0, y1] = log ? [10 ** Math.floor(Math.log10(Math.min(...all))), 10 ** Math.ceil(Math.log10(Math.max(...all)))] : [0, 1];
  const f = log ? Math.log10 : (v) => v, sy = (v) => H - B - (f(v) - f(y0)) / (f(y1) - f(y0)) * (H - T - B);
  const yt = log ? Array.from({ length: Math.log10(y1 / y0) + 1 }, (_, i) => y0 * 10 ** i) : [0, 0.5, 1];
  const grid = yt.map((v) => `<line class="grid" x1="${L}" x2="${W - Rm}" y1="${sy(v)}" y2="${sy(v)}"/>
    <text x="${L - 8}" y="${sy(v) + 4}" text-anchor="end">${log ? num(v) : `${v * 100} %`}</text>`).join('');
  const xt = runs.filter((k) => [10, 20, 40, 80, 160, 240].includes(k)).map((k) =>
    `<text x="${sx(k)}" y="${H - B + 16}" text-anchor="middle">${k}</text>`).join('');
  const ends = ms.map(({ meta }) => sy(value(meta.runs.sweep.at(-1), meta))), order = ends.map((_, j) => j).sort((a, b) => ends[a] - ends[b]);
  order.forEach((j, n) => { if (n) ends[j] = Math.max(ends[j], ends[order[n - 1]] + 13); });  // end labels 13 units apart at least
  const lines = ms.map(({ meta }, j) => {
    const pts = meta.runs.sweep.map((s) => [sx(s.run), sy(value(s, meta)), s]);
    const [ex] = pts.at(-1), ey = ends[j];
    return `<path d="M${pts.map(([x, y]) => `${x.toFixed(1)} ${y.toFixed(1)}`).join('L')}" fill="none" stroke="${REGION_HEX[j]}" stroke-width="2"/>
      ${pts.map(([x, y, s]) => `<g><title>${esc(meta.title)}, min_cluster_size ${s.run}: ${text(s, meta)}</title>
        <circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="3.2" fill="${REGION_HEX[j]}"/><circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="9" fill="transparent"/></g>`).join('')}
      <text class="value" x="${ex + 8}" y="${ey + 4}">${esc(meta.title)}</text>`;
  }).join('');
  return `<figure><figcaption>${title}</figcaption><svg class="rates-chart" viewBox="0 0 ${W} ${H}" role="img"
    aria-label="${title}: ${ms.map(({ meta }) => `${meta.title}, ${meta.runs.sweep.map((s) => `${s.run}: ${text(s, meta)}`).join(', ')}`).join('; ')}">
    ${grid}${xt}<text x="${(L + W - Rm) / 2}" y="${H - 3}" text-anchor="middle">min_cluster_size</text>${lines}</svg></figure>`;
}

// The kiloparsec past the field's 99th percentile and the farthest cluster: some clusters lie beyond the field
const farEdge = (meta) => Math.ceil(Math.max(meta.reach, ...meta.clusters.map((c) => Math.hypot(...c.xyz))) / 1000) * 1000;

// Where the regions sit in the Milky Way, seen from the north Galactic pole: a schematic disc in kiloparsecs,
// the centre at the origin and the Sun 8.2 kpc below it. Longitude l points at (−sin l, −cos l) on screen.
const R0 = 8.2, DISC = 15;
// PARSEC isochrones for the cards' diagrams (export_viewer.py), shared by the regions; none if not exported
const isochrones = fetch('data/isochrones.json').then((r) => r.ok ? r.json() : null).catch(() => null);
let ISO = null;
isochrones.then((v) => { ISO = v; if (R) showCard(hovered ?? pinned); });
const metas = Promise.all([...document.querySelectorAll('[data-region]')].map(async ({ dataset: { region } }) =>
  ({ region, meta: await (await fetch(`data/${region}/meta.json`)).json() })));
const wedges = metas.then((ms) => ms.map(({ region, meta: m }) => ({ region, l: m.box.l, far: farEdge(m) / 1000 })));
metas.then((ms) => {
  $('#rates-key').innerHTML = ms.map(({ meta }, k) => `<span><i style="background:${REGION_HEX[k]}"></i>${esc(meta.title)}</span>`).join('');
  const charts = $('#rates-charts');
  let drawn = 0;
  new ResizeObserver(() => {  // the columns of .rates-charts: as many 300 px ones as fit, at most two, 48 px apart
    const cw = charts.clientWidth, cols = Math.max(1, Math.min(2, Math.floor((cw + 48) / 348)));
    const W = Math.round((cw - 48 * (cols - 1)) / cols);
    if (W === drawn) return;
    drawn = W;
    charts.innerHTML = rateChart(ms, 'distance', 'By the cluster\'s distance', W) + rateChart(ms, 'members', 'By its members in our sample', W);
    if (ms.every(({ meta }) => meta.runs?.sweep.length)) {
      $('#sweep').hidden = false;
      $('#sweep-charts').innerHTML = sweepChart(ms, 'Recovered', W, { value: (s, m) => s.recovered / m.recoverable,
          text: (s, m) => `${s.recovered} of ${m.recoverable} recovered` })
        + sweepChart(ms, 'Groups no Hunt &amp; Reffert cluster owns', W, { log: true, value: (s) => s.unverified,
          text: (s) => `${num(s.unverified)} of ${num(s.groups)} groups` });
    }
  }).observe(charts);
}).catch(() => { $('#rates').hidden = true; });
async function galaxy(name) {
  const at = (l, d) => `${(-Math.sin(l * RAD) * d).toFixed(2)} ${(R0 - Math.cos(l * RAD) * d).toFixed(2)}`;
  const shapes = (await wedges).map(({ region, l: [l0, l1], far }) =>
    `<path class="${region === name ? 'here' : 'other'}" d="M0 ${R0}L${at(l0, far)}A${far} ${far} 0 0 0 ${at(l1, far)}Z"/>`);
  const el = $('#galaxy');
  el.innerHTML = `<circle class="disc" r="${DISC}"/><path class="centre" d="M-0.8 0H0.8M0 -0.8V0.8"/>
    <text x="0" y="-1.6" text-anchor="middle">Galactic centre</text>${shapes.join('')}
    <circle class="sun" cy="${R0}" r="0.55"/><text x="-1.2" y="${R0 + 0.9}" text-anchor="end">Sun</text>
    <text class="title" x="0" y="-15.3" text-anchor="middle">Milky Way, from above</text>`;
  el.setAttribute('aria-label', `Where ${name[0].toUpperCase() + name.slice(1)} sits in the Milky Way, seen from above: a wedge from the Sun, ${R0} kpc from the Galactic centre`);
}

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
  const lab = { el, at, distance, rank: kind.includes('sun') ? 0 : kind === 'named' ? 1 : 2 };
  region.labels.push(lab);
  region.labels.sort((a, b) => a.rank - b.rank);  // placeLabels() hides a label that would cover one placed before it
  return lab;
}

const projected = new THREE.Vector3();
function placeLabels() {
  const w = box.clientWidth, h = box.clientHeight;
  const at = box.getBoundingClientRect();
  const taken = [$('#hint'), $('.zoom')].map((el) => {
    const r = el.getBoundingClientRect();
    return [r.left - at.left, r.top - at.top, r.right - at.left, r.bottom - at.top];
  });
  for (const lab of R.labels) {
    drifted(lab.at, lab.pm, projected).project(camera);
    let off = projected.z > 1 || Math.abs(projected.x) > 1.05 || Math.abs(projected.y) > 1.05;
    if (!off) {
      lab.w ||= lab.el.offsetWidth;
      lab.h ||= lab.el.offsetHeight;
      const px = (projected.x + 1) / 2 * w, y = (1 - projected.y) / 2 * h + (lab.below ? 4 : -15);
      const x = px + 7 + lab.w > w ? px - 7 - lab.w : px + 7;  // to the left where the right edge would cut it
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

let playing = null;  // the time and clock time the play started from
// Where the pinned cluster or group is now, if the view is close enough to be looking at it
function followed() {
  if (!R || !pinned || flight || camera.position.distanceTo(controls.target) > 300) return null;
  if (pinned.kind === 'cluster') return clusterAt(R.meta.clusters[pinned.i]);
  const { xyz, pm } = R.view.data.unverified;
  return drifted(new THREE.Vector3().fromArray(xyz, 3 * pinned.i), [pm[2 * pinned.i], pm[2 * pinned.i + 1]]);
}
// Change where things are drawn, with the view travelling along with the pinned cluster, so it stays in sight
function moving(change) {
  const before = followed();
  change();
  if (before) {
    const by = followed().sub(before);
    controls.target.add(by);
    camera.position.add(by);
  }
  needs = true;
}
function setTime(myr) {
  moving(() => { motion.myr.value = Math.max(-SPAN, Math.min(SPAN, myr)); });
  $('#myr').value = motion.myr.value;
  const t = motion.myr.value;
  $('#myr-out').textContent = Math.abs(t) < 0.005 ? 'now' : `${t > 0 ? '+' : '\u2212'}${Math.abs(t).toFixed(2)} Myr`;
  $('#now').disabled = t === 0;
  lastMove = performance.now();  // draws fewer field stars while the clock moves, as while the view moves
  needs = true;
}
async function play() {
  if (playing) return stop();
  const region = R;
  $('#play').textContent = 'Loading motions…';
  await loadMotions(region);
  if (region !== R) return;
  if (motion.myr.value >= SPAN) setTime(-SPAN);
  playing = { from: motion.myr.value, t0: performance.now() };
  $('#play').textContent = 'Pause';
  $('#play').setAttribute('aria-pressed', 'true');
}
function stop() {
  playing = null;
  $('#play').textContent = 'Play';
  $('#play').setAttribute('aria-pressed', 'false');
}
function tick() {
  const t = playing.from + (performance.now() - playing.t0) / 1000;  // a million years a second
  setTime(t);
  if (t >= SPAN) stop();
}

// The guided tour (tour.js): each stop pins its cluster and flies to it, or shows the whole region
let stop_ = -1;
function tour(k) {
  const stops = TOURS[R?.name] ?? [];
  if (k < 0 || k >= stops.length) return endTour();
  stop_ = k;
  if (R.view !== R.score) { $('#one-run').checked = false; chooseRun(); }  // the stops' numbers are the score's
  const s = stops[k], i = s.cluster ? R.meta.clusters.findIndex((c) => c.name === s.cluster) : -1;
  if (i >= 0) {  // from a few hundred parsecs, so its neighbours show too
    pin({ kind: 'cluster', i });
    const at = centre(R, i);
    flyTo(drifted(at, R.meta.clusters[i].pm), Math.max(4 * viewDistance(R, i, at), 0.6 * R.meta.clusters[i].dist));
    $('#find').value = String(i);
  } else {
    pin(null);
    $('#find').value = '';
    zoomTo('home');
  }
  $('#tour-text').textContent = s.text(R.meta.clusters[i], R.meta, { kpc, num, pct, age, purityMisses }).replace(/\s+/g, ' ');
  $('#tour-step').textContent = `${k + 1} of ${stops.length}`;
  $('#tour-back').disabled = k === 0;
  $('#tour-next').textContent = k === stops.length - 1 ? 'Done' : 'Next';
  $('#tour').hidden = false;
  $('#tour-next').focus({ preventScroll: true });
}
// Missed clusters whose closest group holds most of their members: they fail on purity
const purityMisses = (m) => m.clusters.filter((c) => !c.recovered && c.recovery > 0.5).length;
function endTour() {
  stop_ = -1;
  $('#tour').hidden = true;
}

function clear() {
  if (!R) return;
  endTour();
  stop();
  setTime(0);
  scene.remove(R.group);
  fieldScene.remove(R.field);
  for (const g of [R.group, R.field]) g.traverse((o) => o.geometry?.dispose());
  for (const v of [R.score, ...R.views.values()]) {
    Promise.resolve(v).then((v) => [v.members.points, v.rings, v.unverified].forEach((o) => o.geometry.dispose())).catch(() => {});
  }
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
  for (const lab of region.labels) {  // under the edge across a wide view, beside it up a tall one: off the clusters
    if (lab.distance) lab.at.copy(toward(l, b, lab.distance)), lab.below = !tall;
  }
  needs = true;
}

// Where a cluster is drawn, and how close to come to see it whole
function centre(region, i) {
  const m = region.view.clusters[i].members;
  if (!m) return new THREE.Vector3().fromArray(region.meta.clusters[i].xyz);
  const c = new THREE.Vector3(), p = region.view.members.pos;
  for (let j = m[0]; j < m[0] + m[1]; j++) c.x += p[3 * j], c.y += p[3 * j + 1], c.z += p[3 * j + 2];
  return c.divideScalar(m[1]);
}

function viewDistance(region, i, at) {
  const m = region.view.clusters[i].members;
  if (!m) return 40;
  const p = region.view.members.pos, d = [];
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
  for (const i of [...R.view.found, ...R.view.missed]) {
    const c = R.meta.clusters[i];
    const p = clusterAt(c);
    consider(dist2(p.x, p.y, p.z), screenRadius(c) + 4, { kind: 'cluster', i });
  }
  if (R.view.unverified.visible) {
    const { xyz, pm } = R.view.data.unverified, p = new THREE.Vector3();
    for (let i = 0; i < xyz.length / 3; i++) {
      drifted(p.fromArray(xyz, 3 * i), [pm[2 * i], pm[2 * i + 1]], p);
      consider(dist2(p.x, p.y, p.z), 8, { kind: 'group', i });
    }
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
    flyTo(drifted(at, R.meta.clusters[hit.i].pm), viewDistance(R, hit.i, at));
  } else {
    const { xyz, pm } = R.view.data.unverified;
    flyTo(drifted(new THREE.Vector3().fromArray(xyz, 3 * hit.i), [pm[2 * hit.i], pm[2 * hit.i + 1]]), 40);
  }
}

function showCard(hit) {
  if (!R) return;
  card.innerHTML = !hit ? emptyCard() : hit.kind === 'cluster' ? clusterCard(hit.i) : groupCard(hit.i);
}

function emptyCard() {
  const m = R.meta, cs = R.view.clusters, rec = cs.filter((c) => c.recovered).length;
  const rel = cs.filter((c) => c.reliable), relRec = rel.filter((c) => c.recovered).length;
  return `<div class="info"><p class="kicker">${esc(m.title)}${R.view.run === null ? '' : ` · min_cluster_size ${R.view.run} alone`}</p>
    <p>Point at a blue sphere or a red ring, or tap it, for its card. Every cluster is also in the list under the map.</p>
    <dl class="facts"><dt>Recovered</dt><dd>${rec} of ${cs.length} (${pct(rec / cs.length)})</dd>
    <dt>Highly reliable</dt><dd>${relRec} of ${rel.length} (${pct(relRec / rel.length)})</dd>
    <dt>Missed</dt><dd>${cs.length - rec}</dd><dt>Field stars</dt><dd>${num(m.field.count)}</dd></dl></div>`;
}

function clusterCard(i) {
  const c = R.view.clusters[i], plain = c.name.replaceAll('_', ' ');
  const status = c.recovered ? '<p class="kicker rec">Recovered</p>' : `<p class="kicker miss">Missed · ${esc(c.how.replace(/ \(.*\)$/, ''))}</p>`;
  const others = c.runs.filter((k) => k !== c.run);
  const group = !c.size ? `No group of ${R.view.run === null ? 'the leaf runs' : 'this run'} holds any of its members.`
    : `${c.recovered ? 'Shown' : 'Closest'}: a group of ${num(c.size)} stars in the run with min_cluster_size ${c.run}.`
      + (others.length ? ` Also recovered at ${list(others)}.` : '');
  const figure = c.cmd[1] ? `<div class="plots"><figure>${cmdFigure(c)}<figcaption>Colour and brightness</figcaption></figure>
    <figure>${motionFigure(c)}<figcaption>Motion on the sky</figcaption></figure></div>
    <div class="cmd-key"><span><i style="background:var(--recovered)"></i>in both</span><span><i style="background:var(--missed)"></i>ours only</span><span><i style="background:var(--theirs)"></i>theirs only</span><span><i style="background:var(--near)"></i>field around it</span>${c.fit && ISO ? '<span><b class="iso-key"></b>our isochrone</span>' : ''}${c.hr && ISO ? '<span><b class="iso-key theirs"></b>Hunt &amp; Reffert\'s</span>' : ''}</div>` : '';
  return `<div class="info"><div class="head"><div>${status}<h3>${esc(c.label)}</h3>${c.label !== plain ? `<span class="alias">${esc(plain)} in Hunt &amp; Reffert</span>` : ''}</div>
    ${same(pinned, { kind: 'cluster', i }) ? '<button type="button" class="close" data-act="unpin" aria-label="Close">×</button>' : ''}</div>
    <p class="meta">${kpc(c.dist)} · ${num(c.n)} members in our sample · ${age(c)}</p>
    <div class="scores"><div><span>Recovery</span><b>${pct(c.recovery)}</b></div><div><span>Purity</span><b>${pct(c.purity)}</b></div></div>
    <p>${group} ${c.reliable ? 'Highly reliable in Hunt &amp; Reffert.' : 'Not rated highly reliable by Hunt &amp; Reffert.'}</p>
    ${ages(c)}
    <button type="button" class="go" data-act="go" data-i="${i}">Go to ${esc(c.label)}</button></div>
    ${figure ? `<div class="plot">${figure}</div>` : ''}`;
}

// Our isochrone age beside Hunt & Reffert's (docs/protocol-ages.md), if the cluster was fitted
function ages(c) {
  const f = c.fit;
  if (!f) return '';
  const ours = f.given ? years(f.age, f.range) : `none (${esc(f.why)})`;
  return `<dl class="facts"><dt>Age, Hunt &amp; Reffert</dt><dd>${years(c.age[1], [c.age[0], c.age[2]])}</dd>
    <dt>Age, our isochrone fit</dt><dd>${ours}${f.calibration ? ' · calibration cluster' : ''}</dd></dl>`;
}
const years = (mid, [lo, hi]) => {
  const [a, b, m] = [lo, hi, mid].map((v) => 10 ** v / 1e6), gyr = m >= 1000, f = (v) => gyr ? (v / 1000).toFixed(1) : String(+v.toPrecision(2));
  return `${f(m)} ${gyr ? 'Gyr' : 'Myr'} (${f(a)}–${f(b)})`;
};

// The age in Hunt & Reffert: median, with the 16th and 84th percentiles
const age = (c) => {
  const [lo, mid, hi] = c.age.map((v) => 10 ** v / 1e6), gyr = mid >= 1000, f = (v) => gyr ? (v / 1000).toFixed(1) : String(+v.toPrecision(2));
  return `${f(mid)} ${gyr ? 'Gyr' : 'Myr'} old (${f(lo)}–${f(hi)})`;
};

const list = (xs) => xs.length < 2 ? `${xs[0]}` : `${xs.slice(0, -1).join(', ')} and ${xs.at(-1)}`;

function groupCard(i) {
  const u = R.view.data.unverified, at = new THREE.Vector3().fromArray(u.xyz, 3 * i);
  return `<div class="info"><div class="head"><div><p class="kicker">Unverified group</p><h3>${num(u.n[i])} stars</h3></div></div>
    <p class="meta">${kpc(at.length())}, the median of its stars · min_cluster_size ${R.view.unverifiedRun}</p>
    <p>No Hunt &amp; Reffert cluster holds most of its stars. Hunt &amp; Reffert discard groups like this with a significance
    test and a colour-magnitude classifier, which this project does not have, so it is not counted as a cluster.</p></div>`;
}

// A small scatter plot: `layers` of [class, xs, ys] in drawing order, inside [x0, x1] × [y0, y1]; y grows
// downwards unless `up`
let clips = 0;
function scatter({ x0, x1, y0, y1, up, layers, lines = [], xlabel, ylabel, label, W = 300, H = 250 }) {
  const L = 34, B = 38, T = 6, Rm = 6;
  const sx = (x) => L + (x - x0) / (x1 - x0) * (W - L - Rm);
  const sy = (y) => up ? H - B - (y - y0) / (y1 - y0) * (H - T - B) : T + (y - y0) / (y1 - y0) * (H - T - B);
  const count = layers.reduce((n, [, xs]) => n + xs.length, 0), r = count > 1500 ? 1.1 : count > 300 ? 1.5 : 2;
  const paths = layers.map(([cls, xs, ys]) => {
    let d = '';
    xs.forEach((v, j) => {
      if (v < x0 || v > x1 || ys[j] < y0 || ys[j] > y1) return;
      const x = sx(v), y = sy(ys[j]);
      d += `M${(x - r).toFixed(1)} ${y.toFixed(1)}a${r} ${r} 0 1 0 ${2 * r} 0a${r} ${r} 0 1 0 ${-2 * r} 0`;
    });
    return `<path class="${cls}" d="${d}"/>`;
  }).join('');
  const tx = ticks(x0, x1, 4).map((v) => `<line class="axis" x1="${sx(v)}" x2="${sx(v)}" y1="${H - B}" y2="${H - B + 4}"/>
    <text x="${sx(v)}" y="${H - B + 15}" text-anchor="middle">${+v.toFixed(2)}</text>`).join('');
  const ty = ticks(y0, y1, 5).map((v) => `<line class="axis" x1="${L - 4}" x2="${L}" y1="${sy(v)}" y2="${sy(v)}"/>
    <text x="${L - 7}" y="${sy(v) + 3.5}" text-anchor="end">${+v.toFixed(1)}</text>`).join('');
  const clip = `clip${++clips}`;  // the lines stop at the axes
  const curves = lines.map(([cls, xs, ys]) => `<path class="${cls}" clip-path="url(#${clip})" d="M${xs.map((v, j) =>
    `${sx(v).toFixed(1)} ${sy(ys[j]).toFixed(1)}`).join('L')}"/>`).join('');
  return `<svg class="cmd" viewBox="0 0 ${W} ${H}" role="img" aria-label="${label}">
      <clipPath id="${clip}"><rect x="${L}" y="${T}" width="${W - L - Rm}" height="${H - T - B}"/></clipPath>${paths}${curves}
      <line class="axis" x1="${L}" x2="${L}" y1="${T}" y2="${H - B}"/><line class="axis" x1="${L}" x2="${W - Rm}" y1="${H - B}" y2="${H - B}"/>
      ${tx}${ty}<text x="${(L + W - Rm) / 2}" y="${H - 4}" text-anchor="middle">${xlabel}</text><text x="${L + 6}" y="${T + 9}">${ylabel}</text>
    </svg>`;
}

// The stars of a cluster's diagrams by category, 0 in both, 1 ours only, 2 theirs only (drawn 2, 0, 1)
function byCategory(c, x, y) {
  const [start, count] = c.cmd, out = [[[], []], [[], []], [[], []]];
  for (let j = start; j < start + count; j++) out[R.view.cmd.cat[j]][0].push(x(j)), out[R.view.cmd.cat[j]][1].push(y(j));
  return [['theirs', ...out[2]], ['both', ...out[0]], ['ours', ...out[1]]];
}

// The isochrones of our fit (solid) and of Hunt & Reffert's (dashed), reddened by A_V and shifted by the
// distance modulus as the fit has them: the grid age nearest to theirs
function isochroneLines(c) {
  const e = R.meta.extinction;
  if (!ISO || !e) return [];
  const line = (cls, age, av, dm) => {
    const k = ISO.ages.reduce((b, a, i) => Math.abs(a - age) < Math.abs(ISO.ages[b] - age) ? i : b, 0), p = ISO.curves[k];
    const xs = [], ys = [];
    for (let j = 0; j < p.length; j += 2) xs.push(p[j] + (e.BP - e.RP) * av), ys.push(p[j + 1] + dm + e.G * av);
    return [cls, xs, ys];
  };
  return [c.hr && line('iso-theirs', c.age[1], c.hr.av, c.hr.dm), c.fit && line('iso-ours', c.fit.age, c.fit.av, c.fit.dm)].filter(Boolean);
}

// Colour-magnitude diagram: G against BP − RP, blue in both, red ours only, grey theirs only
function cmdFigure(c) {
  const { bp, g } = R.view.cmd, layers = byCategory(c, (j) => bp[j] / 1000, (j) => g[j] / 1000);
  const xs = layers.flatMap((l) => l[1]), ys = layers.flatMap((l) => l[2]), n = layers.map((l) => l[1].length);
  return scatter({ x0: Math.min(...xs) - 0.1, x1: Math.max(...xs) + 0.1, y0: Math.min(...ys) - 0.3, y1: Math.max(...ys) + 0.3,
    layers, lines: isochroneLines(c), xlabel: 'BP − RP', ylabel: 'G',
    label: `Colour-magnitude diagram of ${esc(c.label)}: ${n[1]} stars in both, ${n[2]} in our group only, ${n[0]} Hunt and Reffert members it misses` });
}

// Proper-motion diagram: the same stars as the colour-magnitude diagram, over the field that shares the cluster's
// patch of sky and distance (light grey). Square, centred on the cluster's motion, wide enough for most of the field
function motionFigure(c) {
  const step = R.meta.pm_step, pm = R.view.cmd.pm, [o, k] = c.near, [cl, cb] = c.pm;
  const field = [[], []];
  for (let j = o; j < o + k; j++) field[0].push(R.near[2 * j] * step), field[1].push(R.near[2 * j + 1] * step);
  const off = field[0].map((v, j) => Math.max(Math.abs(v - cl), Math.abs(field[1][j] - cb))).sort((a, b) => a - b);
  const half = Math.max(1, 1.1 * (off[Math.floor(0.9 * (off.length - 1))] ?? 0));
  const [theirs, both, ours] = byCategory(c, (j) => pm[2 * j] * step, (j) => pm[2 * j + 1] * step);
  const layers = [['near', ...field], theirs, ours, both];  // the field stars in our group move with it: under the rest
  return scatter({ x0: cl - half, x1: cl + half, y0: cb - half, y1: cb + half, up: true, W: 250, H: 250,
    layers, xlabel: 'μl cos b (mas/yr)', ylabel: 'μb',
    label: `Proper-motion diagram of ${esc(c.label)}: its stars against ${k} field stars in the same patch of sky and distance` });
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

function fillFind({ clusters }) {
  const order = (a, b) => clusters[a].label.localeCompare(clusters[b].label, 'en', { numeric: true });
  const all = clusters.map((_, i) => i), keep = $('#find').value;
  const option = (i) => `<option value="${i}">${esc(clusters[i].label)}</option>`;
  const rec = all.filter((i) => clusters[i].recovered).sort(order), miss = all.filter((i) => !clusters[i].recovered).sort(order);
  $('#find').innerHTML = `<option value="">${all.length} clusters</option>
    <optgroup label="Recovered (${rec.length})">${rec.map(option).join('')}</optgroup>
    <optgroup label="Missed (${miss.length})">${miss.map(option).join('')}</optgroup>`;
  $('#find').value = keep;
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
$('#by-age').addEventListener('change', () => R && paint(R));
$('#unverified').addEventListener('change', (e) => {
  if (!R) return;
  R.view.unverified.visible = e.target.checked;
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
$('#play').addEventListener('click', () => play().catch(failed));
$('#now').addEventListener('click', () => { stop(); setTime(0); });
$('#sun-frame').addEventListener('change', (e) => moving(() => { motion.relative.value = e.target.checked ? 0 : 1; }));
$('#myr').min = -SPAN;
$('#myr').max = SPAN;
$('#myr').addEventListener('input', async (e) => {
  stop();
  const region = R, t = +e.target.value;
  await loadMotions(region).catch(failed);
  if (region === R) setTime(t);
});
card.addEventListener('click', (e) => {
  const act = e.target.closest('[data-act]')?.dataset.act;
  if (act === 'unpin') { pin(null); $('#find').value = ''; }
  if (act === 'go') goTo({ kind: 'cluster', i: +e.target.closest('[data-act]').dataset.i });
});
addEventListener('keydown', (e) => {
  if (stop_ >= 0 && e.key === 'Escape') return endTour();
  if (stop_ >= 0 && (e.key === 'ArrowRight' || e.key === 'ArrowLeft') && !e.target.closest('input, select')) {
    return tour(stop_ + (e.key === 'ArrowRight' ? 1 : -1));
  }
  if (e.key === 'Escape' && pinned) { pin(null); $('#find').value = ''; }
});
$('#one-run').addEventListener('change', () => chooseRun().catch(failed));
$('#run').addEventListener('input', () => chooseRun().catch(failed));
$('#tour-start').addEventListener('click', () => tour(0));
$('#tour-next').addEventListener('click', () => tour(stop_ + 1));
$('#tour-back').addEventListener('click', () => tour(stop_ - 1));
$('#tour-close').addEventListener('click', endTour);

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
  if (playing) tick();
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
