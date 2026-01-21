import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

// =====================
// Constants / helpers
// =====================
const R_EARTH_KM = 6378.137;
const MU_EARTH = 398600.4418; // km^3/s^2
const deg2rad = (d) => d * Math.PI / 180;
const rad2deg = (r) => r * 180 / Math.PI;

const kmToScene = (km) => km / R_EARTH_KM;       // Earth radius = 1 scene unit
const rKmToScene = (rKm) => rKm / R_EARTH_KM;

function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }

function safeNum(v, fallback = 0) {
  const n = Number(v);
  return Number.isFinite(n) ? n : fallback;
}

// =====================
// UI wiring
// =====================
const ui = {
  // buttons
  btnStart: document.getElementById("btnStart"),
  btnPause: document.getElementById("btnPause"),
  btnReset: document.getElementById("btnReset"),
  btnApply: document.getElementById("btnApply"),
  btnPlan:  document.getElementById("btnPlan"),

  // speed
  spd: document.getElementById("spd"),
  spdTxt: document.getElementById("spdTxt"),

  // target inputs
  tPreset: document.getElementById("tPreset"),
  tPerAlt: document.getElementById("tPerAlt"),
  tApoAlt: document.getElementById("tApoAlt"),
  tInc: document.getElementById("tInc"),
  tRAAN: document.getElementById("tRAAN"),
  tArgp: document.getElementById("tArgp"),
  tNu0: document.getElementById("tNu0"),

  // chaser inputs
  cPreset: document.getElementById("cPreset"),
  cPerAlt: document.getElementById("cPerAlt"),
  cApoAlt: document.getElementById("cApoAlt"),
  cInc: document.getElementById("cInc"),
  cRAAN: document.getElementById("cRAAN"),
  cArgp: document.getElementById("cArgp"),
  cNu0: document.getElementById("cNu0"),

  // beta
  beta: document.getElementById("beta"),
  betaTxt: document.getElementById("betaTxt"),

  // dv output
  dvOut: document.getElementById("dvOut"),
};

ui.spdTxt.textContent = `${ui.spd.value}x`;
ui.spd.addEventListener("input", () => ui.spdTxt.textContent = `${ui.spd.value}x`);
ui.betaTxt.textContent = ui.beta.value;
ui.beta.addEventListener("input", () => ui.betaTxt.textContent = ui.beta.value);

// =====================
// Scene / camera / renderer
// =====================
const scene = new THREE.Scene();
scene.fog = new THREE.Fog(0x070b12, 3, 25);

const camera = new THREE.PerspectiveCamera(55, innerWidth / innerHeight, 0.01, 500);
camera.position.set(3.8, 2.6, 3.8);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(innerWidth, innerHeight);
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
document.body.appendChild(renderer.domElement);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.target.set(0, 0, 0);

// Lights
scene.add(new THREE.AmbientLight(0xffffff, 0.45));
const key = new THREE.DirectionalLight(0xffffff, 1.0);
key.position.set(5, 8, 3);
scene.add(key);

// Planet (Earth)
const earth = new THREE.Mesh(
  new THREE.SphereGeometry(1, 64, 48),
  new THREE.MeshStandardMaterial({ color: 0x2d6cdf, roughness: 0.85, metalness: 0.0 })
);
earth.renderOrder = 0;
scene.add(earth);

// Equatorial plane grid (visual reference)
const planeRadius = rKmToScene(R_EARTH_KM + 35786) * 1.2;
const eqDisk = new THREE.Mesh(
  new THREE.CircleGeometry(planeRadius, 128),
  new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.05, side: THREE.DoubleSide })
);
eqDisk.rotation.x = 0;
eqDisk.renderOrder = 0;
scene.add(eqDisk);

const grid = new THREE.PolarGridHelper(planeRadius, 24, 12, 128, 0xffffff, 0xffffff);
grid.material.transparent = true;
grid.material.opacity = 0.18;
scene.add(grid);

const axes = new THREE.AxesHelper(2.2);
axes.material.transparent = true;
axes.material.opacity = 0.6;
scene.add(axes);

// =====================
// Orbit math (ECI-ish)
// r_eci = Rz(Ω) * Rx(i) * Rz(ω) * r_pf
// =====================
function rotateOrbitPoint(x, y, iRad, raanRad, argpRad) {
  let X = x, Y = y, Z = 0;

  // Rz(argp)
  const c3 = Math.cos(argpRad), s3 = Math.sin(argpRad);
  const x1 = c3 * X - s3 * Y;
  const y1 = s3 * X + c3 * Y;
  const z1 = Z;

  // Rx(i)
  const c1 = Math.cos(iRad), s1 = Math.sin(iRad);
  const x2 = x1;
  const y2 = c1 * y1 - s1 * z1;
  const z2 = s1 * y1 + c1 * z1;

  // Rz(raan)
  const c2 = Math.cos(raanRad), s2 = Math.sin(raanRad);
  const x3 = c2 * x2 - s2 * y2;
  const y3 = s2 * x2 + c2 * y2;
  const z3 = z2;

  return new THREE.Vector3(x3, y3, z3);
}

function positionOnOrbit(orbit, nuRad) {
  const a = rKmToScene(orbit.a_km);
  const e = orbit.e;
  const iRad = deg2rad(orbit.i_deg);
  const raanRad = deg2rad(orbit.raan_deg);
  const argpRad = deg2rad(orbit.argp_deg);

  const p = a * (1 - e * e);
  const r = p / (1 + e * Math.cos(nuRad));
  const x = r * Math.cos(nuRad);
  const y = r * Math.sin(nuRad);
  return rotateOrbitPoint(x, y, iRad, raanRad, argpRad);
}

function buildOrbitLine(orbit, material, segments = 512) {
  const pts = [];
  for (let k = 0; k <= segments; k++) {
    const nu = (k / segments) * Math.PI * 2;
    pts.push(positionOnOrbit(orbit, nu));
  }
  const geom = new THREE.BufferGeometry().setFromPoints(pts);
  const line = new THREE.Line(geom, material);
  return line;
}

// =====================
// Visual objects (lines, markers)
// =====================
const chaserMat = new THREE.LineBasicMaterial({ color: 0xffd54a, transparent: true, opacity: 0.95 });
const targetMat = new THREE.LineBasicMaterial({ color: 0x3fdcff, transparent: true, opacity: 0.95 });

// Transfer: dashed & on-top
const transferMat = new THREE.LineDashedMaterial({
  color: 0xff4fd8,
  dashSize: 0.12,
  gapSize: 0.08,
  transparent: true,
  opacity: 1.0
});
transferMat.depthTest = false;         // IMPORTANT: don't let Earth hide it

let chaserLine = null;
let targetLine = null;
let transferLine = null;

// Burn markers (white, on-top)
const burnGeo = new THREE.SphereGeometry(0.05, 16, 12);
const burnMat = new THREE.MeshBasicMaterial({ color: 0xffffff });
burnMat.depthTest = false;             // IMPORTANT: don't let Earth hide it
const burn1 = new THREE.Mesh(burnGeo, burnMat);
const burn2 = new THREE.Mesh(burnGeo, burnMat);
burn1.visible = false;
burn2.visible = false;
burn1.renderOrder = 11;
burn2.renderOrder = 11;
scene.add(burn1, burn2);

// Spacecraft markers
const chaserCraft = new THREE.Mesh(
  new THREE.SphereGeometry(0.04, 16, 12),
  new THREE.MeshBasicMaterial({ color: 0xffd54a })
);
const targetCraft = new THREE.Mesh(
  new THREE.SphereGeometry(0.04, 16, 12),
  new THREE.MeshBasicMaterial({ color: 0x3fdcff })
);
scene.add(chaserCraft, targetCraft);

// Sun direction (visual only)
const sunArrow = new THREE.ArrowHelper(new THREE.Vector3(1,0,0), new THREE.Vector3(0,0,0), 3.0, 0xffffaa, 0.18, 0.10);
sunArrow.renderOrder = 2;
scene.add(sunArrow);

// =====================
// Element presets
// =====================
function presetTarget(name) {
  switch (name) {
    case "ISS":
      return { perAlt_km: 400, apoAlt_km: 400, i_deg: 51.6, raan_deg: 0, argp_deg: 0, nu0_deg: 0 };
    case "GEO":
      return { perAlt_km: 35786, apoAlt_km: 35786, i_deg: 0, raan_deg: 0, argp_deg: 0, nu0_deg: 0 };
    case "SSO":
      return { perAlt_km: 600, apoAlt_km: 600, i_deg: 98.0, raan_deg: 40, argp_deg: 0, nu0_deg: 0 };
    case "MOLNIYA":
      // typical-ish: a ~ 26600km, e ~0.74 -> rp ~ 6930km (alt ~552), ra ~ 46270km (alt ~39892)
      return { perAlt_km: 600, apoAlt_km: 39892, i_deg: 63.4, raan_deg: 20, argp_deg: 270, nu0_deg: 0 };
    default:
      return null;
  }
}

function presetChaser(name) {
  switch (name) {
    case "LEO100":
      return { perAlt_km: 100, apoAlt_km: 100, i_deg: 51.6, raan_deg: 0, argp_deg: 0, nu0_deg: -40 };
    case "LEO200":
      return { perAlt_km: 200, apoAlt_km: 200, i_deg: 51.6, raan_deg: 0, argp_deg: 0, nu0_deg: -40 };
    case "GTO":
      return { perAlt_km: 250, apoAlt_km: 35786, i_deg: 28.5, raan_deg: 0, argp_deg: 0, nu0_deg: -20 };
    default:
      return null;
  }
}

// =====================
// Convert UI -> orbit elements
// =====================
function elementsFromAlt(perAlt_km, apoAlt_km, i_deg, raan_deg, argp_deg, nu0_deg) {
  const rp_km = R_EARTH_KM + perAlt_km;
  const ra_km = R_EARTH_KM + apoAlt_km;
  const a_km = 0.5 * (rp_km + ra_km);
  const e = clamp((ra_km - rp_km) / (ra_km + rp_km), 0, 0.999999);

  return {
    a_km,
    e,
    i_deg,
    raan_deg,
    argp_deg,
    nu0_deg
  };
}

function getTargetOrbitFromUI() {
  return elementsFromAlt(
    safeNum(ui.tPerAlt.value, 400),
    safeNum(ui.tApoAlt.value, 400),
    safeNum(ui.tInc.value, 0),
    safeNum(ui.tRAAN.value, 0),
    safeNum(ui.tArgp.value, 0),
    safeNum(ui.tNu0.value, 0)
  );
}

function getChaserOrbitFromUI() {
  return elementsFromAlt(
    safeNum(ui.cPerAlt.value, 200),
    safeNum(ui.cApoAlt.value, 200),
    safeNum(ui.cInc.value, 0),
    safeNum(ui.cRAAN.value, 0),
    safeNum(ui.cArgp.value, 0),
    safeNum(ui.cNu0.value, -40)
  );
}

// =====================
// Draw / update orbits
// =====================
let state = {
  running: true,
  tNu: 0, // rad
  cNu: 0, // rad
  targetOrbit: getTargetOrbitFromUI(),
  chaserOrbit: getChaserOrbitFromUI(),
  transferOrbit: null,
  planned: false
};

function removeIf(line) { if (line) scene.remove(line); }

function applyOrbits() {
  state.targetOrbit = getTargetOrbitFromUI();
  state.chaserOrbit = getChaserOrbitFromUI();

  // Set starting anomalies from UI nu0
  state.tNu = deg2rad(state.targetOrbit.nu0_deg);
  state.cNu = deg2rad(state.chaserOrbit.nu0_deg);

  // Rebuild orbit lines
  removeIf(targetLine);
  removeIf(chaserLine);

  targetLine = buildOrbitLine(state.targetOrbit, targetMat, 720);
  chaserLine = buildOrbitLine(state.chaserOrbit, chaserMat, 720);

  targetLine.renderOrder = 3;
  chaserLine.renderOrder = 3;

  scene.add(targetLine, chaserLine);

  // Hide transfer until re-planned
  if (transferLine) {
    scene.remove(transferLine);
    transferLine = null;
  }
  burn1.visible = false;
  burn2.visible = false;
  state.transferOrbit = null;
  state.planned = false;

  ui.dvOut.textContent = "Orbits applied.\nClick “Plan Maneuvers (ΔV)” to compute transfer + burns.";
}

ui.btnApply.addEventListener("click", applyOrbits);

// =====================
// Transfer planning + ΔV (baseline)
// =====================
// Visual: Hohmann-like between circular proxies (use a proxy radius = a_km)
// ΔV: Hohmann between r1 and r2 + plane-change at apogee (simple baseline)
function makeTransferOrbit(chaserOrbit, targetOrbit) {
  const r1 = chaserOrbit.a_km; // proxy
  const r2 = targetOrbit.a_km;

  const a_t = 0.5 * (r1 + r2);
  const e_t = Math.abs(r2 - r1) / (r2 + r1);

  // Visual plane: blend when large inclination mismatch (helps readability)
  const di = Math.abs(chaserOrbit.i_deg - targetOrbit.i_deg);
  const w = di < 5 ? 0.0 : 0.5;
  const i_deg = (1 - w) * chaserOrbit.i_deg + w * targetOrbit.i_deg;
  const raan_deg = (1 - w) * chaserOrbit.raan_deg + w * targetOrbit.raan_deg;

  return { a_km: a_t, e: e_t, i_deg, raan_deg, argp_deg: 0.0 };
}

function circularSpeed(rKm) {
  return Math.sqrt(MU_EARTH / rKm); // km/s
}

function hohmannDV(r1, r2) {
  // km/s
  const v1 = circularSpeed(r1);
  const v2 = circularSpeed(r2);
  const a_t = 0.5 * (r1 + r2);
  const v_per = Math.sqrt(MU_EARTH * (2 / r1 - 1 / a_t));
  const v_apo = Math.sqrt(MU_EARTH * (2 / r2 - 1 / a_t));
  const dv1 = Math.abs(v_per - v1);
  const dv2 = Math.abs(v2 - v_apo);
  return { dv1, dv2, a_t, v_per, v_apo, v1, v2 };
}

function planeChangeDV(v, deltaI_rad) {
  // simple instantaneous plane change at speed v
  return 2 * v * Math.sin(Math.abs(deltaI_rad) / 2);
}

function planManeuvers() {
  const ch = state.chaserOrbit;
  const tg = state.targetOrbit;

  // Compute transfer orbit (visual)
  const transfer = makeTransferOrbit(ch, tg);
  state.transferOrbit = transfer;

  // Draw transfer (magenta dashed) + force on top
  if (transferLine) scene.remove(transferLine);
  transferLine = buildOrbitLine(transfer, transferMat, 720);
  transferLine.computeLineDistances();     // REQUIRED for dashed lines
  transferLine.renderOrder = 10;           // draw after Earth
  transferLine.visible = true;
  scene.add(transferLine);

  // Burns at perigee/apogee of transfer
  burn1.position.copy(positionOnOrbit(transfer, 0));
  burn2.position.copy(positionOnOrbit(transfer, Math.PI));
  burn1.visible = true;
  burn2.visible = true;

  // ---- Baseline ΔV output ----
  const r1 = ch.a_km;
  const r2 = tg.a_km;

  const { dv1, dv2, v_apo } = hohmannDV(r1, r2);

  // Plane change baseline: at apogee for cheaper plane change (not always optimal, but good baseline)
  const dI = deg2rad(tg.i_deg - ch.i_deg);
  const dvPlane = planeChangeDV(v_apo, dI);

  const total = dv1 + dv2 + dvPlane;

  ui.dvOut.textContent =
`Maneuver Plan (baseline): Chaser → Target
Assumptions: circular proxies, Hohmann + plane change at transfer apogee

ΔV1 (perigee burn): ${dv1.toFixed(3)} km/s
ΔV2 (apogee burn):  ${dv2.toFixed(3)} km/s
ΔV plane-change:    ${dvPlane.toFixed(3)} km/s
------------------------------
Total ΔV:           ${total.toFixed(3)} km/s

Notes:
- This is a baseline estimator (good for intuition).
- True rendezvous requires phasing + timing + targeting, not just ΔV.`;

  state.planned = true;
}

ui.btnPlan.addEventListener("click", () => {
  console.log("Plan Maneuvers clicked ✅");
  planManeuvers();
});

// =====================
// Start / Pause / Reset
// =====================
ui.btnStart.addEventListener("click", () => state.running = true);
ui.btnPause.addEventListener("click", () => state.running = false);

ui.btnReset.addEventListener("click", () => {
  state.running = false;
  applyOrbits(); // resets orbit lines + phases + hides transfer
  state.running = true;
});

// =====================
// Animation / propagation
// =====================
// Very simple: animate true anomaly at a constant rate for visualization.
// If you want physically correct mean motion, we can upgrade to n = sqrt(mu/a^3) and Kepler solve.
let lastT = performance.now();

function updateSun(betaDeg) {
  const b = deg2rad(betaDeg);
  // sun direction tilted "up" in +Z by beta (visual)
  const dir = new THREE.Vector3(Math.cos(b), 0, Math.sin(b)).normalize();
  sunArrow.setDirection(dir);
}

function animate(t) {
  requestAnimationFrame(animate);

  const dt = Math.min((t - lastT) / 1000, 0.05);
  lastT = t;

  controls.update();

  // Always render; only advance time when running
  const timeScale = safeNum(ui.spd.value, 30);

  if (state.running) {
    // visual angular speed (not physically exact)
    const omega = dt * (0.2 + timeScale * 0.01);

    state.cNu = (state.cNu + omega) % (Math.PI * 2);
    state.tNu = (state.tNu + omega * 0.95) % (Math.PI * 2); // slightly different rate for interest
  }

  // Update craft positions
  chaserCraft.position.copy(positionOnOrbit(state.chaserOrbit, state.cNu));
  targetCraft.position.copy(positionOnOrbit(state.targetOrbit, state.tNu));

  // Update sun direction
  updateSun(safeNum(ui.beta.value, 0));

  renderer.render(scene, camera);
}

applyOrbits();
animate(performance.now());

addEventListener("resize", () => {
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
});
