import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

/** =========================
 *  Constants (Earth-centered)
 *  ========================= */
const R_EARTH_KM = 6378.137;         // km
const MU_EARTH   = 398600.4418;      // km^3/s^2
const TWO_PI = Math.PI * 2;

const deg2rad = (d) => d * Math.PI / 180;
const rad2deg = (r) => r * 180 / Math.PI;

// Scene scaling: 1 Earth radius = 1 scene unit
const kmToScene = (km) => km / R_EARTH_KM;

/** =========================
 *  DOM / UI
 *  ========================= */
const $ = (id) => document.getElementById(id);

const ui = {
  btnStart: $("btnStart"),
  btnPause: $("btnPause"),
  btnReset: $("btnReset"),
  speed: $("speed"),
  spdVal: $("spdVal"),

  targetPreset: $("targetPreset"),
  t_rp: $("t_rp"), t_ra: $("t_ra"), t_i: $("t_i"), t_O: $("t_O"), t_w: $("t_w"), t_nu: $("t_nu"),

  chaserPreset: $("chaserPreset"),
  c_rp: $("c_rp"), c_ra: $("c_ra"), c_i: $("c_i"), c_O: $("c_O"), c_w: $("c_w"), c_nu: $("c_nu"),

  beta: $("beta"),
  betaVal: $("betaVal"),

  btnApply: $("btnApply"),
  btnPlan: $("btnPlan"),
  dvOut: $("dvOut"),
};


const planBtn = document.getElementById("planBtn"); // <- adjust id if different
console.log("planBtn:", planBtn);

planBtn?.addEventListener("click", () => {
  console.log("Plan Maneuvers clicked ✅");
});

/** =========================
 *  Three.js Scene
 *  ========================= */
const scene = new THREE.Scene();
scene.fog = new THREE.Fog(0x070b12, 3, 18);

const camera = new THREE.PerspectiveCamera(55, innerWidth / innerHeight, 0.01, 300);
camera.position.set(3.2, 2.2, 3.2);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(innerWidth, innerHeight);
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
document.body.appendChild(renderer.domElement);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.target.set(0, 0, 0);

scene.add(new THREE.AmbientLight(0xffffff, 0.45));
const key = new THREE.DirectionalLight(0xffffff, 1.0);
key.position.set(5, 8, 3);
scene.add(key);

// Planet
scene.add(new THREE.Mesh(
  new THREE.SphereGeometry(1, 64, 48),
  new THREE.MeshStandardMaterial({ color: 0x2d6cdf, roughness: 0.85, metalness: 0.0 })
));

// Equatorial plane visuals (strong)
const planeRadiusScene = kmToScene(R_EARTH_KM + 35786) * 1.15;
const eqDisk = new THREE.Mesh(
  new THREE.CircleGeometry(planeRadiusScene, 128),
  new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.06, side: THREE.DoubleSide })
);
scene.add(eqDisk);

const polarGrid = new THREE.PolarGridHelper(planeRadiusScene, 24, 12, 128, 0xffffff, 0xffffff);
polarGrid.material.transparent = true;
polarGrid.material.opacity = 0.18;
scene.add(polarGrid);

const axes = new THREE.AxesHelper(2.2);
axes.material.transparent = true;
axes.material.opacity = 0.7;
scene.add(axes);

// Sun direction arrow (beta visualization)
const sunDir = new THREE.Vector3(1, 0, 0);
const sunArrow = new THREE.ArrowHelper(sunDir, new THREE.Vector3(0,0,0), planeRadiusScene * 0.8, 0xffffff);
sunArrow.cone.material.transparent = true;
sunArrow.cone.material.opacity = 0.6;
sunArrow.line.material.transparent = true;
sunArrow.line.material.opacity = 0.6;
scene.add(sunArrow);

/** =========================
 *  Orbit utilities
 *  ========================= */

// Build orbit from perigee/apogee altitudes
function orbitFromPeriApoAlt({ rpAlt_km, raAlt_km, i_deg, raan_deg, argp_deg, nu0_deg }) {
  const rp = R_EARTH_KM + rpAlt_km;   // km radius
  const ra = R_EARTH_KM + raAlt_km;

  const a = 0.5 * (rp + ra);
  const e = (ra - rp) / (ra + rp);   // works for ra>=rp
  return {
    a_km: a,
    e,
    i_deg,
    raan_deg,
    argp_deg,
    // phase initialization (true anomaly)
    nu0_deg
  };
}

// Rotation: r_eci = Rz(Ω) * Rx(i) * Rz(ω) * r_pf
function rotateOrbitPoint(x, y, iRad, raanRad, argpRad) {
  // Rz(argp)
  const c3 = Math.cos(argpRad), s3 = Math.sin(argpRad);
  const x1 = c3 * x - s3 * y;
  const y1 = s3 * x + c3 * y;
  const z1 = 0;

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

function orbitalPeriod_s(a_km) {
  return TWO_PI * Math.sqrt((a_km**3) / MU_EARTH);
}

function meanMotion_rad_s(a_km) {
  return Math.sqrt(MU_EARTH / (a_km**3));
}

// Solve Kepler: M = E - e sin E for E
function solveKeplerE(M, e) {
  // normalize M to [-pi, pi] for stability
  M = ((M + Math.PI) % TWO_PI) - Math.PI;

  let E = e < 0.8 ? M : Math.PI; // initial guess
  for (let k = 0; k < 20; k++) {
    const f = E - e * Math.sin(E) - M;
    const fp = 1 - e * Math.cos(E);
    const dE = -f / fp;
    E += dE;
    if (Math.abs(dE) < 1e-10) break;
  }
  return E;
}

// Convert E to true anomaly nu
function E_to_nu(E, e) {
  const c = Math.cos(E);
  const s = Math.sin(E);
  const fac = Math.sqrt((1 + e) / (1 - e));
  return 2 * Math.atan2(fac * s/2, (1 + c)/2); // stable-ish form
}

// Position from time using mean anomaly propagation
function positionAtTime(orbit, t_s) {
  const { a_km, e } = orbit;

  // Initialize from nu0 -> E0 -> M0
  const nu0 = deg2rad(orbit.nu0_deg);
  const E0 = 2 * Math.atan2(Math.tan(nu0/2) * Math.sqrt(1 - e), Math.sqrt(1 + e));
  const M0 = E0 - e * Math.sin(E0);

  const n = meanMotion_rad_s(a_km);
  const M = M0 + n * t_s;

  const E = solveKeplerE(M, e);

  // true anomaly
  const nu = 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E/2), Math.sqrt(1 - e) * Math.cos(E/2));

  // radius
  const r_km = a_km * (1 - e * Math.cos(E));

  // perifocal
  const x_pf = kmToScene(r_km) * Math.cos(nu);
  const y_pf = kmToScene(r_km) * Math.sin(nu);

  // rotate
  const iRad = deg2rad(orbit.i_deg);
  const raanRad = deg2rad(orbit.raan_deg);
  const argpRad = deg2rad(orbit.argp_deg);

  return rotateOrbitPoint(x_pf, y_pf, iRad, raanRad, argpRad);
}

// Orbit polyline
function buildOrbitLine(orbit, material, segments = 512) {
  const { a_km, e } = orbit;
  const iRad = deg2rad(orbit.i_deg);
  const raanRad = deg2rad(orbit.raan_deg);
  const argpRad = deg2rad(orbit.argp_deg);

  const p = a_km * (1 - e * e);

  const pts = [];
  for (let k = 0; k <= segments; k++) {
    const nu = (k / segments) * TWO_PI;
    const r_km = p / (1 + e * Math.cos(nu));
    const x_pf = kmToScene(r_km) * Math.cos(nu);
    const y_pf = kmToScene(r_km) * Math.sin(nu);
    pts.push(rotateOrbitPoint(x_pf, y_pf, iRad, raanRad, argpRad));
  }

  return new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), material);
}

/** =========================
 *  Presets
 *  ========================= */
function applyTargetPreset(name) {
  if (name === "ISS") {
    ui.t_rp.value = 400; ui.t_ra.value = 400;
    ui.t_i.value = 51.6; ui.t_O.value = 0; ui.t_w.value = 0; ui.t_nu.value = 0;
  } else if (name === "GEO") {
    ui.t_rp.value = 35786; ui.t_ra.value = 35786;
    ui.t_i.value = 0; ui.t_O.value = 0; ui.t_w.value = 0; ui.t_nu.value = 0;
  } else if (name === "SSO") {
    ui.t_rp.value = 600; ui.t_ra.value = 600;
    ui.t_i.value = 98; ui.t_O.value = 40; ui.t_w.value = 0; ui.t_nu.value = 0;
  } else if (name === "MOLNIYA") {
    // Molniya typical: a~26600 km, e~0.74; approximate with rp/ra
    // Solve rp/ra from a,e: rp=a(1-e), ra=a(1+e)
    const a = 26600, e = 0.74;
    const rp = a * (1 - e), ra = a * (1 + e);
    ui.t_rp.value = Math.max(0, rp - R_EARTH_KM);
    ui.t_ra.value = Math.max(0, ra - R_EARTH_KM);
    ui.t_i.value = 63.4; ui.t_O.value = 20; ui.t_w.value = 270; ui.t_nu.value = 0;
  }
}

function applyChaserPreset(name) {
  if (name === "LEO200") {
    ui.c_rp.value = 200; ui.c_ra.value = 200;
    ui.c_i.value = 51.6; ui.c_O.value = 0; ui.c_w.value = 0; ui.c_nu.value = -40;
  } else if (name === "DRAGON_TO_ISS") {
    // Demonstration: lower circular orbit with same plane, behind target
    ui.c_rp.value = 200; ui.c_ra.value = 200;
    ui.c_i.value = 51.6; ui.c_O.value = 0; ui.c_w.value = 0; ui.c_nu.value = -60;
  } else if (name === "GTO") {
    // GTO-ish (visual): ~250 km x 35786 km, i ~ 28.5
    ui.c_rp.value = 250; ui.c_ra.value = 35786;
    ui.c_i.value = 28.5; ui.c_O.value = 0; ui.c_w.value = 0; ui.c_nu.value = 0;
  }
}

/** =========================
 *  Materials / Objects
 *  ========================= */
const matTarget = new THREE.LineBasicMaterial({ color: 0x3fdcff, transparent: true, opacity: 0.95 });
const matChaser = new THREE.LineBasicMaterial({ color: 0xffd54a, transparent: true, opacity: 0.95 });

let lineTarget = null;
let lineChaser = null;

const targetMarker = new THREE.Mesh(
  new THREE.SphereGeometry(0.03, 16, 12),
  new THREE.MeshBasicMaterial({ color: 0x3fdcff })
);
const chaserMarker = new THREE.Mesh(
  new THREE.SphereGeometry(0.03, 16, 12),
  new THREE.MeshBasicMaterial({ color: 0xffd54a })
);
scene.add(targetMarker, chaserMarker);

/** =========================
 *  State / simulation loop
 *  ========================= */
let running = false;
let tSim_s = 0;

let targetOrbit = null;
let chaserOrbit = null;

function readOrbitsFromUI() {
  targetOrbit = orbitFromPeriApoAlt({
    rpAlt_km: parseFloat(ui.t_rp.value),
    raAlt_km: parseFloat(ui.t_ra.value),
    i_deg: parseFloat(ui.t_i.value),
    raan_deg: parseFloat(ui.t_O.value),
    argp_deg: parseFloat(ui.t_w.value),
    nu0_deg: parseFloat(ui.t_nu.value),
  });

  chaserOrbit = orbitFromPeriApoAlt({
    rpAlt_km: parseFloat(ui.c_rp.value),
    raAlt_km: parseFloat(ui.c_ra.value),
    i_deg: parseFloat(ui.c_i.value),
    raan_deg: parseFloat(ui.c_O.value),
    argp_deg: parseFloat(ui.c_w.value),
    nu0_deg: parseFloat(ui.c_nu.value),
  });
}

function rebuildLines() {
  if (lineTarget) scene.remove(lineTarget);
  if (lineChaser) scene.remove(lineChaser);

  lineTarget = buildOrbitLine(targetOrbit, matTarget);
  lineChaser = buildOrbitLine(chaserOrbit, matChaser);
  scene.add(lineTarget, lineChaser);
}

function updateSunFromBeta() {
  const beta = deg2rad(parseFloat(ui.beta.value));
  ui.betaVal.textContent = `${ui.beta.value}`;

  // Sun direction: rotate out of equatorial plane by beta about +Y (simple visual convention)
  // (You can later tie this to RAAN/time of year.)
  sunDir.set(Math.cos(beta), 0, Math.sin(beta)).normalize();
  sunArrow.setDirection(sunDir);

  // Also use as a light direction (simple “day/night” cue)
  key.position.set(sunDir.x * 8, sunDir.y * 8 + 4, sunDir.z * 8);
}

ui.beta.addEventListener("input", updateSunFromBeta);

// Buttons
ui.btnStart.onclick = () => { running = true; };
ui.btnPause.onclick = () => { running = false; };
ui.btnReset.onclick = () => { running = false; tSim_s = 0; };

ui.speed.addEventListener("input", () => ui.spdVal.textContent = ui.speed.value);

ui.targetPreset.addEventListener("change", () => applyTargetPreset(ui.targetPreset.value));
ui.chaserPreset.addEventListener("change", () => applyChaserPreset(ui.chaserPreset.value));

ui.btnApply.onclick = () => {
  readOrbitsFromUI();
  rebuildLines();
};

window.addEventListener("resize", () => {
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
});

// init presets + draw
applyTargetPreset(ui.targetPreset.value);
applyChaserPreset(ui.chaserPreset.value);
readOrbitsFromUI();
rebuildLines();
updateSunFromBeta();

/** =========================
 *  ΔV Planner (baseline)
 *  ========================= */

// vis-viva
function speedAtRadius(a_km, r_km) {
  return Math.sqrt(MU_EARTH * (2 / r_km - 1 / a_km));
}

// Circular speed
function vCircular(r_km) {
  return Math.sqrt(MU_EARTH / r_km);
}

// Plane angle between orbital planes from (i, Ω)
function planeAngle_rad(i1, O1, i2, O2) {
  // plane normal unit vectors in inertial frame:
  // n = [sin i * sin Ω, -sin i * cos Ω, cos i]
  const si1 = Math.sin(i1), ci1 = Math.cos(i1);
  const si2 = Math.sin(i2), ci2 = Math.cos(i2);

  const n1 = new THREE.Vector3(si1 * Math.sin(O1), -si1 * Math.cos(O1), ci1).normalize();
  const n2 = new THREE.Vector3(si2 * Math.sin(O2), -si2 * Math.cos(O2), ci2).normalize();

  const dot = THREE.MathUtils.clamp(n1.dot(n2), -1, 1);
  return Math.acos(dot);
}

function planDv() {
  // Build radii in km at "circular proxy" (use a for circular or rp for perigee)
  const r1 = targetOrbit ? null : null;

  // We'll plan: chaser -> target
  // For now treat each as circular at its semi-major axis (good for LEO/GEO circular demos).
  const rCh = chaserOrbit.a_km;
  const rTg = targetOrbit.a_km;

  // Hohmann transfer ΔV (coplanar)
  const aTrans = 0.5 * (rCh + rTg);

  const vCh = vCircular(rCh);
  const vTg = vCircular(rTg);

  const vPeriTrans = speedAtRadius(aTrans, rCh);
  const vApoTrans  = speedAtRadius(aTrans, rTg);

  const dV1 = Math.abs(vPeriTrans - vCh);
  const dV2 = Math.abs(vTg - vApoTrans);

  // Plane change (approx) — do at the slower point (apogee if raising, perigee if lowering)
  const i1 = deg2rad(chaserOrbit.i_deg);
  const O1 = deg2rad(chaserOrbit.raan_deg);
  const i2 = deg2rad(targetOrbit.i_deg);
  const O2 = deg2rad(targetOrbit.raan_deg);

  const dTheta = planeAngle_rad(i1, O1, i2, O2); // includes RAAN difference effect

  const vSlow = (rTg >= rCh) ? vApoTrans : vPeriTrans; // slower end of transfer
  const dVplane = 2 * vSlow * Math.sin(dTheta / 2);

  const total = dV1 + dV2 + dVplane;

  // Output text
  const lines = [];
  lines.push(`Maneuver Plan (baseline): Chaser → Target`);
  lines.push(`Assumptions: circular proxies, Hohmann + plane change at slow point`);
  lines.push(``);
  lines.push(`Burn 1 (transfer inject):  ΔV1 = ${dV1.toFixed(3)} km/s`);
  lines.push(`Burn 2 (circularize):     ΔV2 = ${dV2.toFixed(3)} km/s`);
  lines.push(`Plane change (approx):     ΔV = ${dVplane.toFixed(3)} km/s`);
  lines.push(`----------------------------------------`);
  lines.push(`TOTAL ΔV:                  ${total.toFixed(3)} km/s`);
  lines.push(``);
  lines.push(`Notes:`);
  lines.push(`- If chaser & target planes match, plane change term drops near zero.`);
  lines.push(`- Next upgrade: add phasing (rendezvous) ΔV when same altitude but different phase.`);
  ui.dvOut.textContent = lines.join("\n");
}

ui.btnPlan.onclick = planDv;

/** =========================
 *  Animation loop
 *  ========================= */
let lastT = performance.now();

function loop(t) {
  requestAnimationFrame(loop);
  const dtReal = Math.min((t - lastT) / 1000, 0.05);
  lastT = t;

  controls.update();

  // sim speed multiplier
  const speedMult = parseFloat(ui.speed.value);

  if (running) tSim_s += dtReal * speedMult;

  if (targetOrbit && chaserOrbit) {
    targetMarker.position.copy(positionAtTime(targetOrbit, tSim_s));
    chaserMarker.position.copy(positionAtTime(chaserOrbit, tSim_s));
  }

  renderer.render(scene, camera);
}
loop(performance.now());
