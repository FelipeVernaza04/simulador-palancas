/* Modelo sagital 2D cuasi-estático de palancas (puerto de palancas/*.py).
 * x adelante, y arriba, origen en el suelo bajo el medio pie. Momentos extensores positivos (N·m), carga total.
 * Unidades: m, kg, rad. Los números de de Leva se transcribieron de un PDF: verificar contra la tabla original.
 */
(function (root) {
  "use strict";
  const G = 9.81;
  const LUMBAR_MU = 0.10;
  const KNEE_FRONT = 0.05;  // centro articular de la rodilla -> cara anterior de la rótula, con partes blandas (m) (aproximación)
  const BAR_RADIUS = 0.014; // barra olímpica de 28 mm
  const KNEE_CLEAR = KNEE_FRONT + BAR_RADIUS; // distancia mínima centro de barra - centro de rodilla en el plano sagital

  class Infeasible extends Error {}

  // ---------------------------------------------------------------- antropometría
  const DE_LEVA = {
    M: { thigh: [0.1416, 0.4095], shank: [0.0433, 0.4459], foot: [0.0137, 0.4415], head: 0.0694, trunk: [0.4346, 0.4486],
         upper_arm: [0.0271, 0.5772, 0.2817], forearm: [0.0162, 0.4574, 0.2689], hand: [0.0061, 0.7900, 0.0862] },
    F: { thigh: [0.1478, 0.3612], shank: [0.0481, 0.4416], foot: [0.0129, 0.4014], head: 0.0668, trunk: [0.4257, 0.4151],
         upper_arm: [0.0255, 0.5754, 0.2751], forearm: [0.0138, 0.4559, 0.2643], hand: [0.0056, 0.7474, 0.0780] },
  };
  const SHOULDER_BELOW_ACROMION = { M: 0.0345, F: 0.0337 };
  const DEFAULT_HIP_AP_OFFSET = 0.058;
  const HEAD_COM_ABOVE_SHOULDER_FRAC_H = 0.12 / 1.75;

  function armMassAndComFraction(sex) {
    const p = DE_LEVA[sex];
    const [mUa, cUa, LUa] = p.upper_arm, [mFa, cFa, LFa] = p.forearm, [mH, cH, LH] = p.hand;
    const pos = [cUa * LUa, LUa + cFa * LFa, LUa + LFa + cH * LH];
    const masses = [mUa, mFa, mH];
    let com = 0;
    for (let i = 0; i < 3; i++) com += masses[i] * pos[i];
    com /= masses[0] + masses[1] + masses[2];
    const Lgrip = LUa + LFa + 0.5 * LH;
    return [2 * (masses[0] + masses[1] + masses[2]), com / Lgrip];
  }

  /** raw: {sex, mass, H, Lt_raw, Ls_raw, ha, Lf, heel_to_malleolus, Ltr_raw, La_raw, D_asis} en m / kg */
  function fromRaw(r, hipVertical = 0, hipAp = DEFAULT_HIP_AP_OFFSET) {
    const sh = SHOULDER_BELOW_ACROMION[r.sex];
    return {
      sex: r.sex, mass: r.mass, H: r.H,
      Lt: r.Lt_raw - hipVertical, Ls: r.Ls_raw, ha: r.ha, Lf: r.Lf,
      dam: r.Lf / 2 - r.heel_to_malleolus,
      Ltr: r.Ltr_raw - sh + hipVertical, La: r.La_raw - sh,
      b: 0.66 * r.D_asis + 0.0146,
      hip_ap_offset: hipAp,
      shoulder_w: r.D_acro || 0,
    };
  }

  function closureChecks(r) {
    const out = {};
    if (r.h_troc != null) { const v = r.h_troc - (r.ha + r.Ls_raw + r.Lt_raw); out.legs_cm = v * 100; out.legs_ok = Math.abs(v) <= 0.02; }
    if (r.h_acro != null && r.h_troc != null) { const v = r.h_acro - (r.h_troc + r.Ltr_raw); out.trunk_cm = v * 100; out.trunk_ok = Math.abs(v) <= 0.03; }
    return out;
  }

  function massPoints(a) {
    const p = DE_LEVA[a.sex];
    const [armsFrac, armCom] = armMassAndComFraction(a.sex);
    return {
      m_thigh: p.thigh[0] * a.mass, f_thigh: p.thigh[1],
      m_shank: p.shank[0] * a.mass, f_shank: p.shank[1],
      m_foot: p.foot[0] * a.mass, f_foot: p.foot[1],
      m_head: p.head * a.mass,
      m_trunk: p.trunk[0] * a.mass, f_trunk_from_hip: 1 - p.trunk[1],
      m_arms: armsFrac * a.mass, f_arm_com: armCom,
      head_above_shoulder: HEAD_COM_ABOVE_SHOULDER_FRAC_H * a.H,
    };
  }

  // ---------------------------------------------------------------- técnicas (calibradas con literatura; ver README)
  const BASE = { calibrated: false, source: "", ranges: {}, at_bound: [], lam: 1.0, d: 0.06, shank_top_deg: 3.0, shank_max_deg: 35.0,
    depth_deg: 5.0, bar_x_mode: "midfoot", trunk_top_deg: 0.0, bar_x_offset: 0.0, f_h0: 0.55, eps0: 0.05, eps_top: -0.03, p: 1.0, plate_radius: 0.225, stance: "frac",
    toe_out_deg: 0, knee_out_deg: 0, grip_inside: false };
  function T(o) { return Object.assign({}, BASE, o); }
  const TECHNIQUES = {
    sentadilla_barra_alta: T({ name: "Sentadilla barra alta", exercise: "squat", w_frac: 48.3 / 180.3, toe_out_deg: 0, knee_out_deg: 12, stance: "shoulder", calibrated: true,
      lam: 1.0, d: 0.06, shank_max_deg: 32.8, depth_deg: -1.9, bar_x_offset: 0.040 }),
    sentadilla_barra_baja: T({ name: "Sentadilla barra baja", exercise: "squat", w_frac: 89.6 / 180.3, toe_out_deg: 0, knee_out_deg: 12, stance: "shoulder", calibrated: true,
      lam: 0.93, d: 0.07, shank_max_deg: 20.8, depth_deg: 2.7, bar_x_offset: -0.021 }),
    sentadilla_frontal: T({ name: "Sentadilla frontal", exercise: "squat", w_frac: 0.30, toe_out_deg: 0, knee_out_deg: 12, stance: "shoulder", calibrated: false,
      lam: 1.0, d: -0.05, shank_max_deg: 40, depth_deg: 5, bar_x_offset: 0.0 }),
    peso_muerto_convencional: T({ name: "Peso muerto convencional", exercise: "deadlift", w_frac: 37 / 167, toe_out_deg: 15, knee_out_deg: 15, stance: "shoulder", calibrated: true,
      f_h0: 0.8424, eps0: 0.12, eps_top: 0.0177, p: 1.5633, bar_x_offset: 0 }),
    peso_muerto_sumo: T({ name: "Peso muerto sumo", exercise: "deadlift", w_frac: 73 / 168, toe_out_deg: 45, knee_out_deg: 25, grip_inside: true, calibrated: true,
      f_h0: 0.7762, eps0: 0.102, eps_top: -0.0153, p: 1.2119, bar_x_offset: 0 }),
  };

  // ---------------------------------------------------------------- modelo
  function csqrt(v, tol = 1e-9) {
    if (v < -tol) throw new Infeasible("radicando negativo");
    return Math.sqrt(Math.max(v, 0));
  }
  function stanceC(a, t) {
    const w = (t.stance === "shoulder" && a.shoulder_w > 0) ? a.shoulder_w : t.w_frac * a.H;
    return Math.max((w - a.b) / 2, 0);
  }

  /** Rodilla de la pierna derecha para un tobillo y una cadera 3D, puntos [x adelante, y arriba, z lateral hacia afuera].
   * Las longitudes de tibia y muslo fijan una circunferencia de posiciones posibles alrededor del eje tobillo-cadera. La rodilla se
   * abre en el plano que contiene ese eje y la dirección horizontal (cos th, 0, sin th), con th = outDeg hacia afuera: la rodilla
   * sigue la dirección del pie. Con th = 0 es el plano que contiene al eje x (el modelo plano original). Espejo de model._knee_3d. */
  function knee3d(a, A, H, outDeg) {
    const d = [H[0] - A[0], H[1] - A[1], H[2] - A[2]];
    const D = Math.hypot(d[0], d[1], d[2]);
    if (D > a.Ls + a.Lt + 1e-9 || D < Math.abs(a.Ls - a.Lt) - 1e-9) throw new Infeasible("cadera fuera del alcance de la pierna");
    const w = [d[0] / D, d[1] / D, d[2] / D];
    const aa = (a.Ls * a.Ls - a.Lt * a.Lt + D * D) / (2 * D);
    const h = csqrt(a.Ls * a.Ls - aa * aa);
    const th = outDeg * Math.PI / 180, u = [Math.cos(th), 0, Math.sin(th)];
    const pr = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
    let n = [u[0] - pr * w[0], u[1] - pr * w[1], u[2] - pr * w[2]];
    const nn = Math.hypot(n[0], n[1], n[2]);
    n = nn > 1e-9 ? [n[0] / nn, n[1] / nn, n[2] / nn] : [1, 0, 0];
    return [A[0] + aa * w[0] + h * n[0], A[1] + aa * w[1] + h * n[1], A[2] + aa * w[2] + h * n[2]];
  }

  function momentsOf(a, mb, psi, xa, xk, xh, shoulder, phi, barX, armsX, toeOut) {
    const mp = massPoints(a);
    const sph = Math.sin(phi);
    const xTrunk = xh + mp.f_trunk_from_hip * a.Ltr * sph;
    const xHead = shoulder[0] + mp.head_above_shoulder * sph;
    const xThigh = xh + mp.f_thigh * (xk - xh);
    const xShank = xk + mp.f_shank * (xa - xk);
    const xL = xh + LUMBAR_MU * a.Ltr * sph;
    const xhc = xh + a.hip_ap_offset;

    const above = [[mp.m_trunk, xTrunk], [mp.m_head, xHead], [mp.m_arms, armsX]];
    const bar = [[mb, barX]];
    const thighs = [[2 * mp.m_thigh, xThigh]];
    const shanks = [[2 * mp.m_shank, xShank]];
    const cpsi = Math.cos(psi);
    const fwd = (items, xj) => G * items.reduce((s, [m, x]) => s + m * (x - xj), 0);
    const bwd = (items, xj) => G * items.reduce((s, [m, x]) => s + m * (xj - x), 0);

    const barI = { hip: bar, knee: bar, ankle: bar, lumbar: bar };
    const bodyI = { hip: above, knee: above.concat(thighs), ankle: above.concat(thighs, shanks), lumbar: above };
    const part = (I) => ({ cadera: fwd(I.hip, xhc) * cpsi, rodilla: bwd(I.knee, xk) * cpsi, tobillo: fwd(I.ankle, xa) * cpsi, lumbar: fwd(I.lumbar, xL) });
    const pb = part(barI), pd = part(bodyI);
    const total = {};
    for (const k in pb) total[k] = pb[k] + pd[k];
    const arms = { cadera: barX - xhc, rodilla: xk - barX, tobillo: barX - xa, lumbar: barX - xL };

    const xFoot = -a.Lf / 2 + mp.f_foot * a.Lf;
    const allm = above.concat(bar, thighs, shanks, [[2 * mp.m_foot, xFoot]]);
    let sm = 0, sx = 0;
    for (const [m, x] of allm) { sm += m; sx += m * x; }
    const cop = sx / sm;
    const half = a.Lf * Math.cos(toeOut * Math.PI / 180) / 2;
    return { total, bar: pb, body: pd, arms, cop, copIn: cop >= -half && cop <= half, half };
  }

  function pack(a, t, s, mb, ankle, knee, hip, shoulder, barXY, psi, phi, barX, armsX, lateral) {
    const m = momentsOf(a, mb, psi, ankle[0], knee[0], hip[0], shoulder, phi, barX, armsX, t.toe_out_deg);
    const legs = { tobillo: Math.abs(m.total.tobillo), rodilla: Math.abs(m.total.rodilla), cadera: Math.abs(m.total.cadera) };
    const tot = legs.tobillo + legs.rodilla + legs.cadera || 1;
    return {
      s, joints: { tobillo: ankle, rodilla: knee, cadera: hip, hombro: shoulder }, lateral, bar: barXY,
      psi_deg: psi * 180 / Math.PI, trunk_deg: phi * 180 / Math.PI,
      moments: m.total, bar_moments: m.bar, body_moments: m.body, arms: m.arms,
      shares: { tobillo: legs.tobillo / tot, rodilla: legs.rodilla / tot, cadera: legs.cadera / tot },
      cop_x: m.cop, cop_in_base: m.copIn, base_half: m.half,
    };
  }

  function solveSquat(a, t, load, s) {
    // Pierna parametrizada por el ángulo del muslo respecto a la vertical (beta), lineal en s (ver palancas/model.py).
    const c = stanceC(a, t), xa = -a.dam, ha = a.ha;
    const ph0 = t.shank_top_deg * Math.PI / 180, phm = t.shank_max_deg * Math.PI / 180, al = t.depth_deg * Math.PI / 180;
    const beta = ph0 + s * ((Math.PI / 2 + al) - ph0);
    const phs = ph0 + s * (phm - ph0);
    const uK = xa + a.Ls * Math.sin(phs), vK = a.Ls * Math.cos(phs);
    const vH = vK + a.Lt * Math.cos(beta);
    const uH = uK - a.Lt * Math.sin(beta);
    const yp = csqrt(vH * vH - c * c);
    const psi = Math.atan2(c, yp);
    const A = t.lam * a.Ltr, B = t.d;
    const xb = t.bar_x_offset;
    const C = xb - uH, R = Math.hypot(A, B);
    if (Math.abs(C) > R) throw new Infeasible("la barra no puede quedar sobre la vertical elegida con este tronco");
    const phi = Math.atan2(B, A) + Math.asin(C / R);
    const hip = [uH, ha + yp];
    const zH = a.b / 2, zA = zH + c;
    // la cadera y el tobillo vienen del plano de la pierna; la rodilla se abre según t.knee_out_deg (ver knee3d)
    const K = knee3d(a, [xa, ha, zA], [uH, hip[1], zH], t.knee_out_deg);
    const knee = [K[0], K[1]];
    const shoulder = [uH + a.Ltr * Math.sin(phi), hip[1] + a.Ltr * Math.cos(phi)];
    const barXY = [uH + A * Math.sin(phi) - B * Math.cos(phi), hip[1] + A * Math.cos(phi) + B * Math.sin(phi)];
    return pack(a, t, s, load, [xa, ha], knee, hip, shoulder, barXY, psi, phi, barXY[0], shoulder[0], { tobillo: zA, rodilla: K[2], cadera: zH });
  }

  function solveDeadlift(a, t, load, s) {
    const c = stanceC(a, t), xa = -a.dam, ha = a.ha, xb = t.bar_x_offset;
    const uTop = xb + t.eps_top;
    const vTop = csqrt((a.Ls + a.Lt) ** 2 - (uTop - xa) ** 2);
    const ypTop = csqrt(vTop * vTop - c * c);
    const yp0 = t.f_h0 * ypTop;
    const yp = yp0 + (ypTop - yp0) * Math.pow(s, t.p);
    const eps = t.eps0 + (t.eps_top - t.eps0) * s;
    const yHTop = ha + ypTop;
    const yBTop = yHTop + a.Ltr - csqrt(a.La * a.La - t.eps_top * t.eps_top);
    const yB = t.plate_radius + s * (yBTop - t.plate_radius);
    const yH = ha + yp, psi = Math.atan2(c, yp), zH = a.b / 2, zA = zH + c;
    const xSh = xb + eps; // el hombro queda donde lo pone la ficha de técnica
    // La barra no puede atravesar la rodilla: si su centro pasa a menos de KNEE_CLEAR del centro de la rodilla, la barra se adelanta
    // lo justo respecto del cuerpo (hombro, cadera y rodilla no se mueven: los brazos se inclinan un poco hacia adelante). Con la
    // rodilla abierta en la dirección del pie casi no actúa. Se itera porque la altura del hombro depende de la inclinación del brazo.
    let xbe = xb, ySh = 0, dy = 0, xH = 0, K = null;
    for (let it = 0; it < 12; it++) {
      const e = xSh - xbe;
      ySh = yB + csqrt(a.La * a.La - e * e);
      dy = ySh - yH;
      if (dy < 0 || dy > a.Ltr + 1e-9) throw new Infeasible("el hombro no es alcanzable desde esa altura de cadera");
      xH = xSh - csqrt(a.Ltr * a.Ltr - dy * dy);
      K = knee3d(a, [xa, ha, zA], [xH, yH, zH], t.knee_out_deg);
      const dyk = Math.abs(yB - K[1]);
      const need = dyk < KNEE_CLEAR ? K[0] + Math.sqrt(KNEE_CLEAR * KNEE_CLEAR - dyk * dyk) : -Infinity;
      if (need <= xbe + 1e-7) break;
      xbe = need;
    }
    const phi = Math.atan2(xSh - xH, dy);
    const hip = [xH, yH], knee = [K[0], K[1]], shoulder = [xSh, ySh];
    const fArm = massPoints(a).f_arm_com;
    const armsX = xSh + fArm * (xbe - xSh);
    return pack(a, t, s, load, [xa, ha], knee, hip, shoulder, [xbe, yB], psi, phi, xbe, armsX, { tobillo: zA, rodilla: K[2], cadera: zH });
  }

  function solve(a, t, load, s) { return (t.exercise === "squat" ? solveSquat : solveDeadlift)(a, t, load, s); }

  // ---------------------------------------------------------------- incertidumbre (Monte Carlo)
  const TOL = { Lt: 0.010, Ls: 0.010, ha: 0.010, Ltr: 0.010, La: 0.010, D_asis: 0.010, Lf: 0.005, heel_to_malleolus: 0.005,
    hip_vertical: 0.015, shoulder: 0.005, hip_ap: 0.04, D_acro: 0.010 };
  const JOINTS = ["tobillo", "rodilla", "cadera", "lumbar"];

  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function makeNormal(rng) {
    return function () {
      let u = 0, v = 0;
      while (u === 0) u = rng();
      v = rng();
      return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    };
  }

  /** n antropometrías perturbadas (SD = tolerancia/2). raw con las medidas del atleta; tol por defecto TOL. */
  function perturbedSet(raw, n, seed = 1, tol = TOL) {
    const norm = makeNormal(mulberry32(seed));
    const keys = Object.keys(tol);
    const out = [];
    for (let i = 0; i < n; i++) {
      const d = {};
      for (const k of keys) d[k] = norm() * tol[k] / 2;
      const r = Object.assign({}, raw, {
        Lt_raw: raw.Lt_raw + d.Lt, Ls_raw: raw.Ls_raw + d.Ls, ha: raw.ha + d.ha,
        Ltr_raw: raw.Ltr_raw + d.Ltr + d.shoulder, La_raw: raw.La_raw + d.La,
        D_asis: raw.D_asis + d.D_asis, Lf: raw.Lf + d.Lf, heel_to_malleolus: raw.heel_to_malleolus + d.heel_to_malleolus,
        D_acro: raw.D_acro > 0 ? raw.D_acro + d.D_acro : 0,
      });
      out.push(fromRaw(r, d.hip_vertical, DEFAULT_HIP_AP_OFFSET + d.hip_ap));
    }
    return out;
  }

  function percentile(sorted, q) {
    if (!sorted.length) return NaN;
    const pos = (sorted.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
    return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
  }

  /** Para un conjunto perturbado y una posición s: percentiles 2,5/50/97,5 de cada momento y número de sorteos válidos. */
  function band(set, t, load, s) {
    const vals = { tobillo: [], rodilla: [], cadera: [], lumbar: [] };
    for (const a of set) {
      try {
        const m = solve(a, t, load, s).moments;
        for (const k of JOINTS) vals[k].push(m[k]);
      } catch (e) { if (!(e instanceof Infeasible)) throw e; }
    }
    const out = { n: vals.cadera.length };
    for (const k of JOINTS) {
      const v = vals[k].slice().sort((x, y) => x - y);
      out[k] = [percentile(v, 0.025), percentile(v, 0.5), percentile(v, 0.975)];
    }
    return out;
  }

  /** Comparación pareada B−A con las mismas perturbaciones: mediana de la diferencia y P(B<A). */
  function pairedDiff(set, tA, tB, load, s) {
    const d = { tobillo: [], rodilla: [], cadera: [], lumbar: [] };
    for (const a of set) {
      try {
        const mA = solve(a, tA, load, s).moments, mB = solve(a, tB, load, s).moments;
        for (const k of JOINTS) d[k].push(mB[k] - mA[k]);
      } catch (e) { if (!(e instanceof Infeasible)) throw e; }
    }
    const out = {};
    for (const k of JOINTS) {
      const v = d[k].slice().sort((x, y) => x - y);
      out[k] = { median: percentile(v, 0.5), lo: percentile(v, 0.025), hi: percentile(v, 0.975), pLess: v.length ? v.filter(x => x < 0).length / v.length : NaN };
    }
    return out;
  }

  const api = { G, Infeasible, DE_LEVA, DEFAULT_HIP_AP_OFFSET, KNEE_CLEAR, TECHNIQUES, TOL, JOINTS, fromRaw, closureChecks, massPoints, solve, perturbedSet, band, pairedDiff };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Palancas = api;
})(typeof self !== "undefined" ? self : this);
