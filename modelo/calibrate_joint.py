"""Calibración conjunta: las cuatro técnicas comparten un único desfase antero-posterior del centro de cadera (hx).

Motivo: en el ajuste independiente las cuatro técnicas dieron un brazo de momento de cadera 3-7 cm mayor que el publicado, lo que apunta a
una diferencia de definición del centro de cadera (modelo: vértice basado en el trocánter; literatura: centro articular estimado por marcadores).
hx se estima junto con los parámetros de estilo, con un a priori débil N(0, 3 cm). Es una INFERENCIA a partir de los datos, no una medida directa.
"""
import json
from dataclasses import replace
import math
import numpy as np
from scipy.optimize import least_squares

from anthropometry import from_raw
from techniques import TECHNIQUES
from model import solve, Infeasible
from calibrate import (proxy_athlete, incl_3d, knee_flexion_inplane, knee_passing_s, DL, DL_LABELS, SQ, SQ_LABELS)

HX_PRIOR_SD = 0.03
DL_BOUNDS = ([0.35, -0.02, -0.08, 0.6, -0.25], [0.90, 0.12, 0.02, 1.6, 0.10])
SQ_BOUNDS = ([10.0, -15.0, -0.15, 0.85, 0.03], [50.0, 25.0, 0.15, 1.05, 0.10])
HX_BOUNDS = ([-0.05], [0.10])
LAYOUT = [("peso_muerto_convencional", "dl", 5), ("peso_muerto_sumo", "dl", 5),
          ("sentadilla_barra_alta", "sq", 5), ("sentadilla_barra_baja", "sq", 5)]
DL_NAMES = ("f_h0", "eps0", "eps_top", "p", "bar_x_offset")
SQ_NAMES = ("shank_max_deg", "depth_deg", "bar_x_offset", "lam", "d")


def tech_from(key, kind, x):
    if kind == "dl":
        d = DL[key]
        return replace(TECHNIQUES[key], w_frac=d["w_frac"], toe_out_deg=d["toe"], f_h0=x[0], eps0=x[1], eps_top=x[2], p=x[3],
                       bar_x_offset=x[4], calibrated=True)
    return replace(TECHNIQUES[key], w_frac=SQ[key]["w_frac"], bar_x_mode="midfoot", shank_max_deg=x[0], depth_deg=x[1],
                   bar_x_offset=x[2], lam=x[3], d=x[4], calibrated=True)


def predict(key, kind, x, hx):
    if kind == "dl":
        a = from_raw(proxy_athlete(1.675, 76.0), hip_ap_offset=hx)
        t, load = tech_from(key, kind, x), 130.0
        try:
            r0 = solve(a, t, load, 0.0)
            rk = solve(a, t, load, knee_passing_s(a, t, load))
            r1 = solve(a, t, load, 1.0)
        except Infeasible:
            return None
        return np.array([r0.arms["cadera"] * 100, -r0.arms["rodilla"] * 100, r0.arms["tobillo"] * 100,
                         incl_3d(r0, "rodilla", "cadera", a.ha), incl_3d(rk, "tobillo", "rodilla", a.ha),
                         rk.arms["cadera"] * 100, -rk.arms["rodilla"] * 100, rk.arms["tobillo"] * 100,
                         (r1.bar[1] - r0.bar[1]) * 100])
    a = from_raw(proxy_athlete(1.803, 100.2), hip_ap_offset=hx)
    try:
        r = solve(a, tech_from(key, kind, x), 154.0, 1.0)
    except Infeasible:
        return None
    return np.array([r.arms["cadera"] * 100, -r.arms["rodilla"] * 100, r.arms["tobillo"] * 100,
                     knee_flexion_inplane(r, a), x[0], r.trunk_deg])


def predict_from_technique(key):
    """Evalúa la ficha vigente de TECHNIQUES contra los objetivos publicados (z-scores)."""
    from anthropometry import DEFAULT_HIP_AP_OFFSET
    t = TECHNIQUES[key]
    kind = "dl" if t.exercise == "deadlift" else "sq"
    x = ([t.f_h0, t.eps0, t.eps_top, t.p, t.bar_x_offset] if kind == "dl"
         else [t.shank_max_deg, t.depth_deg, t.bar_x_offset, t.lam, t.d])
    p = predict(key, kind, x, DEFAULT_HIP_AP_OFFSET)
    tb = table(key, kind)
    return p, (p - np.array(tb["mean"], float)) / np.array(tb["sd"], float)


def cop_penalty(key, kind, x, hx):
    """Residuos blandos: el CoP debe quedar dentro del 85 % central de la base de apoyo (equilibrio físicamente posible)."""
    if kind == "dl":
        a = from_raw(proxy_athlete(1.675, 76.0), hip_ap_offset=hx)
        load, ss = 130.0, (0.0, 0.5, 1.0)
    else:
        a = from_raw(proxy_athlete(1.803, 100.2), hip_ap_offset=hx)
        load, ss = 154.0, (0.5, 1.0)
    t = tech_from(key, kind, x)
    half = 0.85 * a.Lf * math.cos(math.radians(t.toe_out_deg)) / 2.0
    out = []
    for sv in ss:
        try:
            cop = solve(a, t, load, sv).cop_x
        except Infeasible:
            out.append(10.0)
            continue
        out.append(max(0.0, abs(cop) - half) / 0.01)
    return np.array(out)


def split(theta):
    parts, i = [], 0
    for key, kind, n in LAYOUT:
        parts.append((key, kind, theta[i:i + n]))
        i += n
    return parts, theta[i]


def table(key, kind):
    return (DL if kind == "dl" else SQ)[key]


def residuals(theta, means=None):
    parts, hx = split(theta)
    out = []
    for key, kind, x in parts:
        tb = table(key, kind)
        mean = np.array(tb["mean"] if means is None else means[key], float)
        p = predict(key, kind, x, hx)
        out.append(np.full(len(mean), 10.0) if p is None else (p - mean) / np.array(tb["sd"], float))
        out.append(cop_penalty(key, kind, x, hx))
    out.append(np.array([hx / HX_PRIOR_SD]))
    return np.concatenate(out)


SQ_FIXED = {}   # key -> (lam, d) fijos (posición de barra no identificable con estos datos)


def bounds():
    lo, hi = [], []
    for key, kind, n in LAYOUT:
        b = DL_BOUNDS if kind == "dl" else SQ_BOUNDS
        bl, bh = list(b[0]), list(b[1])
        if key in SQ_FIXED:
            lam, d = SQ_FIXED[key]
            bl[3], bh[3], bl[4], bh[4] = lam - 1e-4, lam + 1e-4, d - 1e-4, d + 1e-4
        lo += bl; hi += bh
    lo += HX_BOUNDS[0]; hi += HX_BOUNDS[1]
    return np.array(lo), np.array(hi)


def fit(theta0s, means=None):
    lo, hi = bounds()
    best = None
    for t0 in theta0s:
        s = least_squares(residuals, np.clip(t0, lo, hi), bounds=(lo, hi), args=(means,), x_scale=(hi - lo) / 5)
        if best is None or s.cost < best.cost:
            best = s
    return best


if __name__ == "__main__":
    import sys
    if "fixbar" in sys.argv:
        SQ_FIXED.update({"sentadilla_barra_alta": (1.0, 0.06), "sentadilla_barra_baja": (0.93, 0.07)})
    lo, hi = bounds()
    rng = np.random.default_rng(3)
    starts = [lo + (hi - lo) * rng.random(len(lo)) for _ in range(10)]
    sol = fit(starts)
    theta = sol.x
    parts, hx = split(theta)
    print(f"chi² total = {2 * sol.cost:.1f}   (31 objetivos + 10 penalizaciones de equilibrio + 1 a priori; 21 parámetros)")
    print(f"hx (centro de cadera respecto al vértice trocantérico) = {hx * 100:+.1f} cm")

    # rangos por remuestreo de los objetivos publicados dentro de su desviación estándar
    N = 30
    thetas = []
    for k in range(N):
        r = np.random.default_rng(100 + k)
        means = {key: r.normal(table(key, kind)["mean"], table(key, kind)["sd"]) for key, kind, _ in LAYOUT}
        thetas.append(fit([theta], means).x)
    thetas = np.array(thetas)
    plo, phi = np.percentile(thetas, 2.5, axis=0), np.percentile(thetas, 97.5, axis=0)

    result = {"hx": [hx, plo[-1], phi[-1]]}
    i = 0
    for key, kind, x in parts:
        names = DL_NAMES if kind == "dl" else SQ_NAMES
        labels = DL_LABELS if kind == "dl" else SQ_LABELS
        p = predict(key, kind, x, hx)
        tb = table(key, kind)
        print(f"\n== {key}")
        for j, nme in enumerate(names):
            print(f"   {nme:14s} = {x[j]:8.3f}   rango [{plo[i + j]:7.3f}, {phi[i + j]:7.3f}]")
        result[key] = {nme: [float(x[j]), float(plo[i + j]), float(phi[i + j])] for j, nme in enumerate(names)}
        print("   objetivo                                  modelo   publicado ± DE    z")
        for lab, m, mu, sd in zip(labels, p, tb["mean"], tb["sd"]):
            print(f"   {lab:38s} {m:8.1f}   {mu:8.1f} ± {sd:4.1f}   {(m - mu) / sd:+5.1f}{'  <--' if abs((m - mu) / sd) > 1.5 else ''}")
        i += len(names)
    result["hx"] = [float(hx), float(plo[-1]), float(phi[-1])]
    print(f"\nhx rango (P2.5–P97.5): [{plo[-1] * 100:+.1f}, {phi[-1] * 100:+.1f}] cm")

    # validación independiente: momentos pico de la fase 1 (Frontiers 2025: cadera conv/sumo 303/276; rodilla sumo/conv 80.2/52.8)
    a = from_raw(proxy_athlete(1.675, 76.0), hip_ap_offset=hx)
    pk = {}
    for key, kind, x in parts[:2]:
        t = tech_from(key, kind, x)
        skp = knee_passing_s(a, t, 130.0)
        ms = [solve(a, t, 130.0, s).moments for s in np.linspace(0, skp, 20)]
        pk[key] = (max(m["cadera"] for m in ms), max(m["rodilla"] for m in ms))
    c, s = pk["peso_muerto_convencional"], pk["peso_muerto_sumo"]
    print(f"VALIDACION independiente: cadera conv/sumo = {c[0] / s[0]:.2f} (lit. 1.10) ; rodilla sumo/conv = {s[1] / c[1]:.2f} (lit. 1.52)")
    print("\nVALIDACION trunk (no ajustado): ángulo del tronco sobre la horizontal en el despegue (Escamilla 2000: conv 11±7°, sumo 16±6°; definición inferida)")
    for key, kind, x in parts[:2]:
        a = from_raw(proxy_athlete(1.675, 76.0), hip_ap_offset=hx)
        r0 = solve(a, tech_from(key, kind, x), 130.0, 0.0)
        print(f"   {key}: {90 - r0.trunk_deg:.1f}°   CoP en el despegue: {r0.cop_x:+.3f} m")
    json.dump(result, open("calibration_result_fixbar.json" if SQ_FIXED else "calibration_result.json", "w"), indent=1)
