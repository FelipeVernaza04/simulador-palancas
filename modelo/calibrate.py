"""Calibración de las fichas de técnica contra la literatura (ajuste global por mínimos cuadrados ponderados).

Peso muerto: Escamilla y cols. (2000), 40 powerlifters hombres (~167 cm, 76 kg).
Sentadilla:  Swinton y cols. (2012), 12 powerlifters hombres (180,3 cm, 100,2 kg, 1RM 220,2 kg), 70 % 1RM, inicio de la fase concéntrica.

Se simula un "atleta proxy" con las proporciones de SAMPLE_MALE escaladas a la estatura media del estudio, y se buscan los
parámetros de estilo que reproducen TODAS las magnitudes publicadas a la vez, ponderadas por su desviación estándar.
Los rangos de los parámetros salen de repetir el ajuste con los objetivos muestreados dentro de su desviación estándar.

Supuestos de interpretación (a verificar en los artículos originales):
 - "ángulo de muslo" de Escamilla = 180° - inclinación sobre la horizontal.
 - los ángulos de segmento de Escamilla son 3D (respecto al plano horizontal).
 - "brazo de momento" = distancia horizontal de la barra a la articulación, positiva si la barra está por delante.
 - "ancho de postura" se usa como distancia entre centros de tobillos (la definición exacta no está en lo que pude leer).
"""
from dataclasses import replace
import math
import numpy as np
from scipy.optimize import least_squares

from anthropometry import SAMPLE_MALE, from_raw, RawMeasures
from techniques import TECHNIQUES
from model import solve, Infeasible

LEN_FIELDS = ("H", "h_troc", "Lt_raw", "Ls_raw", "ha", "Lf", "heel_to_malleolus", "Ltr_raw", "h_acro", "La_raw", "D_asis")


def proxy_athlete(H: float, mass: float) -> RawMeasures:
    k = H / SAMPLE_MALE.H
    # D_acro del proxy = 0,23·H (ASUNCIÓN propia, sin fuente verificada): fija el ancho de postura estándar (tobillos bajo los hombros)
    return replace(SAMPLE_MALE, mass=mass, D_acro=0.23 * H, **{f: getattr(SAMPLE_MALE, f) * k for f in LEN_FIELDS})


# --- geometría auxiliar ------------------------------------------------------------------------------------------------------
def incl_3d(r, p, q, h_a):
    """Inclinación sobre la horizontal (°) del segmento p->q con las piernas en el plano inclinado."""
    psi = math.radians(r.psi_deg)
    (x1, y1), (x2, y2) = r.joints[p], r.joints[q]
    dz = ((y2 - y1) / math.cos(psi)) * math.sin(psi)
    return math.degrees(math.atan2(y2 - y1, math.hypot(x2 - x1, dz)))


def knee_flexion_inplane(r, a):
    psi = math.radians(r.psi_deg)
    (xa, _), (xk, yk), (xh, yh) = r.joints["tobillo"], r.joints["rodilla"], r.joints["cadera"]
    vk, vh = (yk - a.ha) / math.cos(psi), (yh - a.ha) / math.cos(psi)
    u, w = (xa - xk, -vk), (xh - xk, vh - vk)
    c = (u[0] * w[0] + u[1] * w[1]) / (math.hypot(*u) * math.hypot(*w))
    return 180.0 - math.degrees(math.acos(max(-1.0, min(1.0, c))))


def knee_passing_s(a, tech, load):
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        r = solve(a, tech, load, mid)
        lo, hi = (mid, hi) if r.bar[1] < r.joints["rodilla"][1] else (lo, mid)
    return (lo + hi) / 2


# --- peso muerto -------------------------------------------------------------------------------------------------------------
DL_LABELS = ["brazo cadera LO", "brazo rodilla LO", "brazo tobillo LO", "muslo sobre horizontal LO (°)",
             "tibia sobre horizontal KP (°)", "brazo cadera KP", "brazo rodilla KP", "brazo tobillo KP", "recorrido barra"]
DL = {
    "peso_muerto_convencional": dict(
        mean=[27.9, -3.3, 6.1, 44, 83, 21.1, 3.2, 2.9, 51], sd=[5.2, 5.2, 3.7, 6, 4, 3.6, 3.2, 4.0, 6],
        w_frac=37 / 167, toe=0),
    "peso_muerto_sumo": dict(
        mean=[27.9, -13.2, -4.3, 37, 74, 23.5, -4.8, -5.2, 47], sd=[6.2, 8.0, 5.0, 6, 6, 5.4, 5.4, 5.3, 5],
        w_frac=73 / 168, toe=25),
}
DL_PARAMS = ("f_h0", "eps0", "eps_top", "p", "bar_x_offset")
DL_BOUNDS = ([0.35, -0.02, -0.10, 0.5, -0.25], [0.90, 0.15, 0.02, 2.0, 0.10])


def dl_tech(key, x):
    d = DL[key]
    return replace(TECHNIQUES[key], w_frac=d["w_frac"], toe_out_deg=d["toe"], f_h0=x[0], eps0=x[1], eps_top=x[2], p=x[3],
                   bar_x_offset=x[4], calibrated=True)


def dl_pred(key, x):
    a = from_raw(proxy_athlete(1.675, 76.0))
    t = dl_tech(key, x)
    load = 130.0
    try:
        r0 = solve(a, t, load, 0.0)
        skp = knee_passing_s(a, t, load)
        rk = solve(a, t, load, skp)
        r1 = solve(a, t, load, 1.0)
    except Infeasible:
        return None
    return np.array([r0.arms["cadera"] * 100, -r0.arms["rodilla"] * 100, r0.arms["tobillo"] * 100,
                     incl_3d(r0, "rodilla", "cadera", a.ha), incl_3d(rk, "tobillo", "rodilla", a.ha),
                     rk.arms["cadera"] * 100, -rk.arms["rodilla"] * 100, rk.arms["tobillo"] * 100,
                     (r1.bar[1] - r0.bar[1]) * 100])


# --- sentadilla --------------------------------------------------------------------------------------------------------------
SQ_LABELS = ["brazo cadera", "brazo rodilla", "brazo tobillo", "flexión de rodilla (°)", "flexión de tobillo (°)", "flexión de torso (°)"]
SQ = {
    "sentadilla_barra_alta": dict(mean=[25.2, -10.1, 9.9, 121.1, 37.2, 33.5], sd=[2.9, 1.1, 2.2, 3.4, 3.9, 4.6], w_frac=48.3 / 180.3),
    "sentadilla_barra_baja": dict(mean=[26.2, -8.1, 5.6, 112.1, 26.7, 33.1], sd=[2.1, 0.8, 1.6, 4.3, 5.1, 4.5], w_frac=89.6 / 180.3),
}
SQ_PARAMS = ("shank_max_deg", "depth_deg", "bar_x_offset")
SQ_BOUNDS = ([10.0, -15.0, -0.15], [50.0, 25.0, 0.15])


def sq_tech(key, x):
    return replace(TECHNIQUES[key], w_frac=SQ[key]["w_frac"], bar_x_mode="midfoot", shank_max_deg=x[0], depth_deg=x[1],
                   bar_x_offset=x[2], calibrated=True)


def sq_pred(key, x):
    a = from_raw(proxy_athlete(1.803, 100.2))
    t = sq_tech(key, x)
    try:
        r = solve(a, t, 154.0, 1.0)
    except Infeasible:
        return None
    return np.array([r.arms["cadera"] * 100, -r.arms["rodilla"] * 100, r.arms["tobillo"] * 100,
                     knee_flexion_inplane(r, a), x[0], r.trunk_deg])


# --- ajuste ------------------------------------------------------------------------------------------------------------------
def fit_general(pred, key, table, bounds, starts, mean=None):
    mean = np.array(table[key]["mean"] if mean is None else mean, dtype=float)
    sd = np.array(table[key]["sd"], dtype=float)
    lo, hi = np.array(bounds[0], float), np.array(bounds[1], float)

    def resid(x):
        p = pred(key, x)
        return np.full(len(mean), 10.0) if p is None else (p - mean) / sd

    best = None
    for x0 in starts:
        s = least_squares(resid, np.clip(x0, lo, hi), bounds=(lo, hi), x_scale=(hi - lo) / 5)
        if best is None or s.cost < best.cost:
            best = s
    return best


def starts_for(bounds, n=24, seed=0):
    rng = np.random.default_rng(seed)
    lo, hi = np.array(bounds[0], float), np.array(bounds[1], float)
    return [lo + (hi - lo) * rng.random(len(lo)) for _ in range(n)]


def param_ranges(pred, key, table, bounds, x_best, n=60, seed=1):
    rng = np.random.default_rng(seed)
    mean, sd = np.array(table[key]["mean"], float), np.array(table[key]["sd"], float)
    xs = []
    for _ in range(n):
        m = rng.normal(mean, sd)
        xs.append(fit_general(pred, key, table, bounds, [x_best], mean=m).x)
    xs = np.array(xs)
    return np.percentile(xs, 2.5, axis=0), np.percentile(xs, 97.5, axis=0)


def report(kind, key, pred, table, labels, params, bounds, with_ranges=True):
    sol = fit_general(pred, key, table, bounds, starts_for(bounds))
    p = pred(key, sol.x)
    lo, hi = param_ranges(pred, key, table, bounds, sol.x) if with_ranges else (sol.x, sol.x)
    print(f"\n== {key}   chi² = {2 * sol.cost:.1f}  ({len(p)} objetivos, {len(sol.x)} parámetros)")
    for n, v, a, b in zip(params, sol.x, lo, hi):
        print(f"   {n:14s} = {v:8.3f}   rango (P2.5–P97.5): [{a:7.3f}, {b:7.3f}]")
    print("   objetivo                                  modelo   publicado ± DE    z")
    for lab, m, mu, sd in zip(labels, p, table[key]["mean"], table[key]["sd"]):
        flag = "  <--" if abs((m - mu) / sd) > 1.5 else ""
        print(f"   {lab:38s} {m:8.1f}   {mu:8.1f} ± {sd:4.1f}   {(m - mu) / sd:+5.1f}{flag}")
    return sol.x, lo, hi


if __name__ == "__main__":
    res = {}
    for key in DL:
        res[key] = report("dl", key, dl_pred, DL, DL_LABELS, DL_PARAMS, DL_BOUNDS)
    for key in SQ:
        res[key] = report("sq", key, sq_pred, SQ, SQ_LABELS, SQ_PARAMS, SQ_BOUNDS)

    # validación independiente (no usada en el ajuste): cociente de momentos máximos de la fase 1, Frontiers 2025
    a = from_raw(proxy_athlete(1.675, 76.0))
    pk = {}
    for key in DL:
        t = dl_tech(key, res[key][0])
        skp = knee_passing_s(a, t, 130.0)
        ms = [solve(a, t, 130.0, s).moments for s in np.linspace(0, skp, 20)]
        pk[key] = (max(m["cadera"] for m in ms), max(m["rodilla"] for m in ms))
    c, s = pk["peso_muerto_convencional"], pk["peso_muerto_sumo"]
    print(f"\nVALIDACION independiente (picos de la fase 1): cadera conv/sumo = {c[0] / s[0]:.2f} (Frontiers 2025: 303/276 = 1.10);"
          f" rodilla sumo/conv = {s[1] / c[1]:.2f} (80.2/52.8 = 1.52)")
    np.save("/tmp/claude-0/-home-claude/963f9c7d-cd41-534c-a098-696a25da88b6/scratchpad/calib.npy",
            {k: [np.array(v[0]), np.array(v[1]), np.array(v[2])] for k, v in res.items()}, allow_pickle=True)
