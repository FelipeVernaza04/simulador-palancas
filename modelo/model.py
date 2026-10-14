"""Modelo sagital 2D, cuasi-estático, simétrico, de momentos articulares netos.

Ejes: x hacia adelante, y hacia arriba, origen en el suelo bajo el medio pie.
Las piernas se resuelven en un plano inclinado (postura ancha): ver docs/protocolo y la sección 2.1 del modelo.
Momentos extensores positivos, en N·m, para la carga total (ambas piernas).
"""
from dataclasses import dataclass, field
import math
from anthropometry import Anthro, mass_points, G
from techniques import Technique

LUMBAR_MU = 0.10   # punto lumbar a esta fracción de Ltr sobre la cadera (APROXIMACIÓN)


class Infeasible(Exception):
    """La postura no es alcanzable con estas medidas y esta técnica."""


@dataclass
class Result:
    s: float
    joints: dict            # nombre -> (x, y) en el plano sagital
    bar: tuple
    psi_deg: float
    trunk_deg: float
    moments: dict           # tobillo, rodilla, cadera, lumbar (N·m)
    bar_moments: dict       # contribución de la barra
    body_moments: dict      # contribución del cuerpo
    arms: dict              # brazos de momento de la barra (m)
    shares: dict            # reparto relativo entre tobillo, rodilla y cadera
    cop_x: float
    cop_in_base: bool
    warnings: list = field(default_factory=list)


def _clampsqrt(v: float, tol: float = 1e-9) -> float:
    if v < -tol:
        raise Infeasible("radicando negativo")
    return math.sqrt(max(v, 0.0))


def _stance_c(a: Anthro, t: Technique) -> float:
    w = a.shoulder_w if (t.stance == "shoulder" and a.shoulder_w > 0) else t.w_frac * a.H
    return max((w - a.b) / 2.0, 0.0)


def _two_link(ax, ay, hx, hy, L1, L2):
    """Rodilla (la rama anterior) para un tobillo (ax, ay) y una cadera (hx, hy)."""
    dx, dy = hx - ax, hy - ay
    d = math.hypot(dx, dy)
    if d > L1 + L2 + 1e-9 or d < abs(L1 - L2) - 1e-9:
        raise Infeasible("cadera fuera del alcance de la pierna")
    aa = (L1 * L1 - L2 * L2 + d * d) / (2 * d)
    h = _clampsqrt(L1 * L1 - aa * aa)
    ex, ey = dx / d, dy / d
    cands = [(ax + aa * ex - h * ey, ay + aa * ey + h * ex), (ax + aa * ex + h * ey, ay + aa * ey - h * ex)]
    return max(cands, key=lambda p: p[0])


def _moments(a: Anthro, mb: float, psi: float, x_a, x_k, x_h, hip, shoulder, trunk_phi, bar_x, arms_x, toe_out_deg):
    mp = mass_points(a)
    u_h, u_k = x_h, x_k
    s_ph, c_ph = math.sin(trunk_phi), math.cos(trunk_phi)
    x_trunk = u_h + mp["f_trunk_from_hip"] * a.Ltr * s_ph
    x_head = shoulder[0] + mp["head_above_shoulder"] * s_ph
    x_thigh = u_h + mp["f_thigh"] * (u_k - u_h)
    x_shank = u_k + mp["f_shank"] * (x_a - u_k)
    x_L = u_h + LUMBAR_MU * a.Ltr * s_ph
    x_hc = x_h + a.hip_ap_offset   # centro de cadera: el vértice del modelo se desplaza hacia adelante en hip_ap_offset

    body_above_hip = [(mp["m_trunk"], x_trunk), (mp["m_head"], x_head), (mp["m_arms"], arms_x)]
    bar = [(mb, bar_x)]
    thighs = [(2 * mp["m_thigh"], x_thigh)]
    shanks = [(2 * mp["m_shank"], x_shank)]
    cpsi = math.cos(psi)

    def ext_fwd(items, xj):   # peso por delante de la articulación -> extensor (cadera, tobillo, lumbar)
        return G * sum(m * (x - xj) for m, x in items)

    def ext_bwd(items, xj):   # peso por detrás de la articulación -> extensor (rodilla)
        return G * sum(m * (xj - x) for m, x in items)

    def parts(items_bar, items_body):
        out = {}
        for tag, items in (("bar", items_bar), ("body", items_body)):
            out[tag] = {
                "cadera": ext_fwd(items["hip"], x_hc) * cpsi,
                "rodilla": ext_bwd(items["knee"], x_k) * cpsi,
                "tobillo": ext_fwd(items["ankle"], x_a) * cpsi,
                "lumbar": ext_fwd(items["lumbar"], x_L),
            }
        return out

    bar_items = {"hip": bar, "knee": bar, "ankle": bar, "lumbar": bar}
    body_items = {"hip": body_above_hip, "knee": body_above_hip + thighs,
                  "ankle": body_above_hip + thighs + shanks, "lumbar": body_above_hip}
    pr = parts(bar_items, body_items)
    total = {k: pr["bar"][k] + pr["body"][k] for k in pr["bar"]}
    arms_bar = {"cadera": bar_x - x_hc, "rodilla": x_k - bar_x, "tobillo": bar_x - x_a, "lumbar": bar_x - x_L}

    # equilibrio: CoP = x del centro de masa total (incluye los pies)
    x_foot = -a.Lf / 2.0 + mp["f_foot"] * a.Lf
    allm = body_above_hip + bar + thighs + shanks + [(2 * mp["m_foot"], x_foot)]
    cop = sum(m * x for m, x in allm) / sum(m for m, _ in allm)
    half = a.Lf * math.cos(math.radians(toe_out_deg)) / 2.0
    return total, pr["bar"], pr["body"], arms_bar, cop, (-half <= cop <= half)


def _pack(a, t, s, mb, ankle, knee, hip, shoulder, bar_xy, psi, phi, bar_x, arms_x):
    x_a, x_k, x_h = ankle[0], knee[0], hip[0]
    total, barm, bodym, arms, cop, ok = _moments(a, mb, psi, x_a, x_k, x_h, hip, shoulder, phi, bar_x, arms_x,
                                                 t.toe_out_deg)
    legs = {k: abs(total[k]) for k in ("tobillo", "rodilla", "cadera")}
    tot = sum(legs.values()) or 1.0
    shares = {k: v / tot for k, v in legs.items()}
    warnings = []
    if not ok:
        warnings.append(f"CoP fuera de la base de apoyo ({cop:+.3f} m)")
    return Result(s=s, joints={"tobillo": ankle, "rodilla": knee, "cadera": hip, "hombro": shoulder},
                  bar=bar_xy, psi_deg=math.degrees(psi), trunk_deg=math.degrees(phi), moments=total,
                  bar_moments=barm, body_moments=bodym, arms=arms, shares=shares, cop_x=cop,
                  cop_in_base=ok, warnings=warnings)


def solve_squat(a: Anthro, t: Technique, load_kg: float, s: float) -> Result:
    """Sentadilla. La pierna se parametriza por el ángulo del muslo respecto a la vertical (beta), lineal en s, entre la postura de
    pie (beta = inclinación inicial de la tibia, cadera detrás de la rodilla) y el fondo (muslo a `depth_deg` bajo la horizontal).
    Los extremos s=0 y s=1 coinciden con la parametrización anterior por altura de cadera; solo cambia el interior del recorrido
    (evita el salto de la posición horizontal de la cadera al inicio, donde el muslo está casi vertical)."""
    c = _stance_c(a, t)
    x_a, h_a = -a.dam, a.ha
    ph0, phm = math.radians(t.shank_top_deg), math.radians(t.shank_max_deg)
    al = math.radians(t.depth_deg)

    beta = ph0 + s * ((math.pi / 2 + al) - ph0)
    phs = ph0 + s * (phm - ph0)
    u_k, v_k = x_a + a.Ls * math.sin(phs), a.Ls * math.cos(phs)
    v_h = v_k + a.Lt * math.cos(beta)
    u_h = u_k - a.Lt * math.sin(beta)
    yp = _clampsqrt(v_h ** 2 - c ** 2)
    psi = math.atan2(c, yp)

    A, B = t.lam * a.Ltr, t.d
    if t.bar_x_mode == "erect_top":
        v_top = (a.Ls + a.Lt) * math.cos(ph0)
        u_k0, v_k0 = x_a + a.Ls * math.sin(ph0), a.Ls * math.cos(ph0)
        u_h0 = u_k0 - _clampsqrt(a.Lt ** 2 - (v_top - v_k0) ** 2)
        ph_t = math.radians(t.trunk_top_deg)
        x_b = u_h0 + A * math.sin(ph_t) - B * math.cos(ph_t)
    else:
        x_b = t.bar_x_offset
    C, R = x_b - u_h, math.hypot(A, B)
    if abs(C) > R:
        raise Infeasible("la barra no puede quedar sobre la vertical elegida con este tronco")
    phi = math.atan2(B, A) + math.asin(C / R)

    hip = (u_h, h_a + yp)
    knee = (u_k, h_a + v_k * math.cos(psi))
    shoulder = (u_h + a.Ltr * math.sin(phi), hip[1] + a.Ltr * math.cos(phi))
    bar_xy = (u_h + A * math.sin(phi) - B * math.cos(phi), hip[1] + A * math.cos(phi) + B * math.sin(phi))
    return _pack(a, t, s, load_kg, (x_a, h_a), knee, hip, shoulder, bar_xy, psi, phi, bar_xy[0], shoulder[0])


def solve_deadlift(a: Anthro, t: Technique, load_kg: float, s: float) -> Result:
    c = _stance_c(a, t)
    x_a, h_a = -a.dam, a.ha
    x_b = t.bar_x_offset
    u_top = x_b + t.eps_top
    v_top = _clampsqrt((a.Ls + a.Lt) ** 2 - (u_top - x_a) ** 2)
    yp_top = _clampsqrt(v_top ** 2 - c ** 2)
    yp0 = t.f_h0 * yp_top
    yp = yp0 + (yp_top - yp0) * (s ** t.p)
    eps = t.eps0 + (t.eps_top - t.eps0) * s

    y_h_top = h_a + yp_top
    y_b_top = y_h_top + a.Ltr - _clampsqrt(a.La ** 2 - t.eps_top ** 2)
    y_b = t.plate_radius + s * (y_b_top - t.plate_radius)

    x_sh = x_b + eps
    y_sh = y_b + _clampsqrt(a.La ** 2 - eps ** 2)
    y_h = h_a + yp
    dy = y_sh - y_h
    if dy < 0 or dy > a.Ltr + 1e-9:
        raise Infeasible("el hombro no es alcanzable desde esa altura de cadera")
    x_h = x_sh - _clampsqrt(a.Ltr ** 2 - dy ** 2)
    phi = math.atan2(x_sh - x_h, dy)
    psi = math.atan2(c, yp)
    v_h = math.hypot(yp, c)
    ku, kv = _two_link(x_a, 0.0, x_h, v_h, a.Ls, a.Lt)
    hip = (x_h, y_h)
    knee = (ku, h_a + kv * math.cos(psi))
    shoulder = (x_sh, y_sh)
    from anthropometry import mass_points as _mp
    f_arm = _mp(a)["f_arm_com"]
    arms_x = x_sh + f_arm * (x_b - x_sh)
    return _pack(a, t, s, load_kg, (x_a, h_a), knee, hip, shoulder, (x_b, y_b), psi, phi, x_b, arms_x)


def solve(a: Anthro, t: Technique, load_kg: float, s: float) -> Result:
    return (solve_squat if t.exercise == "squat" else solve_deadlift)(a, t, load_kg, s)
