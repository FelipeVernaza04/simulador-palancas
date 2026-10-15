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
KNEE_FRONT = 0.05  # centro articular de la rodilla -> cara anterior de la rótula, con partes blandas (m) (APROXIMACIÓN)
BAR_RADIUS = 0.014 # barra olímpica de 28 mm
KNEE_CLEAR = KNEE_FRONT + BAR_RADIUS   # distancia mínima centro de barra - centro de rodilla en el plano sagital
SHIN_FRONT = 0.02  # eje tobillo-rodilla -> cara anterior de la espinilla (cresta tibial casi subcutánea) (m) (APROXIMACIÓN)
THIGH_FRONT = 0.07 # eje rodilla-cadera -> cara anterior del muslo a media altura del muslo (m) (APROXIMACIÓN)
BAR_TO_HIP = 0.104 # bloqueo: centro de la cadera -> centro de la barra, que apoya en la cara anterior del muslo (m) (APROXIMACIÓN)
BAR_RETURN_P = 2.0 # la barra vuelve de la vertical sobre el pie hacia su x de bloqueo como s^p (más tarde cuanto mayor)
LOCKOUT_BLEND_S = 0.85  # desde este s el tronco pasa del que dicta la rodilla al del bloqueo (solo con bar_top_from_ankle)


class Infeasible(Exception):
    """La postura no es alcanzable con estas medidas y esta técnica."""


@dataclass
class Result:
    s: float
    joints: dict            # nombre -> (x, y) en el plano sagital
    lateral: dict           # tobillo, rodilla, cadera -> distancia lateral al plano medio (m), pierna derecha (+ hacia afuera)
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


def _knee_3d(a: Anthro, ankle, hip, out_deg):
    """Rodilla de la pierna derecha para un tobillo y una cadera 3D, puntos (x adelante, y arriba, z lateral hacia afuera).

    Las longitudes de tibia y muslo fijan una circunferencia de posiciones posibles alrededor del eje tobillo-cadera. La rodilla se
    abre en el plano que contiene ese eje y la dirección horizontal (cos th, 0, sin th), con th = `out_deg` hacia afuera: la rodilla
    sigue la dirección del pie. Con th = 0 es el plano que contiene al eje x (el modelo plano original: misma rodilla que antes).
    """
    d = [h - t for h, t in zip(hip, ankle)]
    D = math.sqrt(sum(c * c for c in d))
    if D > a.Ls + a.Lt + 1e-9 or D < abs(a.Ls - a.Lt) - 1e-9:
        raise Infeasible("cadera fuera del alcance de la pierna")
    w = [c / D for c in d]
    aa = (a.Ls * a.Ls - a.Lt * a.Lt + D * D) / (2 * D)
    h = _clampsqrt(a.Ls * a.Ls - aa * aa)
    th = math.radians(out_deg)
    u = (math.cos(th), 0.0, math.sin(th))
    pr = sum(ui * wi for ui, wi in zip(u, w))
    n = [ui - pr * wi for ui, wi in zip(u, w)]
    nn = math.sqrt(sum(c * c for c in n))
    n = [c / nn for c in n] if nn > 1e-9 else [1.0, 0.0, 0.0]
    return tuple(ankle[i] + aa * w[i] + h * n[i] for i in range(3))


def _seg_dist(px, py, a_, b_):
    """Distancia del punto (px, py) al segmento a_-b_ y parámetro t de su pie sobre el segmento (0 en a_, 1 en b_)."""
    dx, dy = b_[0] - a_[0], b_[1] - a_[1]
    t = max(0.0, min(1.0, ((px - a_[0]) * dx + (py - a_[1]) * dy) / (dx * dx + dy * dy)))
    return math.hypot(px - a_[0] - t * dx, py - a_[1] - t * dy), t


def _leg_gap(px, py, ankle, knee, hip):
    """Distancia (m) del punto sagital (px, py) a la cara anterior de la pierna: espinilla que engrosa hasta la rótula (radio
    SHIN_FRONT -> KNEE_FRONT, cúbica) y mitad inferior del muslo (KNEE_FRONT -> THIGH_FRONT). Es continua, sin esquinas, y en la
    rodilla vale distancia al centro - KNEE_FRONT. No incluye el muslo alto: ahí la barra va con las manos al costado del muslo."""
    d, t = _seg_dist(px, py, ankle, knee)
    gap = d - (SHIN_FRONT + (KNEE_FRONT - SHIN_FRONT) * t ** 3)
    mid = (0.5 * (knee[0] + hip[0]), 0.5 * (knee[1] + hip[1]))
    d, t = _seg_dist(px, py, knee, mid)
    return min(gap, d - (KNEE_FRONT + (THIGH_FRONT - KNEE_FRONT) * t))


def _bar_x_clearing_leg(x0, y_b, ankle, knee, hip):
    """Menor x >= x0 con la barra (centro en (x, y_b), radio BAR_RADIUS) fuera de la pierna: la barra solo se adelanta."""
    if _leg_gap(x0, y_b, ankle, knee, hip) >= BAR_RADIUS:
        return x0
    lo, hi = x0, x0 + 0.15
    for _ in range(22):      # precisión 0,15 / 2^22 = 4e-8 m, bajo la tolerancia del lazo que la llama
        mid = 0.5 * (lo + hi)
        if _leg_gap(mid, y_b, ankle, knee, hip) >= BAR_RADIUS:
            hi = mid
        else:
            lo = mid
    return hi


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


def _pack(a, t, s, mb, ankle, knee, hip, shoulder, bar_xy, psi, phi, bar_x, arms_x, lateral):
    x_a, x_k, x_h = ankle[0], knee[0], hip[0]
    total, barm, bodym, arms, cop, ok = _moments(a, mb, psi, x_a, x_k, x_h, hip, shoulder, phi, bar_x, arms_x,
                                                 t.toe_out_deg)
    legs = {k: abs(total[k]) for k in ("tobillo", "rodilla", "cadera")}
    tot = sum(legs.values()) or 1.0
    shares = {k: v / tot for k, v in legs.items()}
    warnings = []
    if not ok:
        warnings.append(f"CoP fuera de la base de apoyo ({cop:+.3f} m)")
    return Result(s=s, joints={"tobillo": ankle, "rodilla": knee, "cadera": hip, "hombro": shoulder}, lateral=lateral,
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
    z_h = a.b / 2.0
    z_a = z_h + c
    # la cadera y el tobillo vienen del plano de la pierna; la rodilla se abre según t.knee_out_deg (ver _knee_3d)
    kx, ky, kz = _knee_3d(a, (x_a, h_a, z_a), (u_h, hip[1], z_h), t.knee_out_deg)
    knee = (kx, ky)
    shoulder = (u_h + a.Ltr * math.sin(phi), hip[1] + a.Ltr * math.cos(phi))
    bar_xy = (u_h + A * math.sin(phi) - B * math.cos(phi), hip[1] + A * math.cos(phi) + B * math.sin(phi))
    lateral = {"tobillo": z_a, "rodilla": kz, "cadera": z_h}
    return _pack(a, t, s, load_kg, (x_a, h_a), knee, hip, shoulder, bar_xy, psi, phi, bar_xy[0], shoulder[0], lateral)


def solve_deadlift(a: Anthro, t: Technique, load_kg: float, s: float) -> Result:
    """Peso muerto. El recorrido de la pierna lo fija la posición horizontal de la rodilla (como la sentadilla fija la inclinación de
    la tibia): `knee0` es la rodilla respecto de la barra al despegar (+ = rodilla por delante de la barra) y la rodilla llega al
    bloqueo con la pierna recta; entre ambos, kx(s) = kx_top + (kx_0 - kx_top) (1 - s)^knee_p. La altura de la cadera es consecuencia:
    se busca (bisección) la altura de cadera, con el hombro donde lo pone la ficha, que deja la rodilla en kx(s). Antes la altura de
    la cadera era el parámetro y la rodilla quedaba casi quieta (la tibia cambiaba ~3° en todo el tirón)."""
    c = _stance_c(a, t)
    x_a, h_a = -a.dam, a.ha
    x_b = t.bar_x_offset
    z_h = a.b / 2.0
    z_a = z_h + c
    # Bloqueo. Con bar_top_from_ankle la barra termina en x_bar_top (cerca del tobillo), la cadera una profundidad de muslo detrás
    # de ella y el tronco inclinado `layback_deg` hacia atrás (hombro detrás de la cadera); si no, el de siempre (tronco vertical,
    # hombro a eps_top de la barra, todo sobre bar_x_offset).
    layback = t.bar_top_from_ankle is not None
    x_bar_top = x_a + t.bar_top_from_ankle if layback else x_b
    phi_top = -math.radians(t.layback_deg) if layback else 0.0
    u_top = x_bar_top - BAR_TO_HIP if layback else x_b + t.eps_top
    eps_top = (u_top + a.Ltr * math.sin(phi_top) - x_bar_top) if layback else t.eps_top
    v_top = _clampsqrt((a.Ls + a.Lt) ** 2 - (u_top - x_a) ** 2)
    yp_top = _clampsqrt(v_top ** 2 - c ** 2)
    y_h_top = h_a + yp_top
    y_b_top = y_h_top + a.Ltr * math.cos(phi_top) - _clampsqrt(a.La ** 2 - eps_top ** 2)
    y_b = t.plate_radius + s * (y_b_top - t.plate_radius)
    eps = t.eps0 + (eps_top - t.eps0) * s
    x_bar0 = x_b + (x_bar_top - x_b) * s ** BAR_RETURN_P    # la barra sigue la vertical hasta que la pierna la adelanta
    x_sh = x_bar0 + eps   # el hombro queda donde lo pone la ficha de técnica (hombro por delante de la barra en `eps`)

    # rodilla objetivo: al bloqueo es la de la pierna recta (hombro y cadera en u_top)
    kx_top = _knee_3d(a, (x_a, h_a, z_a), (u_top, y_h_top, z_h), t.knee_out_deg)[0]
    kx_0 = x_b + t.knee0
    kx_t = kx_top + (kx_0 - kx_top) * (1.0 - s) ** t.knee_p

    def trunk_for(y_sh, kx_t):
        """Inclinación del tronco con el hombro en (x_sh, y_sh) y la rodilla en x = kx_t.

        La incógnita es la inclinación del tronco phi respecto de la vertical (cadera = hombro - Ltr (sin phi, cos phi)), no la altura
        de la cadera: con el tronco casi vertical la altura es hipersensible (un milímetro de altura mueve la cadera centímetros) y
        el rango de alturas que acotaba la búsqueda dejaba fuera la raíz justo antes del bloqueo. En phi la raíz está bien
        condicionada hasta phi = 0, y el límite superior es el alcance de la pierna (pierna recta)."""
        def at(phi):
            xh, yh = x_sh - a.Ltr * math.sin(phi), y_sh - a.Ltr * math.cos(phi)
            kn = _knee_3d(a, (x_a, h_a, z_a), (xh, yh, z_h), t.knee_out_deg)
            if kn[1] > yh:   # en peso muerto la cadera queda sobre la rodilla (muslo por encima de la horizontal)
                raise Infeasible("el muslo quedaría bajo la horizontal")
            return xh, yh, kn

        def f(phi):          # kx - kx_t: decrece con phi (la cadera sube y se retrasa)
            return at(phi)[2][0] - kx_t

        lo = None
        N, phi_cap = 48, 1.4
        phi_min = min(0.0, phi_top) - 0.15 if layback else 0.0   # con extensión de espalda la cadera puede quedar delante del hombro
        for i in range(N + 1):
            phi = phi_min + (phi_cap - phi_min) * i / N
            try:
                fv = f(phi)
            except Infeasible:
                if lo is None:
                    continue          # cadera demasiado baja: el muslo bajo la horizontal
                # se salió del alcance de la pierna: el borde (pierna recta) es el último phi alcanzable
                a_, b_ = lo[0], phi
                for _ in range(50):
                    mid = 0.5 * (a_ + b_)
                    try:
                        f(mid)
                        a_ = mid
                    except Infeasible:
                        b_ = mid
                fe = f(a_)
                if fe > 1e-6:
                    raise Infeasible("la rodilla no alcanza esa posición con estas medidas")
                lo, hi = (lo[0], lo[1]), (a_, fe)
                break
            if abs(fv) < 1e-12:
                return phi
            if lo is not None and (lo[1] > 0) != (fv > 0):
                hi = (phi, fv)
                break
            lo = (phi, fv)
        else:
            raise Infeasible("la rodilla no alcanza esa posición con estas medidas")
        a_, fa, b_ = lo[0], lo[1], hi[0]
        for _ in range(60):
            mid = 0.5 * (a_ + b_)
            try:
                fm = f(mid)
            except Infeasible:
                b_ = mid
                continue
            if (fm > 0) == (fa > 0):
                a_, fa = mid, fm
            else:
                b_ = mid
        phi = 0.5 * (a_ + b_)
        # La cadera no sube por encima de su altura al bloqueo. Con la rodilla clavada en su posición de pierna recta, justo antes del
        # bloqueo la raíz estaría con la cadera unos mm más arriba y ~6 cm atrás (tronco inclinado), y al llegar a s = 1 la solución
        # saltaría a la otra rama (tronco vertical, cadera adelante). En su lugar la cadera se queda a esa altura y la rodilla es
        # consecuencia: el tronco llega a la vertical de forma continua.
        if y_sh - a.Ltr * math.cos(phi) > y_h_top:
            phi = math.copysign(math.acos(min(1.0, (y_sh - y_h_top) / a.Ltr)), phi)
        return phi

    def place(y_sh, phi):
        """Cadera (x, y) con el tronco inclinado phi y la rodilla 3D resultante."""
        xh, yh = x_sh - a.Ltr * math.sin(phi), y_sh - a.Ltr * math.cos(phi)
        kn = _knee_3d(a, (x_a, h_a, z_a), (xh, yh, z_h), t.knee_out_deg)
        if kn[1] > yh:   # en peso muerto la cadera queda sobre la rodilla (muslo por encima de la horizontal)
            raise Infeasible("el muslo quedaría bajo la horizontal")
        return xh, yh, kn

    end_shift = []

    def hip_for(y_sh):
        """Cadera (x, y) con el hombro en (x_sh, y_sh) y la rodilla en x = kx_t, y la rodilla 3D resultante."""
        phi = trunk_for(y_sh, kx_t)
        if layback and s > LOCKOUT_BLEND_S:
            # Extensión de la espalda. Con la rodilla en su posición de bloqueo el tronco que sale de la geometría de la pierna es el
            # `phi_end`, no el pedido (phi_top). Esa diferencia se suma de forma gradual (suavizada) desde LOCKOUT_BLEND_S: así el
            # tronco llega a phi_top al bloqueo sin que la cadera se adelante antes que la rodilla se extienda.
            if not end_shift:
                y_sh_end = y_h_top + a.Ltr * math.cos(phi_top)
                end_shift.append(phi_top - trunk_for(y_sh_end, kx_top))
            w = (s - LOCKOUT_BLEND_S) / (1.0 - LOCKOUT_BLEND_S)
            phi += w * w * (3.0 - 2.0 * w) * end_shift[0]
            if y_sh - a.Ltr * math.cos(phi) > y_h_top:     # la cadera no sube por encima de su altura al bloqueo
                phi = math.copysign(math.acos(min(1.0, (y_sh - y_h_top) / a.Ltr)), phi)
        return place(y_sh, phi)

    # La barra no puede atravesar la pierna (espinilla, rótula, muslo bajo: ver _leg_gap): si la vertical la cruza, la barra se
    # adelanta lo justo respecto del cuerpo (hombro, cadera y rodilla no se mueven: los brazos se inclinan un poco hacia adelante).
    # El contorno es continuo (la espinilla se une a la rótula y esta al muslo con pendientes suaves), así que el desvío no tiene
    # esquinas. La altura del hombro depende de la inclinación del brazo, y con ella la cadera, así que se itera (la rodilla en x
    # ya es kx_t).
    xbe = x_bar0
    for _ in range(12):
        e = x_sh - xbe
        y_sh = y_b + _clampsqrt(a.La ** 2 - e ** 2)
        x_h, y_h, (kx, ky, kz) = hip_for(y_sh)
        need = _bar_x_clearing_leg(xbe, y_b, (x_a, h_a), (kx, ky), (x_h, y_h))
        if need <= xbe + 1e-7:
            break
        xbe = need
    psi = math.atan2(c, y_h - h_a)
    phi = math.atan2(x_sh - x_h, y_sh - y_h)
    hip = (x_h, y_h)
    knee = (kx, ky)
    shoulder = (x_sh, y_sh)
    from anthropometry import mass_points as _mp
    f_arm = _mp(a)["f_arm_com"]
    arms_x = x_sh + f_arm * (xbe - x_sh)
    lateral = {"tobillo": z_a, "rodilla": kz, "cadera": z_h}
    return _pack(a, t, s, load_kg, (x_a, h_a), knee, hip, shoulder, (xbe, y_b), psi, phi, xbe, arms_x, lateral)


def solve(a: Anthro, t: Technique, load_kg: float, s: float) -> Result:
    return (solve_squat if t.exercise == "squat" else solve_deadlift)(a, t, load_kg, s)
