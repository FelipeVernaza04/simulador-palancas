import math
from dataclasses import replace
import numpy as np
import pytest

from anthropometry import (SAMPLE_MALE, from_raw, closure_checks, total_mass_check, mass_points, G, DE_LEVA)
from techniques import TECHNIQUES
from model import solve, Infeasible, KNEE_CLEAR, SHIN_FRONT, BAR_RADIUS, _knee_3d, _two_link, _seg_dist
from uncertainty import draw_perturbations, monte_carlo

A = from_raw(SAMPLE_MALE)
S_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]


def _p3(r, name):
    """Punto 3D (x adelante, y arriba, z lateral hacia afuera) de una articulación de la pierna derecha."""
    x, y = r.joints[name]
    return x, y, r.lateral[name]


def _dist(p, q):
    return math.dist(p, q)


@pytest.mark.parametrize("key", list(TECHNIQUES))
def test_segment_lengths_preserved(key):
    t = TECHNIQUES[key]
    for s in S_GRID:
        r = solve(A, t, 100, s)
        # tibia y muslo conservan su longitud real en 3D (la rodilla puede abrirse fuera del plano tobillo-cadera-eje x)
        assert _dist(_p3(r, "tobillo"), _p3(r, "rodilla")) == pytest.approx(A.Ls, abs=1e-6)
        assert _dist(_p3(r, "rodilla"), _p3(r, "cadera")) == pytest.approx(A.Lt, abs=1e-6)
        xs, ys = r.joints["hombro"]
        xh, yh = r.joints["cadera"]
        assert math.hypot(xs - xh, ys - yh) == pytest.approx(A.Ltr, abs=1e-6)   # el tronco está en el plano sagital


def test_knee_3d_with_straight_toes_equals_planar_model():
    """Con knee_out = 0 la rodilla es la del modelo plano original (rodilla en el plano de la pierna que contiene al eje x)."""
    t = replace(TECHNIQUES["peso_muerto_convencional"], knee_out_deg=0.0)
    for s in S_GRID:
        r = solve(A, t, 100, s)
        psi = math.radians(r.psi_deg)
        xa, ya = r.joints["tobillo"]
        xh, yh = r.joints["cadera"]
        vh = (yh - A.ha) / math.cos(psi)
        ku, kv = _two_link(xa, 0.0, xh, vh, A.Ls, A.Lt)
        assert r.joints["rodilla"][0] == pytest.approx(ku, abs=1e-6)
        assert r.joints["rodilla"][1] == pytest.approx(A.ha + kv * math.cos(psi), abs=1e-6)


def test_knee_3d_keeps_hip_ankle_distance_and_branch_forward():
    ankle, hip = (-0.07, 0.08, 0.38), (-0.10, 0.70, 0.09)
    for out in (0.0, 20.0, 45.0, 60.0):
        k = _knee_3d(A, ankle, hip, out)
        assert _dist(k, ankle) == pytest.approx(A.Ls, abs=1e-9)
        assert _dist(k, hip) == pytest.approx(A.Lt, abs=1e-9)
    # más apertura -> la rodilla se desplaza hacia afuera y su proyección sagital queda más atrás
    k0, k45 = _knee_3d(A, ankle, hip, 0.0), _knee_3d(A, ankle, hip, 45.0)
    assert k45[2] > k0[2] and k45[0] < k0[0]


def test_deadlift_arm_length_and_bar_only_moves_forward_to_clear_the_knee():
    t = TECHNIQUES["peso_muerto_convencional"]
    shifted = 0
    for s in [i / 20 for i in range(21)]:
        r = solve(A, t, 150, s)
        xs, ys = r.joints["hombro"]
        assert math.hypot(r.bar[0] - xs, r.bar[1] - ys) == pytest.approx(A.La, abs=1e-6)
        assert r.bar[0] >= t.bar_x_offset - 1e-12            # solo se adelanta, nunca se atrasa
        shifted += r.bar[0] > t.bar_x_offset + 1e-9
    assert shifted < 21                                       # fuera del cruce con la rodilla la barra sigue la vertical
    assert solve(A, t, 150, 0.0).bar[0] == pytest.approx(t.bar_x_offset, abs=1e-12)
    assert solve(A, t, 150, 1.0).bar[0] == pytest.approx(t.bar_x_offset, abs=1e-12)


@pytest.mark.parametrize("key", ["peso_muerto_convencional", "peso_muerto_sumo"])
@pytest.mark.parametrize("athlete", ["SAMPLE", "SMALL", "TALL"])
def test_bar_never_through_the_knee(key, athlete):
    from dataclasses import replace as rep
    raws = {"SAMPLE": SAMPLE_MALE,
            "SMALL": rep(SAMPLE_MALE, H=1.62, h_troc=0.85, Lt_raw=0.40, Ls_raw=0.37, ha=0.07, Lf=0.24, Ltr_raw=0.47, La_raw=0.55),
            "TALL": rep(SAMPLE_MALE, H=1.92, h_troc=1.01, Lt_raw=0.48, Ls_raw=0.45, ha=0.09, Lf=0.29, Ltr_raw=0.55, La_raw=0.66)}
    a = from_raw(raws[athlete])
    t = TECHNIQUES[key]
    for i in range(41):
        try:
            r = solve(a, t, 130, i / 40)
        except Infeasible:
            continue
        d = math.hypot(r.bar[0] - r.joints["rodilla"][0], r.bar[1] - r.joints["rodilla"][1])
        assert d >= KNEE_CLEAR - 1e-6, (i / 40, d)


@pytest.mark.parametrize("key", ["peso_muerto_convencional", "peso_muerto_sumo"])
def test_bar_path_is_smooth_and_stays_out_of_the_shin(key):
    """El desvío de la barra alrededor de la rodilla no tiene esquinas (la circunferencia sola daba 54-67° respecto de la vertical;
    ahora ~44°, la inclinación del propio muslo) y la barra nunca queda dentro de la espinilla."""
    t = TECHNIQUES[key]
    prev, worst = None, 0.0
    for i in range(101):
        r = solve(A, t, 150, i / 100)
        d, _ = _seg_dist(r.bar[0], r.bar[1], r.joints["tobillo"], r.joints["rodilla"])
        assert d >= SHIN_FRONT + BAR_RADIUS - 1e-6, (i / 100, d)
        if prev is not None:
            worst = max(worst, abs(math.degrees(math.atan2(r.bar[0] - prev[0], r.bar[1] - prev[1]))))
        prev = r.bar
    assert worst < 50.0, worst


@pytest.mark.parametrize("key", ["peso_muerto_convencional", "peso_muerto_sumo"])
def test_deadlift_is_solvable_and_continuous_into_lockout(key):
    """Todo s en [0, 1] tiene solución (antes había un hueco en s ~ 0,985-0,994) y el tronco y la cadera llegan al bloqueo sin
    saltos: la solución con la rodilla clavada en su posición de pierna recta saltaba ~7° y ~6 cm al llegar a s = 1."""
    t = TECHNIQUES[key]
    prev, hip_top = None, None
    n = 400
    for i in range(n + 1):
        r = solve(A, t, 150, i / n)           # no debe lanzar Infeasible
        cur = (r.trunk_deg, r.joints["cadera"][0], r.joints["cadera"][1])
        if prev is not None:
            assert abs(cur[0] - prev[0]) < 5.5, (i / n, cur[0], prev[0])        # tronco: <= 4,5° por paso de 1/400
            assert abs(cur[1] - prev[1]) < 0.045, (i / n, cur[1], prev[1])      # cadera horizontal: <= 3,8 cm por paso
        prev = cur
    # la cadera nunca queda por encima de su altura al bloqueo
    top = solve(A, t, 150, 1.0).joints["cadera"][1]
    for i in range(n + 1):
        assert solve(A, t, 150, i / n).joints["cadera"][1] <= top + 1e-9


def test_sumo_knees_track_over_the_feet_and_come_in_as_the_legs_extend():
    t = TECHNIQUES["peso_muerto_sumo"]
    r0, r1 = solve(A, t, 100, 0.0), solve(A, t, 100, 1.0)
    za = r0.lateral["tobillo"]
    assert abs(r0.lateral["rodilla"] - za) < 0.10                       # rodilla en torno al pie al despegar (no 18 cm adentro)
    assert r0.lateral["rodilla"] - r1.lateral["rodilla"] > 0.05         # y se acerca al plano medio al extender (sobre todo al final)
    assert solve(A, t, 100, 0.75).lateral["rodilla"] > r1.lateral["rodilla"]
    assert r1.lateral["rodilla"] > r1.lateral["cadera"]                  # sin cruzar la cadera


@pytest.mark.parametrize("key", ["sentadilla_barra_alta", "sentadilla_barra_baja", "sentadilla_frontal"])
def test_squat_knees_do_not_collapse_inside_the_feet(key):
    r = solve(A, TECHNIQUES[key], 100, 1.0)
    assert abs(r.lateral["rodilla"] - r.lateral["tobillo"]) < 0.05      # pies rectos, rodilla sobre el pie en el fondo


def test_bar_only_hip_moment_identity():
    t = TECHNIQUES["peso_muerto_convencional"]
    r = solve(A, t, 150, 0.3)
    expected = G * 150 * (r.bar[0] - (r.joints["cadera"][0] + A.hip_ap_offset)) * math.cos(math.radians(r.psi_deg))
    assert r.bar_moments["cadera"] == pytest.approx(expected, rel=1e-12)
    assert r.arms["cadera"] == pytest.approx(r.bar[0] - (r.joints["cadera"][0] + A.hip_ap_offset))


def test_total_equals_bar_plus_body():
    r = solve(A, TECHNIQUES["sentadilla_barra_alta"], 100, 0.6)
    for k in r.moments:
        assert r.moments[k] == pytest.approx(r.bar_moments[k] + r.body_moments[k])


def test_ankle_moment_matches_ground_reaction():
    """Con pierna en el plano sagital (c = 0): M_tobillo = Fy*(x_CoP - x_a) - 2*m_pie*g*(x_pie - x_a)."""
    t = replace(TECHNIQUES["sentadilla_barra_alta"], w_frac=0.1, stance="frac")   # w < b  =>  c = 0
    r = solve(A, t, 100, 0.7)
    assert r.psi_deg == pytest.approx(0.0, abs=1e-9)
    mp = mass_points(A)
    mtot = A.mass + 100
    x_a = r.joints["tobillo"][0]
    x_foot = -A.Lf / 2 + mp["f_foot"] * A.Lf
    expected = G * mtot * (r.cop_x - x_a) - 2 * mp["m_foot"] * G * (x_foot - x_a) - 0.0
    # mtot incluye la masa de los pies; la identidad se cumple con la masa total
    assert r.moments["tobillo"] == pytest.approx(expected, rel=1e-9, abs=1e-6)


def test_de_leva_masses_sum_to_one():
    assert total_mass_check("M") == pytest.approx(1.0, abs=2e-3)
    assert total_mass_check("F") == pytest.approx(1.0, abs=2e-3)


def test_closure_checks_on_sample():
    c = closure_checks(SAMPLE_MALE)
    assert c["legs_ok"] and c["trunk_ok"]
    bad = replace(SAMPLE_MALE, h_troc=SAMPLE_MALE.h_troc + 0.05)
    assert not closure_checks(bad)["legs_ok"]


def test_sumo_vs_conventional_direction_matches_literature():
    """Frontiers 2025 (momentos pico fase 1, por pierna): cadera conv > sumo (303 vs 276); rodilla sumo > conv (80 vs 53)."""
    import calibrate as cal
    a = from_raw(cal.proxy_athlete(1.675, 76.0))
    peaks = {}
    for key in ("peso_muerto_convencional", "peso_muerto_sumo"):
        t = TECHNIQUES[key]
        skp = cal.knee_passing_s(a, t, 130.0)
        ms = [solve(a, t, 130.0, s).moments for s in np.linspace(0, skp, 20)]
        peaks[key] = (max(m["cadera"] for m in ms), max(m["rodilla"] for m in ms))
    conv, sumo = peaks["peso_muerto_convencional"], peaks["peso_muerto_sumo"]
    assert 1.0 < conv[0] / sumo[0] < 1.35
    assert sumo[1] > conv[1]


def test_sumo_has_shorter_bar_travel():
    """Escamilla 2000: recorrido vertical de la barra 47 cm (sumo) < 51 cm (convencional)."""
    a = from_raw(__import__("calibrate").proxy_athlete(1.675, 76.0))
    d = {}
    for key in ("peso_muerto_convencional", "peso_muerto_sumo"):
        t = TECHNIQUES[key]
        d[key] = solve(a, t, 130, 1.0).bar[1] - solve(a, t, 130, 0.0).bar[1]
    assert d["peso_muerto_sumo"] < d["peso_muerto_convencional"]


@pytest.mark.parametrize("key", ["peso_muerto_convencional", "peso_muerto_sumo", "sentadilla_barra_alta", "sentadilla_barra_baja"])
def test_calibration_reproduces_published_targets(key):
    """Cada objetivo publicado queda dentro de 2 DE; excepción conocida: la flexión de torso de la sentadilla (definición distinta)."""
    from calibrate_joint import predict_from_technique
    from calibrate import SQ_LABELS, DL_LABELS
    _, z = predict_from_technique(key)
    if key.startswith("sentadilla"):
        torso = SQ_LABELS.index("flexión de torso (°)")
        assert abs(z[torso]) < 6.0          # desajuste conocido y documentado (README)
        z = np.delete(z, torso)
    if key == "peso_muerto_sumo":
        # Excepción documentada (README, «Dirección de pies y rodillas»): los brazos de momento publicados para sumo de rodilla
        # (-13,2 y -4,8 cm) y de tobillo (-4,3 y -5,2 cm) ponen la barra DETRÁS del centro de la articulación, lo que no es
        # compatible con que la barra pase por delante de la rótula y de la tibia. Se prioriza la geometría físicamente posible;
        # esos cuatro objetivos pueden llegar a 3 DE y los demás siguen dentro de 2 DE.
        bar_arms = [DL_LABELS.index(n) for n in ("brazo rodilla LO", "brazo rodilla KP", "brazo tobillo LO", "brazo tobillo KP")]
        assert np.all(np.abs(z[bar_arms]) < 3.0), z
        z = np.delete(z, bar_arms)
    assert np.all(np.abs(z) < 2.0), z


def test_hips_rise_early_increases_hip_demand_midrange():
    base = TECHNIQUES["peso_muerto_convencional"]
    early = replace(base, p=0.6)
    assert solve(A, early, 150, 0.3).trunk_deg > solve(A, base, 150, 0.3).trunk_deg


def test_infeasible_when_arms_too_short():
    short = replace(A, La=0.2)
    with pytest.raises(Infeasible):
        solve(short, TECHNIQUES["peso_muerto_convencional"], 150, 0.0)


def test_monte_carlo_reproducible_and_wider_than_point():
    d1 = draw_perturbations(40, seed=7)
    d2 = draw_perturbations(40, seed=7)
    t = TECHNIQUES["peso_muerto_convencional"]
    a1 = monte_carlo(SAMPLE_MALE, t, 150, [0.0, 0.5], d1)
    a2 = monte_carlo(SAMPLE_MALE, t, 150, [0.0, 0.5], d2)
    assert np.allclose(a1, a2, equal_nan=True)
    assert np.nanstd(a1[:, 0, 2]) > 0


def test_shoulder_width_stance_sets_ankle_separation():
    """Postura estándar: la separación de tobillos es el ancho de hombros medido; pies rectos."""
    from model import _stance_c
    t = TECHNIQUES["sentadilla_barra_alta"]
    assert t.stance == "shoulder" and t.toe_out_deg == 0
    a = from_raw(SAMPLE_MALE)
    assert _stance_c(a, t) == pytest.approx((SAMPLE_MALE.D_acro - a.b) / 2)
    wider = from_raw(replace(SAMPLE_MALE, D_acro=SAMPLE_MALE.D_acro + 0.10))
    assert _stance_c(wider, t) == pytest.approx(_stance_c(a, t) + 0.05)
    # el sumo conserva su ancho de la literatura (no depende del ancho de hombros)
    sumo = TECHNIQUES["peso_muerto_sumo"]
    assert _stance_c(wider, sumo) == pytest.approx(_stance_c(a, sumo))
