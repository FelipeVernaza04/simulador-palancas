import math
from dataclasses import replace
import numpy as np
import pytest

from anthropometry import (SAMPLE_MALE, from_raw, closure_checks, total_mass_check, mass_points, G, DE_LEVA)
from techniques import TECHNIQUES
from model import solve, Infeasible
from uncertainty import draw_perturbations, monte_carlo

A = from_raw(SAMPLE_MALE)
S_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]


def _v(a, y, psi):
    return (y - a.ha) / math.cos(psi)


@pytest.mark.parametrize("key", list(TECHNIQUES))
def test_segment_lengths_preserved(key):
    t = TECHNIQUES[key]
    for s in S_GRID:
        r = solve(A, t, 100, s)
        psi = math.radians(r.psi_deg)
        xa, ya = r.joints["tobillo"]
        xk, yk = r.joints["rodilla"]
        xh, yh = r.joints["cadera"]
        # las piernas viven en el plano inclinado: longitudes reales con v = (y - ha)/cos(psi)
        assert math.hypot(xk - xa, _v(A, yk, psi)) == pytest.approx(A.Ls, abs=1e-6)
        assert math.hypot(xh - xk, _v(A, yh, psi) - _v(A, yk, psi)) == pytest.approx(A.Lt, abs=1e-6)
        xs, ys = r.joints["hombro"]
        assert math.hypot(xs - xh, ys - yh) == pytest.approx(A.Ltr, abs=1e-6)   # el tronco está en el plano sagital


def test_deadlift_arm_length_and_vertical_bar():
    t = TECHNIQUES["peso_muerto_convencional"]
    for s in S_GRID:
        r = solve(A, t, 150, s)
        xs, ys = r.joints["hombro"]
        assert math.hypot(r.bar[0] - xs, r.bar[1] - ys) == pytest.approx(A.La, abs=1e-6)
        assert r.bar[0] == pytest.approx(t.bar_x_offset, abs=1e-12)


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
    from calibrate import SQ_LABELS
    _, z = predict_from_technique(key)
    if key.startswith("sentadilla"):
        torso = SQ_LABELS.index("flexión de torso (°)")
        assert abs(z[torso]) < 6.0          # desajuste conocido y documentado (README)
        z = np.delete(z, torso)
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
