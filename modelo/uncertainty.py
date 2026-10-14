"""Propagación de incertidumbre por Monte Carlo (presupuesto del documento del protocolo).

Tolerancias (± en m) -> desviación estándar = tol / 2 (la tolerancia se toma como ~95 %).
Son HIPÓTESIS de diseño pendientes de validar con el piloto. Cadera: error vertical de modo común (alarga el muslo y acorta
el tronco) y desfase antero-posterior (centro DEFAULT_HIP_AP_OFFSET inferido de la calibración, tolerancia ±4 cm).
"""
from dataclasses import replace
import math
import numpy as np
from anthropometry import RawMeasures, from_raw, DEFAULT_HIP_AP_OFFSET
from techniques import Technique
from model import solve, Infeasible

TOL = dict(Lt=0.010, Ls=0.010, ha=0.010, Ltr=0.010, La=0.010, D_asis=0.010, Lf=0.005, heel_to_malleolus=0.005,
           hip_vertical=0.015, shoulder=0.005, hip_ap=0.04, D_acro=0.010)
JOINTS = ("tobillo", "rodilla", "cadera", "lumbar")


def draw_perturbations(n: int, seed: int = 0, tol: dict = TOL):
    rng = np.random.default_rng(seed)
    sd = {k: v / 2.0 for k, v in tol.items()}
    return {k: rng.normal(0.0, sd[k], n) for k in sd}


def perturbed(raw: RawMeasures, d: dict, i: int):
    r = replace(raw,
                Lt_raw=raw.Lt_raw + d["Lt"][i], Ls_raw=raw.Ls_raw + d["Ls"][i], ha=raw.ha + d["ha"][i],
                Ltr_raw=raw.Ltr_raw + d["Ltr"][i] + d["shoulder"][i], La_raw=raw.La_raw + d["La"][i],
                D_asis=raw.D_asis + d["D_asis"][i], Lf=raw.Lf + d["Lf"][i],
                heel_to_malleolus=raw.heel_to_malleolus + d["heel_to_malleolus"][i],
                D_acro=(raw.D_acro + d["D_acro"][i]) if raw.D_acro > 0 else 0.0)
    return from_raw(r, hip_vertical_offset=d["hip_vertical"][i], hip_ap_offset=DEFAULT_HIP_AP_OFFSET + d["hip_ap"][i])


def monte_carlo(raw: RawMeasures, tech: Technique, load_kg: float, s_grid, draws: dict):
    """Devuelve un arreglo [n, len(s_grid), 4] de momentos (N·m); NaN donde la postura es inalcanzable."""
    n = len(draws["Lt"])
    out = np.full((n, len(s_grid), len(JOINTS)), np.nan)
    for i in range(n):
        a = perturbed(raw, draws, i)
        for j, s in enumerate(s_grid):
            try:
                m = solve(a, tech, load_kg, s).moments
            except Infeasible:
                continue
            out[i, j, :] = [m[k] for k in JOINTS]
    return out


def band(arr: np.ndarray, q=(2.5, 50, 97.5)):
    """Percentiles por (s, articulación) ignorando NaN."""
    return np.nanpercentile(arr, q, axis=0)   # [3, len(s), 4]
