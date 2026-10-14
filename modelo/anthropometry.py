"""Antropometría: medidas crudas -> parámetros del modelo, y parámetros de segmentos (de Leva 1996).

Los números de de Leva se transcribieron de texto extraído de un PDF: VERIFICAR contra la tabla original.
Unidades: metros, kilogramos, radianes.
"""
from dataclasses import dataclass, replace
import math

G = 9.81

# --- de Leva (1996), Zatsiorsky-Seluyanov ajustado --------------------------------------------
# Masa en fracción de la masa corporal; com en fracción desde el extremo proximal (pie: desde el talón;
# tronco: desde la escotadura supraesternal hacia las caderas). Longitudes (m) de los segmentos del brazo.
DE_LEVA = {
    "M": dict(thigh=(0.1416, 0.4095), shank=(0.0433, 0.4459), foot=(0.0137, 0.4415),
              head=0.0694, trunk=(0.4346, 0.4486),
              upper_arm=(0.0271, 0.5772, 0.2817), forearm=(0.0162, 0.4574, 0.2689), hand=(0.0061, 0.7900, 0.0862)),
    "F": dict(thigh=(0.1478, 0.3612), shank=(0.0481, 0.4416), foot=(0.0129, 0.4014),
              head=0.0668, trunk=(0.4257, 0.4151),
              upper_arm=(0.0255, 0.5754, 0.2751), forearm=(0.0138, 0.4559, 0.2643), hand=(0.0056, 0.7474, 0.0780)),
}
# Centro glenohumeral distal al acromion (de Leva, tabla 2), m
SHOULDER_BELOW_ACROMION = {"M": 0.0345, "F": 0.0337}
# Desfase antero-posterior del centro de cadera respecto al vértice trocantérico del modelo (m, + = más anterior).
# INFERIDO del ajuste conjunto a Escamilla 2000 y Swinton 2012 (calibrate_joint.py), no medido; rango P2.5-P97.5: 0.032 a 0.084 m.
DEFAULT_HIP_AP_OFFSET = 0.058
# Cabeza: desplazamiento del centro de masa sobre el centro del hombro, como fracción de H (APROXIMACIÓN propia)
HEAD_COM_ABOVE_SHOULDER_FRAC_H = 0.12 / 1.75


def arm_mass_and_com_fraction(sex: str):
    """Masa de ambos brazos (fracción) y posición de su centro de masa como fracción de La (hombro -> centro del agarre)."""
    p = DE_LEVA[sex]
    m_ua, c_ua, L_ua = p["upper_arm"]
    m_fa, c_fa, L_fa = p["forearm"]
    m_h, c_h, L_h = p["hand"]
    pos = [c_ua * L_ua, L_ua + c_fa * L_fa, L_ua + L_fa + c_h * L_h]
    masses = [m_ua, m_fa, m_h]
    com = sum(m * x for m, x in zip(masses, pos)) / sum(masses)
    L_grip = L_ua + L_fa + 0.5 * L_h
    return 2 * sum(masses), com / L_grip


@dataclass(frozen=True)
class RawMeasures:
    """Las 11 medidas del protocolo (las dos de control se usan solo en `closure_checks`)."""
    sex: str                 # "M" | "F"
    mass: float              # kg
    H: float                 # estatura
    h_troc: float            # altura del trocánter (control)
    Lt_raw: float            # trocánter -> epicóndilo lateral
    Ls_raw: float            # epicóndilo -> maléolo lateral
    ha: float                # altura del maléolo lateral
    Lf: float                # longitud del pie
    heel_to_malleolus: float # talón -> eje del maléolo (horizontal)
    Ltr_raw: float           # trocánter -> acromion
    h_acro: float            # altura del acromion (control)
    La_raw: float            # acromion -> centro del agarre
    D_asis: float            # distancia entre EIAS
    D_acro: float = 0.0      # distancia acromion-acromion (0 = no medida); fija el ancho de postura estándar (tobillos bajo los hombros)


@dataclass(frozen=True)
class Anthro:
    """Parámetros del modelo, ya corregidos a centros articulares."""
    sex: str
    mass: float
    H: float
    Lt: float
    Ls: float
    ha: float
    Lf: float
    dam: float   # distancia tobillo -> medio pie (horizontal, positiva hacia adelante)
    Ltr: float   # cadera -> centro del hombro
    La: float    # centro del hombro -> centro del agarre
    b: float     # distancia entre centros de cadera
    hip_ap_offset: float = 0.0   # centro de cadera respecto al vértice del modelo (m, + = más anterior)
    shoulder_w: float = 0.0      # ancho de hombros (acromion-acromion): distancia entre tobillos en la postura estándar


def from_raw(r: RawMeasures, hip_vertical_offset: float = 0.0, hip_ap_offset=None) -> Anthro:
    """Aplica las correcciones del protocolo.

    hip_vertical_offset (m, + = centro de cadera por encima de la medida del trocánter): desplaza el centro de
    cadera. Alarga/acorta el muslo y acorta/alarga el tronco en sentido opuesto (error de modo común).
    """
    sh = SHOULDER_BELOW_ACROMION[r.sex]
    if hip_ap_offset is None:
        hip_ap_offset = DEFAULT_HIP_AP_OFFSET
    return Anthro(
        sex=r.sex, mass=r.mass, H=r.H,
        Lt=r.Lt_raw - hip_vertical_offset,
        Ls=r.Ls_raw, ha=r.ha, Lf=r.Lf,
        dam=r.Lf / 2.0 - r.heel_to_malleolus,
        Ltr=r.Ltr_raw - sh + hip_vertical_offset,
        La=r.La_raw - sh,
        b=0.66 * r.D_asis + 0.0146,   # Harrington 2007, ec. 1 (transcrita de HAS-Motion)
        hip_ap_offset=hip_ap_offset,
        shoulder_w=r.D_acro,
    )


def closure_checks(r: RawMeasures) -> dict:
    """Cierres geométricos del protocolo (tolerancias propuestas: 2 cm y 3 cm)."""
    legs = r.h_troc - (r.ha + r.Ls_raw + r.Lt_raw)
    trunk = r.h_acro - (r.h_troc + r.Ltr_raw)
    return {"legs_cm": legs * 100, "trunk_cm": trunk * 100,
            "legs_ok": abs(legs) <= 0.02, "trunk_ok": abs(trunk) <= 0.03}


def mass_points(a: Anthro):
    """Masas absolutas (kg) de los segmentos del modelo."""
    p = DE_LEVA[a.sex]
    arms_frac, arm_com = arm_mass_and_com_fraction(a.sex)
    return dict(
        m_thigh=p["thigh"][0] * a.mass, f_thigh=p["thigh"][1],
        m_shank=p["shank"][0] * a.mass, f_shank=p["shank"][1],
        m_foot=p["foot"][0] * a.mass, f_foot=p["foot"][1],
        m_head=p["head"] * a.mass,
        m_trunk=p["trunk"][0] * a.mass, f_trunk_from_hip=1.0 - p["trunk"][1],
        m_arms=arms_frac * a.mass, f_arm_com=arm_com,
        head_above_shoulder=HEAD_COM_ABOVE_SHOULDER_FRAC_H * a.H,
    )


def total_mass_check(sex: str) -> float:
    p = DE_LEVA[sex]
    arms, _ = arm_mass_and_com_fraction(sex)
    return p["head"] + p["trunk"][0] + arms + 2 * (p["thigh"][0] + p["shank"][0] + p["foot"][0])


SAMPLE_MALE = RawMeasures(sex="M", mass=80.0, H=1.75, h_troc=0.92, Lt_raw=0.42, Ls_raw=0.43, ha=0.07, Lf=0.26,
                          heel_to_malleolus=0.05, Ltr_raw=0.52, h_acro=1.44, La_raw=0.62, D_asis=0.25, D_acro=0.40)
