"""Fichas de técnica (parámetros de estilo).

Estado de calibración (ver calibrate_joint.py y README):
 - peso muerto convencional/sumo y sentadilla barra alta/baja: AJUSTADAS a Escamilla y cols. (2000) y Swinton y cols. (2012),
   con un desfase antero-posterior de cadera compartido (anthropometry.DEFAULT_HIP_AP_OFFSET) inferido del mismo ajuste.
 - sentadilla frontal: SIN datos; valores ilustrativos (calibrated=False).
 - `ranges` = percentiles 2,5-97,5 de repetir el ajuste con los objetivos publicados remuestreados dentro de su desviación estándar.
 - `at_bound` = parámetros que quedaron en el borde del rango de búsqueda: el dato los empuja más allá de lo plausible; leerlos con cautela.
"""
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Technique:
    name: str
    exercise: str            # "squat" | "deadlift"
    w_frac: float            # distancia entre centros de tobillos / estatura
    toe_out_deg: float       # rotación externa del pie respecto a la dirección de avance (reduce la base de apoyo antero-posterior)
    calibrated: bool = False
    stance: str = "frac"     # "shoulder": tobillos bajo los hombros (ancho = Anthro.shoulder_w; si no hay medida, w_frac) | "frac": w_frac * estatura
    source: str = ""
    ranges: dict = field(default_factory=dict)
    at_bound: tuple = ()
    # --- sentadilla ---
    lam: float = 1.0         # posición de la barra a lo largo del tronco (0 = cadera, 1 = hombro)
    d: float = 0.06          # desplazamiento de la barra hacia atrás, perpendicular al tronco (m); negativo = frontal
    shank_top_deg: float = 3.0
    shank_max_deg: float = 35.0   # inclinación de tibia en el fondo (dorsiflexión disponible)
    depth_deg: float = 5.0        # muslo bajo la horizontal en el fondo (0 = paralelo)
    bar_x_mode: str = "midfoot"    # "midfoot": barra vertical sobre bar_x_offset | "erect_top": vertical con tronco erguido arriba
    trunk_top_deg: float = 0.0
    bar_x_offset: float = 0.0     # posición horizontal de la barra respecto al medio pie (m, + = adelante)
    # --- peso muerto ---
    f_h0: float = 0.55       # altura inicial de cadera sobre el tobillo / altura final
    eps0: float = 0.05       # hombro por delante de la barra al inicio (m)
    eps_top: float = -0.03   # hombro por detrás de la barra arriba (m)
    p: float = 1.0           # sincronía: <1 cadera sube antes, >1 cadera se retrasa
    plate_radius: float = 0.225
    # --- rodilla en 3D (ver model._knee_3d) ---
    knee_out_deg: float = 0.0     # dirección horizontal en que se abre la rodilla, hacia afuera respecto al eje x (0 = hacia adelante).
                                  # Rodilla sobre el pie: ver README, «Dirección de pies y rodillas».
    grip_inside: bool = False     # peso muerto: agarre por dentro de las piernas (sumo)


_ESC = "Escamilla y cols. (2000), 40 powerlifters hombres (~167 cm, 76 kg); stance y ángulo de pie publicados"
_SWI = "Swinton y cols. (2012), 12 powerlifters hombres (180,3 cm, 100,2 kg), 70 % 1RM, inicio de la fase concéntrica"

TECHNIQUES = {
    "sentadilla_barra_alta": Technique(
        "Sentadilla barra alta (squat 'tradicional')", "squat", w_frac=48.3 / 180.3, toe_out_deg=0, stance="shoulder", calibrated=True,
        source=_SWI + ". toe_out, knee_out, lam y d: sin fuente (ilustrativos).",
        knee_out_deg=12.0, lam=1.0, d=0.06, shank_max_deg=32.8, depth_deg=-1.9, bar_x_offset=0.040,
        ranges=dict(shank_max_deg=(28.8, 36.5), depth_deg=(-7.6, 6.6), bar_x_offset=(0.015, 0.060))),
    "sentadilla_barra_baja": Technique(
        "Sentadilla barra baja (squat 'powerlifting')", "squat", w_frac=89.6 / 180.3, toe_out_deg=0, stance="shoulder", calibrated=True,
        source=_SWI + ". toe_out, knee_out, lam y d: sin fuente (ilustrativos).",
        knee_out_deg=12.0, lam=0.93, d=0.07, shank_max_deg=20.8, depth_deg=2.7, bar_x_offset=-0.021,
        ranges=dict(shank_max_deg=(16.9, 25.6), depth_deg=(-15.0, 8.7), bar_x_offset=(-0.044, 0.011))),
    "sentadilla_frontal": Technique(
        "Sentadilla frontal", "squat", w_frac=0.30, toe_out_deg=0, stance="shoulder", calibrated=False,
        source="Sin datos: valores ilustrativos.",
        knee_out_deg=12.0, lam=1.0, d=-0.05, shank_max_deg=40, depth_deg=5, bar_x_offset=0.0),
    "peso_muerto_convencional": Technique(
        "Peso muerto convencional", "deadlift", w_frac=37 / 167, toe_out_deg=15, stance="shoulder", calibrated=True, source=_ESC,
        # barra vertical sobre el medio pie (bar_x_offset = 0, FIJO por criterio de entrenador; ajuste previo con la barra libre:
        # -3,2 cm [-5,3, +1,7]); los otros cuatro parámetros se reajustaron con la barra ahí y la rodilla 3D.
        knee_out_deg=15.0, f_h0=0.8424, eps0=0.12, eps_top=0.0177, p=1.5633, bar_x_offset=0.0,
        at_bound=("eps0",),
        ranges=dict(f_h0=(0.794, 0.853), eps0=(0.12, 0.12), eps_top=(0.012, 0.02), p=(1.127, 1.573))),
    "peso_muerto_sumo": Technique(
        "Peso muerto sumo", "deadlift", w_frac=73 / 168, toe_out_deg=45, calibrated=True, source=_ESC,
        # barra vertical sobre el medio pie (bar_x_offset = 0, FIJO por criterio de entrenador; el ajuste previo con la barra libre
        # daba -7,5 cm, sobre el tobillo, que dejaba la barra dentro de la tibia y la rodilla).
        knee_out_deg=25.0, grip_inside=True, f_h0=0.7762, eps0=0.102, eps_top=-0.0153, p=1.2119, bar_x_offset=0.0,
        at_bound=(),
        ranges=dict(f_h0=(0.776, 0.827), eps0=(0.102, 0.107), eps_top=(-0.015, 0.010), p=(1.212, 1.440))),
}
