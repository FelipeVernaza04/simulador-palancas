# Prototipo de cálculo: palancas en sentadilla y peso muerto

Modelo sagital 2D, cuasi-estático, simétrico. Sin interfaz. Python 3 + numpy + scipy (pytest para las pruebas).

    python3 -m pytest -q                     # 20 pruebas
    python3 run_demo.py                      # bandas de incertidumbre y comparación pareada sumo vs convencional
    python3 calibrate_joint.py fixbar        # reproduce la calibración (tarda unos minutos)

## Archivos
- `anthropometry.py`: 11 medidas del protocolo -> parámetros corregidos a centros articulares; masas de de Leva (1996); cierres geométricos; `DEFAULT_HIP_AP_OFFSET` (+5,8 cm, inferido).
- `techniques.py`: fichas de técnica. Cuatro calibradas con literatura (con rangos); sentadilla frontal ilustrativa.
- `model.py`: cinemática (plano de pierna inclinado para postura ancha), momentos tobillo/rodilla/cadera/lumbar, CoP.
- `uncertainty.py`: Monte Carlo con el presupuesto de incertidumbre (incluye el desfase de cadera).
- `calibrate.py`, `calibrate_joint.py`: tablas de objetivos publicados y ajuste conjunto por mínimos cuadrados ponderados.
- `calibration_result_fixbar.json`: parámetros ajustados y rangos (P2,5–P97,5).
- `test_model.py`: pruebas (geometría, identidades de momento, equilibrio, direcciones sumo/convencional, reproducción de objetivos publicados).

## Convenciones
x adelante, y arriba, origen bajo el medio pie. Momentos extensores positivos (N·m), carga total de ambas piernas.
Slider `s` ∈ [0, 1]: sentadilla de pie -> fondo; peso muerto barra en el suelo -> bloqueo.

## Interfaz web (carpeta `app/`)
`app/model.js` es el modelo portado a JavaScript (se comprobó contra Python en 360 casos, diferencia relativa máxima 4·10⁻⁷). `app/app.src.html` es la interfaz del modo usuario y `app/build.py` genera `index.html` con el modelo incrustado. Modo investigador pendiente.
La pestaña «Qué pasaría si» cambia tres parámetros de la técnica (peso muerto: posición de la barra, altura de la cadera, distancia entre tobillos; sentadilla: distancia entre tobillos, profundidad, inclinación de la tibia) y muestra el efecto en el brazo de palanca y el pico de momento. Los rangos de los controles son rangos prácticos, no datos publicados (la distancia entre tobillos se limita a lo observado entre una técnica y su variante). Los textos de las indicaciones son borradores para validar con un entrenador.

## Figura: piernas a longitud real (cambio de dibujo, no de cálculo)
Con postura ancha las piernas yacen en un plano inclinado el ángulo ψ respecto al plano sagital. En una vista lateral pura la tibia y el muslo se ven más cortos (factor cos ψ): en sentadilla barra baja la tibia llegaba a verse al 54 % de su longitud en el fondo (ψ = 64,5°), en sumo entre 90 % y 94 %, en convencional ≈ 99 %. El tronco siempre está en el plano sagital y no se acorta.
La figura ofrece ahora dos vistas (selector bajo el dibujo): «Piernas reales» (por defecto) abre el plano de la pierna sobre el sagital, de modo que tibia, muslo y tronco miden lo mismo en todo el recorrido (comprobado: 1,00 ± 0,002 en las cinco técnicas); la x de cada articulación no cambia, así que los brazos de palanca y los momentos son los mismos. Como consecuencia, cadera, tronco, brazos y barra se dibujan más arriba que su altura real (hasta 22 cm en el fondo de la barra baja, 7 cm en sumo); un texto bajo la figura lo indica. En peso muerto el disco se dibuja apoyado en el suelo. «Vista lateral» muestra la proyección real, con alturas reales y piernas acortadas.
## Postura estándar: tobillos bajo los hombros, pies rectos
Tras ver ψ = 64,5° en el fondo de la barra baja (ancho de Swinton, 89,6 cm, de definición no verificada), la postura de sentadilla (alta, baja, frontal) y de peso muerto convencional se estandarizó: **separación de centros de tobillo = ancho de hombros medido (acromion–acromion, medida `Wa`/`D_acro`, nueva) y pies rectos (`toe_out_deg` = 0)**. El sumo conserva su ancho de la literatura (73/168 de la estatura) y 25° de rotación. Si no se mide `Wa` (solo posible en sumo, donde es opcional) la app usa 0,23·estatura (estimación propia, sin fuente verificada).
Efecto: ψ máximo en sentadilla baja de 64,5° a ≈ 17° (tibia en vista lateral al 96 % en el fondo, no 54 %). Las magnitudes horizontales (brazos de palanca, ángulos en el plano, tronco) no dependen de la separación lateral, así que **la calibración no cambia** (se repitió: mismos valores, rangos algo distintos). Lo que sí cambia es el factor cos ψ de los momentos de pierna: en barra baja en el fondo los momentos de pierna eran ≈ 55 % menores que ahora por ese factor (artefacto del ancho anterior).
Los valores de ancho de los estudios (48,3, 89,6 cm; stance 37/167) quedan en `w_frac` solo como respaldo cuando no hay `D_acro`. El proxy de calibración usa `D_acro` = 0,23·H (asunción).

## Modo explicar (para que el entrenador enseñe en vivo)
Botón «Explicar con el atleta» en Resultados. Pantalla pensada para celular o tablet: figura grande, tres momentos clave (sentadilla: de pie, mitad, fondo; peso muerto: despegue, la barra pasa la rodilla, bloqueo), botón de reproducir (con «prefiere menos movimiento» salta de momento en momento) y control de recorrido. En cada momento muestra una sola idea: el brazo de palanca de la articulación con más demanda, resaltado en la figura y dicho en una frase generada con los números del atleta, y «Dónde se siente la carga» con barras por zona (glúteos e isquios, cuádriceps, pantorrillas, espalda baja).
Criterios: la zona principal se elige entre cadera, rodilla y tobillo (momentos extensores positivos); la espalda baja se muestra como referencia aproximada y no entra en la elección, porque su punto de cálculo es aproximado y con el desfase de cadera siempre sale mayor que la cadera. Las barras usan una escala común (el mayor momento de pierna de todo el recorrido), así crecen y bajan al mover el recorrido. Si el mayor momento de pierna es menor al 25 % del pico, el texto dice que casi no hay palanca (cuerpo apilado, carga en estructuras pasivas). Textos son borradores para validar con un entrenador. No incluye indicaciones (cues) en pantalla ni comparación lado a lado (se descartaron por ahora).

## Cambio en la sentadilla tras la calibración
La pierna se parametriza ahora por el ángulo del muslo (lineal en s) y no por la altura de la cadera. Los extremos (de pie y fondo) no cambian, así que la calibración es la misma; solo cambia el interior del recorrido, donde antes la cadera saltaba hacia atrás en el primer 5 %.

## Calibración (método)
Fuentes: Escamilla y cols. (2000) — peso muerto convencional y sumo, 40 powerlifters hombres (stance, ángulos de rodilla/cadera/tronco, brazos de momento, recorrido de la barra); Swinton y cols. (2012) — sentadilla barra alta y baja, 12 powerlifters hombres al 70 % 1RM; Frontiers (2025) — momentos pico por pierna, usados solo como razones/diferencias de validación.
Se construyó un atleta «proxy» por estudio (estatura y masa publicadas, proporciones de la muestra de referencia). Las cuatro técnicas se ajustaron **a la vez** minimizando χ² entre modelo y objetivos publicados (peso = 1/DE), con:
- un desfase antero-posterior de cadera `hx` **compartido** (a priori N(0, 3 cm));
- una penalización blanda para que el CoP quede dentro del 85 % central de la base de apoyo;
- rangos por remuestreo de los objetivos dentro de su DE (30 reajustes).
Resultado: χ² = 55,1 (31 objetivos + 10 penalizaciones + 1 a priori; 21 parámetros).

## Parámetros calibrados (valor [rango 95 %])
| Técnica | Parámetros |
|---|---|
| Peso muerto convencional | f_h0 0,814 [0,745–0,842]; barra x −3,2 cm [−5,3, +1,7]; eps0, eps_top, p **en el límite** |
| Peso muerto sumo | f_h0 0,762 [0,635–0,794]; barra x −7,5 cm [−8,2, −6,7]; eps0, eps_top, p **en el límite** |
| Sentadilla barra alta | tibia máx. 32,8° [28,8–36,5]; profundidad −1,9° [−7,6, 6,6]; barra x +4,0 cm [1,5, 6,0] |
| Sentadilla barra baja | tibia máx. 20,8° [16,9–25,6]; profundidad +2,7° [−15, 8,7]; barra x −2,1 cm [−4,4, +1,1] |
| Desfase de cadera hx | **+5,8 cm [3,2, 8,4]** (inferido) |

## Validación independiente (no ajustada)
- Ángulo del tronco sobre la horizontal al despegue: convencional 6,0° (lit. 11 ± 7), sumo 15,2° (lit. 16 ± 6).
- Momento de cadera pico conv/sumo: 1,15 (lit. 1,10). Diferencia por pierna en cadera conv−sumo: 31 N·m (lit. +27).
- Rodilla sumo−conv por pierna: 35 N·m (lit. +27). El recorrido de barra del sumo es más corto que el convencional (como en Escamilla: 47 vs 51 cm).

## Limitaciones conocidas (leer antes de usar)
1. **`hx` = +5,8 cm es una inferencia**, no una medida: el ajuste lo necesita para igualar los brazos de momento de cadera publicados y está confundido con diferencias de definición del centro de cadera. Anatómicamente plausible pero mayor de lo esperado. Pendiente de verificar con imagen (piloto).
2. **`eps0`, `eps_top`, `p` quedan en el borde del rango de búsqueda** en ambos pesos muertos (hombros ~12 cm por delante de la barra al inicio; cadera que se retrasa). Los datos vienen de levantadores de menor nivel con tronco casi horizontal; no extrapolar a levantadores avanzados.
3. Sentadilla: la inclinación del tronco queda **sobreestimada** (barra alta 45° vs 33,5 ± 4,6; barra baja 54° vs 33,1 ± 4,5). La posición de la barra en el tronco (`lam`, `d`) no es identificable con estos datos: se fijó a valores ilustrativos. Prueba con tolerancia explícita.
4. Brazo de momento de rodilla: sumo LO 0 vs −13,2 ± 8 (z +1,6); barra baja z −1,9.
5. Definiciones de ángulo (tronco respecto a la horizontal, muslo, ancho de postura ≈ distancia entre tobillos, instante «inicio de fase concéntrica») fueron **inferidas** de los textos; verificar con los artículos originales.
6. Rotación externa del pie (`toe_out_deg`) y sentadilla frontal: sin fuente.
7. No usados por tablas ilegibles o definiciones ambiguas: Glassbrook y cols. (2017) y otros estudios de sentadilla. Fuentes de PubMed/PMC bloqueadas por reCAPTCHA (no se eludió): ISAK, ANSUR II, Wretenberg, Escamilla 2001 (sentadilla) no verificadas.
8. Cuasi-estático: subestima picos; los momentos de la literatura incluyen dinámica y son por pierna, solo se compararon razones/diferencias. Momentos de tobillo no comparados.
9. Lumbar: punto a 10 % de Ltr sobre la cadera (aproximación). Cabeza: ~6,9 % de H sobre el hombro (aproximación propia).
10. Los números de de Leva se transcribieron de un PDF: verificar contra la tabla original.
11. Sumo solo en el plano sagital; sin demanda frontal (abductores, rotadores).
