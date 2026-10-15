# Prototipo de cálculo: palancas en sentadilla y peso muerto

Modelo sagital 2D, cuasi-estático, simétrico. Sin interfaz. Python 3 + numpy + scipy (pytest para las pruebas).

    python3 -m pytest -q                     # 33 pruebas
    python3 run_demo.py                      # bandas de incertidumbre y comparación pareada sumo vs convencional
    python3 calibrate_joint.py fixbar        # reproduce la calibración (tarda unos minutos)

## Archivos
- `anthropometry.py`: 11 medidas del protocolo -> parámetros corregidos a centros articulares; masas de de Leva (1996); cierres geométricos; `DEFAULT_HIP_AP_OFFSET` (+5,8 cm, inferido).
- `techniques.py`: fichas de técnica. Cuatro calibradas con literatura (con rangos); sentadilla frontal ilustrativa.
- `model.py`: cinemática (plano de pierna inclinado para postura ancha; rodilla en 3D según la dirección del pie), momentos tobillo/rodilla/cadera/lumbar, CoP.
- `uncertainty.py`: Monte Carlo con el presupuesto de incertidumbre (incluye el desfase de cadera).
- `calibrate.py`, `calibrate_joint.py`: tablas de objetivos publicados y ajuste conjunto por mínimos cuadrados ponderados.
- `calibration_result_fixbar.json`: parámetros ajustados y rangos (P2,5–P97,5).
- `test_model.py`: pruebas (geometría, identidades de momento, equilibrio, direcciones sumo/convencional, reproducción de objetivos publicados).

## Convenciones
x adelante, y arriba, origen bajo el medio pie. Momentos extensores positivos (N·m), carga total de ambas piernas.
Slider `s` ∈ [0, 1]: sentadilla de pie -> fondo; peso muerto barra en el suelo -> bloqueo.

## Interfaz web (carpeta `app/`)
`app/model.js` es el modelo portado a JavaScript (se comprobó contra Python tras añadir la rodilla 3D: 420 puntos de dos atletas y las cinco técnicas, diferencia máxima 3·10⁻¹² fuera del bloqueo y ≤ 10⁻⁴ N·m en el bloqueo, donde la pierna recta hace numéricamente singular la rodilla). `app/app.src.html` es la interfaz del modo usuario y `app/build.py` genera `index.html` con el modelo incrustado. Modo investigador pendiente.
La pestaña «Qué pasaría si» cambia tres parámetros de la técnica (peso muerto: posición de la barra, altura de la cadera, distancia entre tobillos; sentadilla: distancia entre tobillos, profundidad, inclinación de la tibia) y muestra el efecto en el brazo de palanca y el pico de momento. Los rangos de los controles son rangos prácticos, no datos publicados (la distancia entre tobillos se limita a lo observado entre una técnica y su variante). Los textos de las indicaciones son borradores para validar con un entrenador.

## Figura: piernas a longitud real (cambio de dibujo, no de cálculo)
Con postura ancha las piernas yacen en un plano inclinado el ángulo ψ respecto al plano sagital. En una vista lateral pura la tibia y el muslo se ven más cortos (factor cos ψ): en sentadilla barra baja la tibia llegaba a verse al 54 % de su longitud en el fondo (ψ = 64,5°), en sumo entre 90 % y 94 %, en convencional ≈ 99 %. El tronco siempre está en el plano sagital y no se acorta.
La figura ofrece dos vistas (selector bajo el dibujo): «Vista lateral» (por defecto) muestra la proyección real, con alturas reales y piernas acortadas; «Piernas reales» gira cada punto de la pierna alrededor del eje antero-posterior del tobillo hasta ponerlo en el plano sagital, de modo que la tibia mide lo mismo en todo el recorrido; la x de cada articulación no cambia, así que los brazos de palanca y los momentos son los mismos. Como consecuencia, cadera, tronco, brazos y barra se dibujan más arriba que su altura real (hasta 22 cm en el fondo de la barra baja, 7 cm en sumo) y el muslo ya no queda siempre a su longitud real (96 % en sumo, porque la rodilla se abre hacia afuera); un texto bajo la figura lo indica. En peso muerto el disco se dibuja apoyado en el suelo.
**Por qué la vista lateral es ahora la predeterminada:** en «Piernas reales» la cadera y la barra suben más que la rodilla, lo que cambia en qué punto del recorrido se cruzan y deja la distancia barra-rodilla dibujada en 1,6–2,0 cm en sumo, aunque la real sea ≥ 6,4 cm. La vista lateral conserva esa distancia. La distancia mínima dibujada se midió en las dos vistas y dos atletas.
## Postura estándar: tobillos bajo los hombros, pies rectos
Tras ver ψ = 64,5° en el fondo de la barra baja (ancho de Swinton, 89,6 cm, de definición no verificada), la postura de sentadilla (alta, baja, frontal) y de peso muerto convencional se estandarizó: **separación de centros de tobillo = ancho de hombros medido (acromion–acromion, medida `Wa`/`D_acro`, nueva) y pies rectos (`toe_out_deg` = 0)**. El sumo conserva su ancho de la literatura (73/168 de la estatura). (La rotación del pie en sumo se cambió después de 25° a 45°: ver «Dirección de pies y rodillas».) Si no se mide `Wa` (solo posible en sumo, donde es opcional) la app usa 0,23·estatura (estimación propia, sin fuente verificada).
Efecto: ψ máximo en sentadilla baja de 64,5° a ≈ 17° (tibia en vista lateral al 96 % en el fondo, no 54 %). Las magnitudes horizontales (brazos de palanca, ángulos en el plano, tronco) no dependen de la separación lateral, así que **la calibración no cambia** (se repitió: mismos valores, rangos algo distintos). Lo que sí cambia es el factor cos ψ de los momentos de pierna: en barra baja en el fondo los momentos de pierna eran ≈ 55 % menores que ahora por ese factor (artefacto del ancho anterior).
Los valores de ancho de los estudios (48,3, 89,6 cm; stance 37/167) quedan en `w_frac` solo como respaldo cuando no hay `D_acro`. El proxy de calibración usa `D_acro` = 0,23·H (asunción).

## Dirección de pies y rodillas (modelo 3D de la rodilla)
**Qué cambió.** Antes la rodilla se resolvía dentro del plano de la pierna que contiene al eje x (pies rectos, cualquiera fuera la técnica). Ahora la tibia y el muslo (longitudes reales) dejan una circunferencia de posiciones posibles de la rodilla alrededor del eje tobillo-cadera, y la rodilla se coloca en el plano que contiene ese eje y la dirección horizontal `knee_out_deg` (hacia afuera respecto al eje x): la rodilla sigue la dirección del pie. Con `knee_out_deg = 0` es exactamente el modelo anterior (prueba `test_knee_3d_with_straight_toes_equals_planar_model`). La cadera, el tobillo y el ángulo de rodilla no cambian; solo gira la rodilla sobre el eje. El modelo devuelve además la distancia lateral de tobillo, rodilla y cadera (`Result.lateral`), que usa la vista frontal.

**Valores usados (sin fuente; criterio de entrenador y coherencia geométrica):**
| Técnica | `toe_out_deg` (pie) | `knee_out_deg` (rodilla) |
|---|---|---|
| Sentadilla (alta, baja, frontal) | 0° (pies rectos) | 12° |
| Peso muerto convencional | 15° | 15° |
| Peso muerto sumo | 45° | 25° |

Con pies rectos y los tobillos al ancho de hombros la rodilla, forzada hacia adelante, queda ≈ 10 cm por dentro del pie en el fondo (parece valgo; las caderas están a ±9 cm). Con 11–15° hacia afuera queda sobre el pie en cinco cuerpos probados (ejemplo, alto, bajo, caderas anchas, mujer). En sumo, la rodilla que sigue exactamente el pie (45–50°) queda por detrás del tobillo, porque la cadera calibrada está ~28 cm detrás de la barra, y el momento extensor de rodilla pasa a negativo; con 25° queda unos 9 cm por dentro del pie al despegar (aproximación entre la apariencia frontal y el signo del momento).

**Barra vertical sobre el medio pie (decisión de entrenador).** Las fichas calibradas ponían la barra en −3,2 cm (convencional) y −7,5 cm (sumo) respecto del medio pie, es decir, sobre el tobillo o detrás de él (el tobillo está ~7,7 cm detrás del medio pie en el atleta de ejemplo). Eso dejaba la barra dentro de la tibia y de la rótula, y obligaba a la barra a desviarse al cruzar la rodilla (el «bulto» del recorrido). Ahora `bar_x_offset = 0` en ambos pesos muertos y la barra sube en vertical sobre el medio pie durante todo el recorrido.
Mover solo la barra no basta: el hombro y la cadera cuelgan de ella, así que la rodilla avanza con ellos y el bulto seguía en 3,1 cm (convencional) y 5,8 cm (sumo). Se **reajustaron los otros cuatro parámetros de estilo** (`f_h0`, `eps0`, `eps_top`, `p`) contra los mismos objetivos publicados con la barra fija en 0, la rodilla 3D y la salvaguarda de abajo, con rangos por remuestreo (30 reajustes). El resultado es un desplazamiento de barra de 0,0 cm en todo el recorrido y una distancia mínima barra-rodilla de 6,4 cm (justo la holgura exigida).

**Salvaguarda de holgura (sigue activa).** Si el centro de la barra queda a menos de `KNEE_CLEAR` = 6,4 cm del centro de la rodilla (rótula 5 cm + radio de barra 1,4 cm, aproximación), la barra se adelanta lo justo respecto del cuerpo: hombro, cadera y rodilla no se mueven y los brazos se inclinan algo hacia adelante. Con los parámetros actuales no actúa para el atleta proxy ni el de ejemplo, pero protege medidas o ajustes del «Qué pasaría si» que la activen. Con la rodilla plana original, en el atleta de ejemplo el centro de la barra pasaba a 2,0 cm de la rodilla en sumo y a 4,7 cm en convencional. Se probó que la distancia mínima es ≥ 6,4 cm en tres atletas (ejemplo, bajo, alto).

**Bloqueo con extensión de espalda (peso muerto).** Dos parámetros de ficha, `bar_top_from_ankle` y `layback_deg` (por defecto `None` / 0: el bloqueo de antes, con el tronco vertical y todo sobre `bar_x_offset`). Con ellos: la barra termina a `bar_top_from_ankle` por delante del tobillo (10 cm en ambos: es la profundidad del muslo, así que la cadera queda sobre el tobillo y las piernas verticales; con la barra más cerca del tobillo y apoyada en el muslo las piernas quedaban inclinadas hacia atrás y la figura se veía cayendo hacia atrás); la cadera queda una profundidad de muslo (`BAR_TO_HIP` = 10,4 cm, aproximación) detrás de la barra; y el tronco termina `layback_deg` (10°) más allá de la vertical, con el hombro detrás de la cadera. La barra baja desde la vertical sobre el medio pie hacia ese x como `s^BAR_RETURN_P` (la pierna la adelanta si hace falta, como antes). El tronco lo dicta la rodilla (ver arriba) hasta `LOCKOUT_BLEND_S` = 0,85 y desde ahí se le suma de forma suavizada la diferencia hasta el tronco pedido, así que no hay saltos. El valor de `eps_top` no se usa en este modo (el hombro queda donde lo dejan la barra y el tronco). Son valores ilustrativos, sin fuente, y los parámetros de estilo ajustados (`eps0`, etc.) no se reajustaron con este bloqueo. En la figura el tronco se dibuja recto, inclinado hacia atrás, y se rotula la extensión.

**Efecto en la calibración y en los momentos.**
- Convencional: los nueve objetivos publicados siguen a ≤ 2 DE (máximo z +1,7, brazo de rodilla al despegar).
- Sumo: los brazos de momento publicados de rodilla (−13,2 y −4,8 cm) y de tobillo (−4,3 y −5,2 cm) ponen la barra **detrás** del centro de la articulación, lo que no es compatible con una barra sobre el medio pie que pasa por delante de la rótula y la tibia. Cuatro objetivos quedan entre 2,1 y 2,4 DE (rodilla al despegar +2,2, rodilla al pasar +2,1, tobillo al despegar +2,4, tobillo al pasar +2,4); los otros cinco, dentro de 2 DE. La prueba usa un umbral de 3 DE solo para esos cuatro.
- **Momento de rodilla.** Con la barra sobre el medio pie y la rodilla ≥ 6,4 cm detrás de ella, el momento extensor de rodilla del modelo es bajo o negativo en la primera fase (pico por pierna, atleta proxy: convencional −22 N·m, sumo −12 N·m; antes de estos cambios +2 y +32). La literatura da +53 y +80 N·m por pierna. Se conserva el orden sumo > convencional, pero **la magnitud del momento de rodilla no es fiable** en este modelo (ver limitación 13).
- Validación independiente actualizada (atleta proxy, 130 kg): tronco al despegar, convencional 3,6° (lit. 11 ± 7), sumo 14,5° (lit. 16 ± 6); momento de cadera pico conv/sumo 1,09 (lit. 1,10); diferencia por pierna en cadera conv−sumo 20 N·m (lit. +27); rodilla sumo−conv por pierna +10 N·m (lit. +27); recorrido de barra 54,8 vs 50,9 cm.

**Momentos.** Siguen siendo los del plano sagital (palancas horizontales × cos ψ). No se aplicó la inclinación del eje de flexión por la dirección del pie: las componentes laterales dependen de cómo se reparte la carga entre piernas, que el modelo no determina (probado: añadirlas daba un momento de tobillo de −475 N·m en sumo, no interpretable). El modelo sigue siendo simétrico y cuasi-estático.

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
| Peso muerto convencional | f_h0 0,842 [0,794–0,853]; barra x **0 (fija, medio pie)**; eps_top 0,018 [0,012–0,020]; p 1,56 [1,13–1,57]; eps0 0,12 **en el límite** (reajustado con la rodilla 3D; el ajuste anterior con barra libre dio f_h0 0,814 y barra −3,2 cm) |
| Peso muerto sumo | f_h0 0,776 [0,776–0,827]; barra x **0 (fija, medio pie)**; eps0 0,102 [0,102–0,107]; eps_top −0,015 [−0,015, 0,010]; p 1,21 [1,21–1,44] (reajustado con la rodilla 3D; el ajuste anterior con barra libre dio f_h0 0,762 y barra −7,5 cm) |
| Sentadilla barra alta | tibia máx. 32,8° [28,8–36,5]; profundidad −1,9° [−7,6, 6,6]; barra x +4,0 cm [1,5, 6,0] |
| Sentadilla barra baja | tibia máx. 20,8° [16,9–25,6]; profundidad +2,7° [−15, 8,7]; barra x −2,1 cm [−4,4, +1,1] |
| Desfase de cadera hx | **+5,8 cm [3,2, 8,4]** (inferido) |

## Validación independiente (no ajustada) — valores del ajuste original, anteriores a la rodilla 3D y a la barra sobre el medio pie
(Los valores vigentes están en «Dirección de pies y rodillas».)
- Ángulo del tronco sobre la horizontal al despegue: convencional 6,0° (lit. 11 ± 7), sumo 15,2° (lit. 16 ± 6).
- Momento de cadera pico conv/sumo: 1,15 (lit. 1,10). Diferencia por pierna en cadera conv−sumo: 31 N·m (lit. +27).
- Rodilla sumo−conv por pierna: 35 N·m (lit. +27). El recorrido de barra del sumo es más corto que el convencional (como en Escamilla: 47 vs 51 cm).

## Limitaciones conocidas (leer antes de usar)
1. **`hx` = +5,8 cm es una inferencia**, no una medida: el ajuste lo necesita para igualar los brazos de momento de cadera publicados y está confundido con diferencias de definición del centro de cadera. Anatómicamente plausible pero mayor de lo esperado. Pendiente de verificar con imagen (piloto).
2. **`eps0`, `eps_top`, `p` quedan en el borde del rango de búsqueda** en ambos pesos muertos (hombros ~12 cm por delante de la barra al inicio; cadera que se retrasa). Los datos vienen de levantadores de menor nivel con tronco casi horizontal; no extrapolar a levantadores avanzados.
3. Sentadilla: la inclinación del tronco queda **sobreestimada** (barra alta 45° vs 33,5 ± 4,6; barra baja 54° vs 33,1 ± 4,5). La posición de la barra en el tronco (`lam`, `d`) no es identificable con estos datos: se fijó a valores ilustrativos. Prueba con tolerancia explícita.
4. Brazos de momento de rodilla y tobillo en sumo: incompatibles con la barra por delante de la rótula y la tibia (ver «Dirección de pies y rodillas»); tres de ellos entre 2 y 2,5 DE. Barra baja, rodilla: z −1,9.
5. Definiciones de ángulo (tronco respecto a la horizontal, muslo, ancho de postura ≈ distancia entre tobillos, instante «inicio de fase concéntrica») fueron **inferidas** de los textos; verificar con los artículos originales.
6. Rotación externa del pie (`toe_out_deg`), dirección de la rodilla (`knee_out_deg`) y sentadilla frontal: sin fuente. El sumo a 45° es criterio de entrenador, no el ángulo publicado.
7. No usados por tablas ilegibles o definiciones ambiguas: Glassbrook y cols. (2017) y otros estudios de sentadilla. Fuentes de PubMed/PMC bloqueadas por reCAPTCHA (no se eludió): ISAK, ANSUR II, Wretenberg, Escamilla 2001 (sentadilla) no verificadas.
8. Cuasi-estático: subestima picos; los momentos de la literatura incluyen dinámica y son por pierna, solo se compararon razones/diferencias. Momentos de tobillo no comparados.
9. Lumbar: punto a 10 % de Ltr sobre la cadera (aproximación). Cabeza: ~6,9 % de H sobre el hombro (aproximación propia).
10. Los números de de Leva se transcribieron de un PDF: verificar contra la tabla original.
11. Sumo: momentos solo en el plano sagital; sin demanda frontal (abductores, rotadores). La posición lateral de la rodilla se calcula, pero no genera momentos frontales.
12. **Barra contra el muslo al bloqueo.** La barra ya no atraviesa la rodilla ni la tibia (sube sobre el medio pie), pero al bloqueo el modelo la deja ~1 cm detrás del eje del muslo, cuando en la realidad apoya en su cara anterior (unos 8–9 cm por delante del eje). Corregirlo cambia los momentos de cadera cerca del bloqueo y exige decidir cuánto se prioriza la literatura.
13. **Momento de rodilla en peso muerto.** La trayectoria de la rodilla es casi plana (la tibia casi no cambia de inclinación durante el tirón: ~3° al despegar, ~7° al pasar la rodilla), cuando en un levantamiento real la rodilla empieza por delante de la barra y retrocede al subir ella. Con la barra sobre el medio pie y la holgura de rótula, eso deja el momento extensor de rodilla bajo o negativo. Hace falta parametrizar el recorrido por la inclinación de la tibia (como en la sentadilla) en lugar de por la altura de la cadera; está pendiente.
