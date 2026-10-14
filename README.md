# Simulador de palancas — sentadilla y peso muerto

Herramienta didáctica para entrenadores: a partir de medidas con cinta métrica de un atleta, simula las palancas (brazos de momento y momentos articulares) en sentadilla (barra alta, baja, frontal) y peso muerto (convencional, sumo), y las muestra con una figura que se mueve con el recorrido. Incluye un **modo explicar** para enseñar el movimiento en vivo y una pestaña «Qué pasaría si» para ver el efecto de ajustes (barra, cadera, postura, profundidad).

**Prototipo para pruebas piloto. No mide la técnica real del atleta ni sustituye la observación del entrenador.**

## Probar la app
Abrir `index.html` en el navegador (es un solo archivo, sin instalación; los datos se calculan en el navegador y no se envían a ningún sitio). Pensada para celular. Con GitHub Pages activado, la dirección es la del repositorio.

## Contenido
- `index.html`: la app, ya construida.
- `app/`: fuentes de la interfaz (`app.src.html`), el modelo en JavaScript (`model.js`) y `build.py` (une ambos y escribe `../index.html`: `cd app && python3 build.py`).
- `modelo/`: prototipo de cálculo en Python (referencia del modelo, calibración y pruebas: `cd modelo && python3 -m pytest -q`). Ver `modelo/README.md` para el método, la calibración y las limitaciones conocidas.

## Qué hay que saber antes de usarla
- Modelo sagital 2D, cuasi-estático y simétrico; calibrado con datos publicados de levantadores hombres (Escamilla y cols., 2000; Swinton y cols., 2012). No se extrapola a mujeres ni a levantadores avanzados.
- Postura estándar: tobillos bajo los hombros y pies rectos (sumo conserva su postura ancha).
- Algunas cifras son inferidas (desfase de cadera de +5,8 cm) o aproximadas (zona lumbar), y los textos de las indicaciones son borradores para validar con entrenadores. Detalle en `modelo/README.md`.
- Las medidas se toman con el protocolo de medición con cinta (12 medidas).
