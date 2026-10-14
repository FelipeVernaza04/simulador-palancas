"""Demo: bandas de incertidumbre y comparación pareada sumo vs convencional (mismas perturbaciones de medida)."""
import numpy as np
from anthropometry import SAMPLE_MALE, from_raw, closure_checks
from techniques import TECHNIQUES
from model import solve
from uncertainty import draw_perturbations, monte_carlo, band, JOINTS

LOAD = 150.0
S_GRID = [0.0, 0.25, 0.5, 0.75]       # s=1 (bloqueo) se omite: el momento ahí lo sostienen estructuras pasivas
N = 1000

print("Cierres:", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in closure_checks(SAMPLE_MALE).items()})
draws = draw_perturbations(N, seed=1)
res = {}
for key in ("peso_muerto_convencional", "peso_muerto_sumo"):
    res[key] = monte_carlo(SAMPLE_MALE, TECHNIQUES[key], LOAD, S_GRID, draws)
    lo, mid, hi = band(res[key])
    print(f"\n{TECHNIQUES[key].name} ({LOAD:.0f} kg) — momento extensor, N·m  [mediana (P2.5–P97.5)]")
    print("   s   " + "".join(f"{j:>24}" for j in JOINTS))
    for i, s in enumerate(S_GRID):
        print(f"{s:5.2f}  " + "".join(f"{mid[i, k]:8.0f} ({lo[i, k]:4.0f}–{hi[i, k]:4.0f})  " for k in range(4)))
    print("   inalcanzables:", int(np.isnan(res[key][:, :, 0]).sum()), "de", res[key].shape[0] * len(S_GRID))

print("\nComparación pareada (sumo − convencional), misma persona y mismas perturbaciones:")
diff = res["peso_muerto_sumo"] - res["peso_muerto_convencional"]
for i, s in enumerate(S_GRID):
    row = []
    for k, j in enumerate(JOINTS):
        d = diff[:, i, k]
        d = d[~np.isnan(d)]
        row.append(f"{j}: mediana {np.median(d):+5.0f}, P(sumo<conv)={np.mean(d < 0):.2f}")
    print(f"s={s:.2f}  " + " | ".join(row))
