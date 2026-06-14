#!/usr/bin/env python3
"""outer_end 行から全 rms_* 残差を片対数でプロットする汎用ツール。"""
import sys, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else ""
rows = [r for r in csv.DictReader(open(f"{run}/residual_history.csv")) if r["phase"] == "outer_end"]
steps = np.array([int(r["step"]) for r in rows])
cols = [c for c in rows[0] if c.startswith("rms_") and not c.startswith("rms_dq")]
plt.figure(figsize=(7,5))
for c in cols:
    v = np.array([float(r[c]) for r in rows])
    if np.nanmax(v) <= 0: continue
    plt.semilogy(steps, np.clip(v, 1e-20, None), label=c)
plt.xlabel("pseudo-time step"); plt.ylabel("rms residual")
plt.title(f"Residual history {label}"); plt.grid(alpha=0.3, which="both"); plt.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"{run}/residual_history.png", dpi=130)
print("wrote", f"{run}/residual_history.png")
