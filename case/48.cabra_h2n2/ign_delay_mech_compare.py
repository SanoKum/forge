#!/usr/bin/env python3
"""Cabra 混合線上の着火遅れ機構比較 (P0: Jachimowski の 1000-1100 K 妥当性検証)。

ジェット (25%H2/75%N2, 305 K) と coflow (T_c, X_O2 0.1474 / X_H2O 0.0989 / N2 bal.) を
混合分率 ξ (質量基準) で断熱混合し、定圧 (1 atm) 断熱リアクタの着火遅れ τ_ign (max dT/dt) を
ξ 掃引で計算する。機構: case/48 mech.yaml (Jachimowski 9sp20r) vs Burke 2012 / Li 2004 /
GRI3.0-H2 subset (既知の不良例)。coflow 温度感度 (1015/1045/1075 K) も出す。

  .venv-chem/bin/python ign_delay_mech_compare.py    # → ign_delay_mech_compare.{csv,png}
"""
import pathlib
import numpy as np
import cantera as ct
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).parent
MECH_DIR = HERE / "../../solver_density_cuda/tools/mechanisms"
MECHS = {
    "Jachimowski 9sp20r (forge)": str(HERE / "mech.yaml"),
    "Burke 2012": str(MECH_DIR / "h2_burke2012_cantera.yaml"),
    "Li 2004": str(MECH_DIR / "h2co_li2004_cantera.yaml"),
    "GRI3.0 H2 subset": "h2o2.yaml",
}
P = ct.one_atm
X_JET = "H2:0.25, N2:0.75"
T_JET = 305.0
X_COF = "O2:0.1474, H2O:0.0989, N2:0.7537"   # balance N2 (2026-09-05 codex 指摘: 旧値 0.8537 は合計 1.1 で Cantera が正規化していた)
TAU_MAX = 2.0  # [s] これ以上は「実質不着火」

def mix_state(gas, xi, t_cof):
    jet = ct.Quantity(gas, constant="HP"); jet.TPX = T_JET, P, X_JET; jet.mass = xi
    cof = ct.Quantity(gas, constant="HP"); cof.TPX = t_cof, P, X_COF; cof.mass = 1.0 - xi
    m = jet + cof   # 断熱等圧混合
    return m.T, m.Y.copy()

def tau_ign(mech, xi, t_cof):
    gas = ct.Solution(mech)
    T0, Y0 = mix_state(gas, xi, t_cof)
    gas.TPY = T0, P, Y0
    r = ct.IdealGasConstPressureReactor(gas)
    net = ct.ReactorNet([r])
    t, tau, dTdt_max = 0.0, None, 0.0
    t_prev, T_prev = 0.0, r.T
    while t < TAU_MAX:
        t = net.step()
        dTdt = (r.T - T_prev) / max(t - t_prev, 1e-30)
        if dTdt > dTdt_max:
            dTdt_max, tau = dTdt, t
        t_prev, T_prev = t, r.T
        if r.T > T_prev + 0.0:  # noop guard
            pass
        if r.T - T0 > 150.0 and t > 3 * (tau or t):  # 着火後は十分回ったら打ち切り
            break
    if r.T - T0 < 50.0:
        return np.nan, T0  # 不着火扱い
    return tau, T0

XIS = np.round(np.geomspace(0.004, 0.5, 25), 5)
T_COFS = [1015.0, 1045.0, 1075.0]

rows = []
for name, mech in MECHS.items():
    for t_cof in T_COFS:
        for xi in XIS:
            tau, T0 = tau_ign(mech, xi, t_cof)
            rows.append((name, t_cof, xi, T0, tau))
            print(f"{name:28s} Tc={t_cof:.0f} xi={xi:.4f} T0={T0:6.1f} tau={tau if np.isnan(tau) else f'{tau*1e3:.3f} ms'}")

import csv
with open(HERE / "ign_delay_mech_compare.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["mech", "T_coflow", "xi", "T_mix", "tau_ign_s"]); w.writerows(rows)

# --- 図: (a) tau(xi) @1045K 全機構, (b) tau_min(T_c) 感度 ---
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
colors = {"Jachimowski 9sp20r (forge)": "tab:red", "Burke 2012": "tab:blue",
          "Li 2004": "tab:green", "GRI3.0 H2 subset": "tab:gray"}
for name in MECHS:
    d = [(xi, tau) for (n, tc, xi, T0, tau) in rows if n == name and tc == 1045.0]
    xs, ts = zip(*d)
    ax[0].loglog(xs, np.array(ts) * 1e3, "o-", ms=3.5, color=colors[name], label=name)
ax[0].set_xlabel("mixture fraction xi"); ax[0].set_ylabel("tau_ign [ms]")
ax[0].set_title("Cabra mixing line, T_c = 1045 K, 1 atm"); ax[0].grid(alpha=.3, which="both"); ax[0].legend(fontsize=8)
ax[0].axhline(1.0, color="k", lw=0.8, ls=":"); ax[0].text(0.3, 1.1, "~conv. time to H/d=10", fontsize=7)
for name in MECHS:
    tmins, ximins = [], []
    for tc in T_COFS:
        d = [(xi, tau) for (n, t, xi, T0, tau) in rows if n == name and t == tc and np.isfinite(tau)]
        if d:
            xi_m, tau_m = min(d, key=lambda p: p[1])
            tmins.append(tau_m * 1e3); ximins.append(xi_m)
        else:
            tmins.append(np.nan); ximins.append(np.nan)
    ax[1].semilogy(T_COFS, tmins, "o-", color=colors[name], label=name)
    for tc, tm, xm in zip(T_COFS, tmins, ximins):
        if np.isfinite(tm): ax[1].annotate(f"xi={xm:.3f}", (tc, tm), fontsize=6, xytext=(3, 3), textcoords="offset points")
ax[1].set_xlabel("T_coflow [K]"); ax[1].set_ylabel("min tau_ign over xi [ms]")
ax[1].set_title("most-reactive ignition delay vs coflow T"); ax[1].grid(alpha=.3, which="both"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(HERE / "ign_delay_mech_compare.png", dpi=130)
print("wrote ign_delay_mech_compare.{csv,png}")
