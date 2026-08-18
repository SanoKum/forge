"""NASA CEA2 (凍結流, 5 成分固定) と forge node Euler TP 軸解の比較 — va3 新条件 (Pt 11.39 bar, Tt 1161 K, 新組成; M4.19 va3 L_c8 dry run_0091).
実行: python3 case/44.vitiated_air_wt/cea_check_va3.py [run_dir]
前提: cea/va3_cea.out (NASA CEA2 の出力; cea/ で `echo va3_cea | ../../../.venv-cea/nasa_cea/FCEA2`; 入力 cea/va3_cea.inp)."""
import re, sys, json
from pathlib import Path
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE = Path(__file__).resolve().parent
RUN = CASE / (sys.argv[1] if len(sys.argv) > 1 else "run_0091_va3_M4.19_Lc8_dry")
PT, TT = 1.139e6, 1161.0     # va3 新条件
S = 0.205711  # r_t [m] (M4.19 va3)
CEA_OUT, CEA_CASE, SUPAR = "va3_cea.out", "v3_fz", 15.066   # supar: CFD 実出口 A/A* (設計 15.124 も va3_cea.inp に入れてある)

# ---------- CEA 出力の解析 ----------
def cea_num(tok):
    # "3.7615 0" / "8.7543-1" 形式 (RHO 行) と通常の数値の両方
    m = re.fullmatch(r"([0-9.]+)([+-]\d+)", tok)
    return float(m.group(1)) * 10.0 ** int(m.group(2)) if m else float(tok)

def parse_cea(path, case_prefix):
    """case_prefix で始まる CASE の全テーブルから列を集める (CHAMBER/THROAT 列は重複するので 1 度だけ)."""
    txt = path.read_text().splitlines()
    rows = {"Pinf/P": [], "P": [], "T": [], "RHO": [], "H": [], "Cp": [], "GAMMAs": [], "a": [], "M": [], "AeAt": [], "u": []}
    i = 0; first = True
    while i < len(txt):
        L = txt[i]
        if L.strip().startswith("CASE =") and L.split("=")[1].strip().startswith(case_prefix):
            # テーブル本体まで進む
            j = i
            while not txt[j].lstrip().startswith("Pinf/P"): j += 1
            blk = {}
            for k in range(j, j + 40):
                s = txt[k]
                key = None
                for kk, pat in [("Pinf/P", "Pinf/P"), ("P", "P, BAR"), ("T", "T, K"), ("RHO", "RHO, KG/CU M"), ("H", "H, KJ/KG"),
                                ("Cp", "Cp, KJ/(KG)(K)"), ("GAMMAs", "GAMMAs"), ("a", "SON VEL,M/SEC"), ("M", "MACH NUMBER"),
                                ("AeAt", "Ae/At"), ("u", "Isp, M/SEC")]:
                    if s.lstrip().startswith(pat): key = kk; body = s.split(pat, 1)[1]; break
                if key is None: continue
                toks = body.split()
                if key == "RHO":  # "3.7615 0" は 2 トークンに割れる → 結合
                    vals = []; t = 0
                    while t < len(toks):
                        if re.fullmatch(r"[0-9.]+", toks[t]) and t + 1 < len(toks) and re.fullmatch(r"[+-]?\d+", toks[t + 1]) and "." not in toks[t + 1]:
                            vals.append(float(toks[t]) * 10.0 ** int(toks[t + 1])); t += 2
                        else:
                            vals.append(cea_num(toks[t])); t += 1
                else:
                    vals = [cea_num(t) for t in toks]
                blk[key] = vals
            n = len(blk["Pinf/P"])
            for key in rows:
                v = blk[key]
                if key in ("AeAt", "u"):  # CHAMBER 列が無い
                    v = [np.nan] + v
                v = v[:n]
                rows[key].extend(v if first else v[2:])  # 2 度目以降は CHAMBER/THROAT を落とす
            first = False
            i = j + 40
        else:
            i += 1
    d = {k: np.array(v, float) for k, v in rows.items()}
    o = np.argsort(d["Pinf/P"])
    return {k: v[o] for k, v in d.items()}

cea = parse_cea(CASE / "cea" / CEA_OUT, CEA_CASE)
# 出口 supar=14.443 の行 (Pinf/P=217.77 の重複) を分離
i_sup = np.where(np.isfinite(cea["AeAt"]) & (np.abs(cea["AeAt"] - SUPAR) < 1e-3))[0]
cea_exit = {k: v[i_sup[0]] for k, v in cea.items()}
print(f"CEA frozen: chamber T={cea['T'][0]:.2f} K rho={cea['RHO'][0]:.4f} cp={cea['Cp'][0]:.4f} gamma={cea['GAMMAs'][0]:.4f}")
it = np.argmin(np.abs(cea["M"] - 1.0))
print(f"CEA frozen: throat Pc/P*={cea['Pinf/P'][it]:.4f} T*={cea['T'][it]:.2f} rho*={cea['RHO'][it]:.4f} a*={cea['a'][it]:.1f} gamma*={cea['GAMMAs'][it]:.4f}")
print(f"CEA frozen exit (Ae/At={SUPAR}): M={cea_exit['M']:.3f} P={cea_exit['P']*1e5:.0f} Pa T={cea_exit['T']:.2f} K rho={cea_exit['RHO']:.5f} u={cea_exit['u']:.1f} a={cea_exit['a']:.1f}")

# ---------- forge 軸データ ----------
info = json.loads((RUN / "prepare_info.json").read_text())
ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]
with h5py.File(RUN / "nozzle.h5") as n:
    xy = n["MESH/COORD"][:].reshape(-1, 3)[:, :2]
X = xy[:, 0].reshape(ni, nj); Y = xy[:, 1].reshape(ni, nj)
assert np.all(Y[:, 0] < 1e-9), "j=0 が軸ではない"
res = sorted(RUN.glob("res_*.h5"), key=lambda p: int(p.stem.split("_")[-1]))
res = [p for p in res if p.stem.count("_") == 1][-1]
with h5py.File(res) as f:
    g = lambda k: f["VALUE/" + k][:].reshape(ni, nj)[:, 0].astype(float)
    P, T, ro, Ux, Uy, a, roe = g("P"), g("T"), g("ro"), g("Ux"), g("Uy"), g("sonic"), g("roe")
xa = X[:, 0]; u = np.hypot(Ux, Uy); M = u / a
h0 = roe / ro + P / ro  # 全エンタルピー (forge datum)
print(f"forge {res.name}: axis nodes {ni}, x∈[{xa[0]:.3f},{xa[-1]:.3f}] m, P∈[{P.min():.0f},{P.max():.0f}], T∈[{T.min():.1f},{T.max():.1f}], M_max={M.max():.4f}")
print(f"forge axis total enthalpy h0: mean {h0.mean():.1f} J/kg, spread (max-min)/|mean| = {(h0.max()-h0.min())/abs(h0.mean()):.2e}, rel std {h0.std()/abs(h0.mean()):.2e}")
# 有効 R (P/(rho T)) と CEA の R
Rf = P / (ro * T); R_cea = 8314.462618 / 29.146
print(f"forge R=P/(rho T): {Rf.mean():.3f} ± {Rf.std():.3f}  (CEA MW 29.146 → R={R_cea:.3f}; design R=285.270)")

# 単調区間 (軸 M が最大となる点まで) で Pinf/P をパラメータに補間
imax = int(np.argmax(M))
pr = PT / P[: imax + 1]; ok = np.r_[True, np.diff(pr) > 0]
# 単調でない点は落とす (稀な数値ノイズ対策)
idx = [0]
for i in range(1, imax + 1):
    if pr[i] > pr[idx[-1]]: idx.append(i)
idx = np.array(idx)
def fi(arr): return np.interp(np.log(cea["Pinf/P"]), np.log(pr[idx]), arr[idx], left=np.nan, right=np.nan)
cmp = {"Pinf/P": cea["Pinf/P"], "M_cea": cea["M"], "T_cea": cea["T"], "rho_cea": cea["RHO"], "u_cea": cea["u"], "a_cea": cea["a"],
       "T_forge": fi(T), "rho_forge": fi(ro), "u_forge": fi(u), "a_forge": fi(a), "M_forge": fi(M), "x_forge": fi(xa)}
cmp["u_cea"][0] = 0.0  # chamber
rows = []
print("\n  Pc/P     x[m]   M_cea  M_forge  dM%    T_cea  T_forge  dT%    rho_cea rho_forge drho%   u_cea  u_forge  du%    a_cea  a_forge  da%")
for k in range(len(cea["Pinf/P"])):
    if not np.isfinite(cmp["T_forge"][k]): continue
    r = {q: cmp[q][k] for q in cmp}
    r["dT"] = 100 * (r["T_forge"] / r["T_cea"] - 1); r["drho"] = 100 * (r["rho_forge"] / r["rho_cea"] - 1)
    r["du"] = 100 * (r["u_forge"] / r["u_cea"] - 1) if r["u_cea"] > 0 else np.nan
    r["da"] = 100 * (r["a_forge"] / r["a_cea"] - 1); r["dM"] = 100 * (r["M_forge"] - r["M_cea"]) / max(r["M_cea"], 1e-9) if r["M_cea"] > 0 else np.nan
    rows.append(r)
    print(f"{r['Pinf/P']:8.3f} {r['x_forge']:7.3f}  {r['M_cea']:5.3f}  {r['M_forge']:6.4f} {r['dM']:6.3f}  {r['T_cea']:7.2f} {r['T_forge']:7.2f} {r['dT']:6.3f}  "
          f"{r['rho_cea']:7.4f} {r['rho_forge']:8.4f} {r['drho']:6.3f}  {r['u_cea']:6.1f} {r['u_forge']:7.1f} {r['du']:6.3f}  {r['a_cea']:6.1f} {r['a_forge']:6.1f} {r['da']:6.3f}")
import csv
with open(CASE / "cea" / f"cea_vs_forge_{RUN.name}.csv", "w", newline="") as fcsv:
    w = csv.DictWriter(fcsv, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
sup = [r for r in rows if r["Pinf/P"] > 1.05 and r["M_cea"] > 1.0]
for q in ("dT", "drho", "du", "da", "dM"):
    v = np.array([r[q] for r in sup]); print(f"supersonic stations (M>1, n={len(v)}): {q} max|.| = {np.nanmax(np.abs(v)):.3f} %  mean = {np.nanmean(v):+.3f} %")

# ---------- 出口断面: 質量流量・Cd・面平均 ----------
i_E = ni - 1  # 出口断面 (壁出口 x_F); x_E は軸 M 法則の終端で出口ではない
with h5py.File(res) as f:
    G = lambda k: f["VALUE/" + k][:].reshape(ni, nj).astype(float)
    P2, T2, ro2, Ux2, Uy2, a2 = G("P"), G("T"), G("ro"), G("Ux"), G("Uy"), G("sonic")
def mdot_at(i):
    r = Y[i]; return np.trapezoid(ro2[i] * Ux2[i] * 2 * np.pi * r, r)
A_star = np.pi * S**2
mdot_1d = cea["RHO"][it] * cea["a"][it] * A_star
md_E, md_in, md_last = mdot_at(i_E), mdot_at(0), mdot_at(ni - 1)
print(f"\nmass flow: CEA 1D rho* a* A* = {mdot_1d:.2f} kg/s (A*={A_star:.5f} m2); forge inlet {md_in:.2f}, exit {md_E:.2f} kg/s")
print(f"  Cd = forge/CEA-1D: inlet {md_in/mdot_1d:.4f}, exit {md_E/mdot_1d:.4f}   (design Hall/series Cd = {info['cd_series']:.4f})")
r = Y[i_E]; w = ro2[i_E] * Ux2[i_E] * r
def mavg(q): return np.trapezoid(q * w, r) / np.trapezoid(w, r)
print(f"exit plane x={xa[i_E]:.3f} m: r_wall={r[-1]:.4f} m A/A*={ (r[-1]/S)**2:.3f}; forge mass-avg M={mavg(np.hypot(Ux2[i_E],Uy2[i_E])/a2[i_E]):.4f} T={mavg(T2[i_E]):.2f} P={mavg(P2[i_E]):.0f} u={mavg(Ux2[i_E]):.1f}; axis M={M[i_E]:.4f} T={T[i_E]:.2f} P={P[i_E]:.0f}")
print(f"  CEA frozen 1D @A/A*={SUPAR}: M={cea_exit['M']:.3f} T={cea_exit['T']:.2f} P={cea_exit['P']*1e5:.0f} u={cea_exit['u']:.1f}")

# ---------- 図 ----------
fig, axs = plt.subplots(2, 3, figsize=(16, 8.5))
mk = dict(ls="none", marker="o", mfc="none", ms=6, color="k", label="NASA CEA2 (frozen, 5 species)")
prf = PT / P
for ax, (q_f, q_c, lab) in zip(axs.flat, [(T, "T", "T [K]"), (ro, "RHO", "ρ [kg/m³]"), (u, "u", "u [m/s]"), (a, "a", "sound speed [m/s]"), (M, "M", "Mach")]):
    ax.plot(prf, q_f, "r-", lw=1.2, label=f"forge axis ({RUN.name})")
    ax.plot(cea["Pinf/P"], cea[q_c], **mk); ax.set_xscale("log"); ax.set_xlabel("Pt / P"); ax.set_ylabel(lab); ax.grid(alpha=0.3, which="both")
axs.flat[0].legend(fontsize=8)
ax = axs.flat[5]
for q, lab in [("dT", "T"), ("drho", "ρ"), ("du", "u"), ("da", "a"), ("dM", "M")]:
    ax.plot([r["M_cea"] for r in rows], [r[q] for r in rows], marker="o", ms=4, lw=1, label=lab)
ax.axhline(0, c="k", lw=0.5); ax.set_xlabel("M (CEA)"); ax.set_ylabel("(forge − CEA)/CEA [%]"); ax.set_ylim(-0.5, 0.5); ax.grid(alpha=0.3); ax.legend(fontsize=8, ncol=5)
fig.suptitle(f"case/44 vitiated air va3 Pt=11.39 bar Tt=1161 K: forge node Euler TP axis vs NASA CEA2 frozen expansion ({RUN.name}, {res.name})", fontsize=11)
fig.tight_layout(); out = CASE / "cea" / f"cea_vs_forge_{RUN.name}.png"; fig.savefig(out, dpi=120); print("saved", out)

# ---------- 追加出力: CEA 凍結流テーブル (全 Pc/P 点) と forge 軸分布 (全ノード) を別 CSV に ----------
with open(CASE / "cea" / f"cea_frozen_table_{CEA_CASE}.csv", "w", newline="") as fcsv:
    keys = ["Pinf/P","P","T","RHO","H","Cp","GAMMAs","a","M","AeAt","u"]
    w = csv.writer(fcsv); w.writerow(["# NASA CEA2 frozen (nfz=1) 5-species, Pt=%.0f Pa Tt=%.0f K; P in bar, u=Isp[m/s]" % (PT, TT)]); w.writerow(keys)
    for k in range(len(cea["Pinf/P"])): w.writerow([cea[q][k] for q in keys])
with open(CASE / "cea" / f"forge_axis_{RUN.name}.csv", "w", newline="") as fcsv:
    w = csv.writer(fcsv); w.writerow(["# forge %s axis (j=0) values; Pt/P for matching with CEA" % res.name])
    w.writerow(["x_m","x_over_rt","Pt_over_P","P","T","rho","u","a","M","h0"])
    for k in range(ni): w.writerow([xa[k], xa[k]/S, PT/P[k], P[k], T[k], ro[k], u[k], a[k], M[k], h0[k]])
print("wrote", CASE / "cea" / f"cea_frozen_table_{CEA_CASE}.csv", "and", CASE / "cea" / f"forge_axis_{RUN.name}.csv")
