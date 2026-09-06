#!/usr/bin/env python3
"""CMC の PDF 積分診断: cmc_Q.bin (Q(η) 全 node) と res_*.h5 (ξ̃, ξ''²) から β-PDF 重み Ω_k を cmc_pdf_weights と同じ規約で作り、
T_pdf = Σ Ω_k T(η_k) (条件付き T の PDF 平均) と輸送平均 T̃ を実験 (半径 z/d 9/11/14/26, 軸) と比較する。
文献 RANS-CMC の標準 (平均スカラーは Q から診断) で、平均場への熱の受け渡し (couple) の遅れを分離して CMC 自体の妥当性を見る。
  python3 cmc_pdf_mean.py run_dir cmc_Q_file res_step [--neta 41 --etapow 1.5]"""
import sys, argparse, pathlib, h5py, numpy as np, yaml, re
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("qfile"); ap.add_argument("step", type=int)
ap.add_argument("--neta", type=int, default=41); ap.add_argument("--etapow", type=float, default=1.5); ap.add_argument("--href", type=float, default=298.15)
a = ap.parse_args(); run = pathlib.Path(a.run); HERE = pathlib.Path(__file__).parent; D = 0.00457
cfg = open(run / "solverConfig.yaml").read(); names = re.search(r"species:\s*\[([^\]]*)\]", cfg).group(1).replace(" ", "").split(",")
db = yaml.safe_load(open(run / "species_db.yaml")); ns = len(names); RU = 8.314462618
MW = np.array([db[s]["MW"] for s in names]); lo = np.array([db[s]["nasa9_low"] for s in names]); hi = np.array([db[s]["nasa9_high"] for s in names])
Tmid = np.array([db[s]["Tmid"] for s in names]); Tlo = np.array([db[s]["Tlo"] for s in names]); Thi = np.array([db[s]["Thi"] for s in names])
def cph(T):   # T: (M,) → cp, h [J/mol] (M, ns) NASA-9 (範囲外は端でクランプ + h 線形外挿)
    Tc = np.clip(T[:, None], Tlo[None], Thi[None]); a9 = np.where((Tc < Tmid[None])[..., None], lo[None], hi[None])
    Ti = 1 / Tc; L = np.log(Tc)
    cp = RU * (a9[..., 0]*Ti**2 + a9[..., 1]*Ti + a9[..., 2] + a9[..., 3]*Tc + a9[..., 4]*Tc**2 + a9[..., 5]*Tc**3 + a9[..., 6]*Tc**4)
    h = RU * Tc * (-a9[..., 0]*Ti**2 + a9[..., 1]*L*Ti + a9[..., 2] + a9[..., 3]*Tc/2 + a9[..., 4]*Tc**2/3 + a9[..., 5]*Tc**3/4 + a9[..., 6]*Tc**4/5 + a9[..., 7]*Ti)
    return cp, h + cp * (T[:, None] - Tc)
_, href = cph(np.array([a.href])); href = href[0]
def h_mass(T, Y): _, h = cph(T); return ((h - href[None]) / MW[None] * Y).sum(1)
def T_from_h(h, Y, T0):
    T = T0.copy()
    for _ in range(30):
        cp, hh = cph(T); hm = ((hh - href[None]) / MW[None] * Y).sum(1); cpm = (cp / MW[None] * Y).sum(1)
        dT = np.clip((hm - h) / np.maximum(cpm, 1e-3), -0.5 * T, 0.5 * T); T = np.clip(T - dT, 200, 4000)
        if np.abs(dT).max() < 1e-3: break
    return T
# η 格子と重み (cmc_d.cu と同じ)
ne = a.neta; eta = (np.arange(ne) / (ne - 1)) ** a.etapow; eta[0] = 0; eta[-1] = 1
w = np.empty(ne); w[0] = 0.5 * (eta[0] + eta[1]); w[-1] = 1 - 0.5 * (eta[-2] + eta[-1]); w[1:-1] = 0.5 * (eta[2:] - eta[:-2])
def weights(xi, var):
    xi = np.clip(xi, 0, 1); vmax = xi * (1 - xi); Om = np.zeros((len(xi), ne))
    delta = (var < 1e-4 * vmax + 1e-10) | (vmax < 1e-12)
    k = np.clip(np.searchsorted(eta, xi, side="right") - 1, 0, ne - 2); f = np.clip((xi - eta[k]) / (eta[k + 1] - eta[k]), 0, 1)
    Om[np.arange(len(xi)), k] = 1 - f; Om[np.arange(len(xi)), k + 1] += f
    b_ = ~delta
    if b_.any():
        x, v, vm = xi[b_], var[b_], vmax[b_]; gam = np.maximum(vm / np.minimum(v, 0.999 * vm) - 1, 0.05); A = x * gam; B = (1 - x) * gam
        ln = np.full((b_.sum(), ne), -1e300)
        ln[:, 1:-1] = np.log(w[1:-1])[None] + (A[:, None] - 1) * np.log(eta[1:-1])[None] + (B[:, None] - 1) * np.log(1 - eta[1:-1])[None]
        ln[:, 0] = A * np.log(0.5 * eta[1]) - np.log(A); ln[:, -1] = B * np.log(0.5 * (1 - eta[-2])) - np.log(B)
        O = np.exp(ln - ln.max(1, keepdims=True)); Om[b_] = O / O.sum(1, keepdims=True)
    return Om
with h5py.File(run / "cabra.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
with h5py.File(run / f"res_{a.step}.h5") as h: V = h["VALUE"]; T = V["T"][:]; xi = V["xi"][:]; var = V["xiVar"][:]; Ym = np.stack([V[f"Y{s}"][:] for s in range(ns)], 1)
N = len(T); x = xyz[:N, 0]; y = xyz[:N, 1]
raw = np.fromfile(a.qfile, dtype=np.int32, count=3); nE, nS, nAll = raw; assert nE == ne and nS == ns, (nE, nS)
Q = np.fromfile(a.qfile, dtype=np.float32, offset=12).reshape(ns + 1, ne, nAll)[:, :, :N].astype(np.float64)   # [var][k][node]
Om = weights(xi.astype(np.float64), var.astype(np.float64))
# 条件付き T(η) を全 (node, η) で
Yq = np.clip(Q[:ns], 0, None); Yq /= np.maximum(Yq.sum(0, keepdims=True), 1e-30); hq = Q[ns]
Tq = np.empty((ne, N))
for k in range(ne): Tq[k] = T_from_h(hq[k], Yq[:, k].T, np.full(N, 1000.0))
Tpdf = (Om * Tq.T).sum(1); Ypdf = np.einsum("nk,skn->ns", Om, Yq); hpdf = (Om * hq.T).sum(1); T_hpdf = T_from_h(hpdf, Ypdf, Tpdf.copy())
print(f"step {a.step}: mean-field Tmax {T.max():.0f} K, T_pdf max {Tpdf.max():.0f} K, T(h_pdf,Y_pdf) max {T_hpdf.max():.0f} K, cond Tq max {Tq.max():.0f} K")
ax = np.where(np.abs(y) < 1e-9)[0]; o = np.argsort(x[ax]); ip = lambda f, zd: np.interp(zd * D, x[ax][o], f[ax][o])
print("axis z/d :  " + "  ".join(f"{zd:5d}" for zd in (9, 11, 14, 20, 26, 30)))
exc = np.loadtxt(HERE / "exp_centerline.csv", delimiter=",", comments="#"); ex = lambda zd: np.interp(zd, exc[:, 0], exc[:, 1])
for lab, f in (("mean T", T), ("T_pdf", Tpdf), ("T(h_pdf)", T_hpdf)): print(f"{lab:9s}:  " + "  ".join(f"{ip(f, zd):5.0f}" for zd in (9, 11, 14, 20, 26, 30)))
print("exp      :  " + "  ".join(f"{ex(zd):5.0f}" for zd in (9, 11, 14, 20, 26, 30)))
exr = np.loadtxt(HERE / "exp_radial_T.csv", delimiter=",", comments="#")
for zd, col in ((9, 2), (11, 3), (14, 4), (26, 5)):
    m = np.where(np.abs(x - zd * D) < 0.5 * D)[0]; r = y[m] * 1e3; oo = np.argsort(r); msk = exr[:, 0] <= 30
    out = []
    for lab, f in (("mean", T), ("pdf", Tpdf), ("h_pdf", T_hpdf)):
        Ti = np.interp(exr[:, 0], r[oo], f[m][oo]); d = Ti[msk] - exr[msk, col]; out.append(f"{lab} mean diff {d.mean():+5.0f} K (max |{np.abs(d).max():4.0f}|, r=0 {Ti[0]:.0f})")
    print(f"radial z/d={zd:2d} (exp r=0 {exr[0, col]:.0f}): " + " | ".join(out))
# 図
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 5, figsize=(17, 3.6))
zz = np.linspace(1, 40, 200); axs[0].plot(zz, [ip(T, z) for z in zz], label="mean T̃"); axs[0].plot(zz, [ip(Tpdf, z) for z in zz], label="Σ Ω T(η)"); axs[0].plot(zz, [ip(T_hpdf, z) for z in zz], "--", label="T(h_pdf,Y_pdf)")
axs[0].plot(exc[:, 0], exc[:, 1], "ko", ms=3, label="exp"); axs[0].set_title("centerline"); axs[0].set_xlabel("z/d"); axs[0].legend(fontsize=7)
for ai, (zd, col) in zip(axs[1:], ((9, 2), (11, 3), (14, 4), (26, 5))):
    m = np.where(np.abs(x - zd * D) < 0.5 * D)[0]; r = y[m] * 1e3; oo = np.argsort(r)
    ai.plot(r[oo], T[m][oo], label="mean"); ai.plot(r[oo], Tpdf[m][oo], label="pdf"); ai.plot(r[oo], T_hpdf[m][oo], "--", label="h_pdf"); ai.plot(exr[:, 0], exr[:, col], "ko", ms=3); ai.set_xlim(0, 30); ai.set_title(f"z/d={zd}"); ai.set_xlabel("r [mm]")
axs[1].legend(fontsize=7); fig.tight_layout(); fig.savefig(run / f"cmc_pdf_mean_{a.step}.png", dpi=110); print("fig:", run / f"cmc_pdf_mean_{a.step}.png")
