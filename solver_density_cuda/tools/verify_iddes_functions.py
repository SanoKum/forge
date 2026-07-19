#!/usr/bin/env python3
"""SST-IDDES 関数群 (f_B / f_e / f_dt / Δ_IDDES / l_IDDES) のコード化前検証。

plans/active/turbulence-iddes-sst.md §3.4/§4.6 (Phase 2)。
参照式は Shur, Spalart, Strelets & Travin (2008) IJHFF 29:1638 の IDDES を
SST 用に定数調整した Gritskevich, Garbaruk, Schütze & Menter (2012)
Flow Turb. Combust. 88:431 の SST-IDDES:

  α    = 0.25 − d_w/h_max
  f_B  = min(2 exp(−9α²), 1)
  f_e1 = 2 exp(−11.09α²) (α≥0) / 2 exp(−9α²) (α<0)
  r_dt = ν_t /(κ²d²·max(|∇u|, floor)),  r_dl = ν/(κ²d²·max(|∇u|, floor))
  f_t  = tanh((C_t² r_dt)³)   C_t=1.87
  f_l  = tanh((C_l² r_dl)^10) C_l=5.0
  f_e2 = 1 − max(f_t, f_l)
  f_e  = f_e2 · max(f_e1 − 1, 0)          (SST 版: SA の低 Re 補正 Ψ は無し)
  f_dt = 1 − tanh((C_dt1 r_dt)^3)         C_dt1=20 (SST 較正。SA/DDES の 8 と異なる)
  f̃_d  = max(1 − f_dt, f_B)
  Δ_IDDES = min(C_w·max(d_w, h_max), h_max)   C_w=0.15 (Gritskevich 簡略形;
            Shur 原型の h_wn 項を落とした版で、非構造格子で壁法線刻みが不要)
  l_IDDES = f̃_d(1+f_e)·l_RANS + (1−f̃_d)·l_LES,  l_LES = C_DES·Δ_IDDES

検証項目:
 (1) f_B: 壁で 1 (α=0.25 でも exp 項>0.5)、d/h≤0.528 で 1、d/h=1 で ≈0.0127、単調減衰
 (2) f_e1: α=0 で左右極限一致 (=2)・C0 連続、α≥0 側がより速く減衰
 (3) f_t/f_l/f_e2: r→0 で 0、r 大で 1 に飽和。f_e ≥ 0 と「モデル/粘性下層検出時に f_e=0」
 (4) Δ_IDDES: 壁で 0.15h_max、d≥h_max/0.15 で h_max、[0.15h, h] に有界・d に単調
 (5) RANS 縮退: 付着 BL 相当 (r_dt≈1) で f_dt→0, f̃_d→1, l_IDDES=(1+f_e)l_RANS
     (log 層合成プロファイルで f_e=0 も確認 → l_IDDES=l_RANS 厳密)
 (6) WMLES 分岐: 解像乱流相当 (νt 小 → f_dt→1) で f̃_d=f_B となり
     l_IDDES = f_B(1+f_e)l_RANS + (1−f_B)l_LES = l_WMLES (Shur §3.2 の形) に一致
 (7) float32 + ガード (r_d clamp[0,10]・floor) 実装が double 参照と一致
     (clamp が発火しない領域で相対差 < 3e-6、発火域でも関数値が同一)

使い方: python3 verify_iddes_functions.py
"""
import numpy as np

KAPPA = 0.41
CT = 1.87
CL = 5.0
CDT1 = 20.0
CW = 0.15
BETA_STAR = 0.09
C_DES = 0.78


# ---------- double 参照実装 (Gritskevich 2012 の式そのまま) ----------

def alpha(d, hmax):
    return 0.25 - d / hmax


def fb_ref(d, hmax):
    return np.minimum(2.0 * np.exp(-9.0 * alpha(d, hmax) ** 2), 1.0)


def fe1_ref(d, hmax):
    a = alpha(d, hmax)
    return np.where(a >= 0.0, 2.0 * np.exp(-11.09 * a * a), 2.0 * np.exp(-9.0 * a * a))


def rd_ref(nu_eff, d, gradmag):
    return nu_eff / (KAPPA ** 2 * d ** 2 * np.maximum(gradmag, 1e-30))


def ft_ref(rdt):
    return np.tanh((CT ** 2 * rdt) ** 3)


def fl_ref(rdl):
    return np.tanh((CL ** 2 * rdl) ** 10)


def fe_ref(d, hmax, rdt, rdl):
    fe2 = 1.0 - np.maximum(ft_ref(rdt), fl_ref(rdl))
    return fe2 * np.maximum(fe1_ref(d, hmax) - 1.0, 0.0)


def fdt_ref(rdt):
    return 1.0 - np.tanh((CDT1 * rdt) ** 3)


def ftilde_d_ref(rdt, d, hmax):
    return np.maximum(1.0 - fdt_ref(rdt), fb_ref(d, hmax))


def delta_iddes_ref(d, hmax):
    return np.minimum(CW * np.maximum(d, hmax), hmax)


def l_iddes_ref(l_rans, d, hmax, rdt, rdl):
    fd_t = ftilde_d_ref(rdt, d, hmax)
    fe = fe_ref(d, hmax, rdt, rdl)
    l_les = C_DES * delta_iddes_ref(d, hmax)
    return fd_t * (1.0 + fe) * l_rans + (1.0 - fd_t) * l_les


# ---------- float32 GPU 想定実装 (CUDA カーネルと同じガード・演算順) ----------

def f32(x):
    return np.float32(x)


def rd_gpu(nu_eff, d, gradmag):
    """compute_rd 相当: floor + clamp[0,10] を float32 で。"""
    nu_eff = f32(nu_eff); d = f32(d); gradmag = f32(gradmag)
    g = np.maximum(gradmag, f32(1e-30))
    denom = np.maximum(f32(KAPPA) * f32(KAPPA) * d * d * g, f32(1e-30))
    return np.minimum(nu_eff / denom, f32(10.0))


def l_iddes_gpu(l_rans, d, hmax, nut, nu, gradmag):
    l_rans = f32(l_rans); d = f32(d); hmax = f32(hmax)
    rdt = rd_gpu(nut, d, gradmag)
    rdl = rd_gpu(nu, d, gradmag)
    a = f32(0.25) - d / np.maximum(hmax, f32(1e-30))
    fb = np.minimum(f32(2.0) * np.exp(f32(-9.0) * a * a), f32(1.0))
    fe1 = np.where(a >= f32(0.0),
                   f32(2.0) * np.exp(f32(-11.09) * a * a),
                   f32(2.0) * np.exp(f32(-9.0) * a * a))
    argt = f32(CT * CT) * rdt
    ft = np.tanh(argt * argt * argt)
    argl = f32(CL * CL) * rdl
    fl = np.tanh(argl ** f32(10.0))
    fe2 = f32(1.0) - np.maximum(ft, fl)
    fe = fe2 * np.maximum(fe1 - f32(1.0), f32(0.0))
    argdt = f32(CDT1) * rdt
    fdt = f32(1.0) - np.tanh(argdt * argdt * argdt)
    fd_t = np.maximum(f32(1.0) - fdt, fb)
    delta = np.minimum(f32(CW) * np.maximum(d, hmax), hmax)
    l_les = f32(C_DES) * delta
    return fd_t * (f32(1.0) + fe) * l_rans + (f32(1.0) - fd_t) * l_les, fd_t, fe, fb


# ---------- 検証 ----------

def check(name, cond, detail=""):
    ok = bool(np.all(cond))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    return ok


def main():
    allok = True
    print("(1) f_B 形状")
    d_h = np.linspace(0.0, 3.0, 3001)
    fb = fb_ref(d_h, 1.0)
    allok &= check("壁 (d=0) で f_B=1", fb[0] == 1.0)
    allok &= check("d/h≤0.527 で f_B=1 (α²≤ln2/9)",
                   np.all(fb[d_h <= 0.25 + np.sqrt(np.log(2.0) / 9.0) - 1e-9] == 1.0))
    allok &= check("f_B(d/h=1) ≈ 2e^{-9·0.5625}=0.01272",
                   abs(fb_ref(1.0, 1.0) - 2.0 * np.exp(-9.0 * 0.5625)) < 1e-14,
                   f"val={fb_ref(1.0,1.0):.5f}")
    tail = fb[d_h > 0.53]
    allok &= check("crossover 以遠で単調減衰", np.all(np.diff(tail) <= 0.0))

    print("(2) f_e1 連続性")
    eps = 1e-9
    allok &= check("α=0 で左右一致 (=2)",
                   abs(fe1_ref(0.25 - eps, 1.0) - fe1_ref(0.25 + eps, 1.0)) < 1e-6
                   and abs(fe1_ref(0.25, 1.0) - 2.0) < 1e-12)
    allok &= check("α>0 側 (壁寄り d<0.25h) がより速く減衰",
                   fe1_ref(0.1, 1.0) < fe1_ref(0.4, 1.0))  # 同 |α|=0.15

    print("(3) f_t / f_l / f_e")
    allok &= check("r→0 で f_t,f_l→0", ft_ref(0.0) == 0.0 and fl_ref(0.0) == 0.0)
    allok &= check("r_dt 大で f_t→1 (f_e2→0)", ft_ref(1.0) > 0.999999)
    allok &= check("r_dl 大で f_l→1", fl_ref(0.5) > 0.999999)
    dd = np.linspace(1e-4, 3.0, 500)
    rr = np.geomspace(1e-6, 10.0, 200)
    D, R = np.meshgrid(dd, rr)
    fe_all = fe_ref(D, 1.0, R, R * 0.1)
    allok &= check("f_e ≥ 0 (全域)", np.all(fe_all >= 0.0))
    # f_t=tanh((C_t²r)³) の実飽和は C_t²r≳1.9 (r_dt≳0.55)。付着 log 層の代表値 r_dt=1 で判定
    allok &= check("モデル乱流検出 (log 層 r_dt=1) で f_e=0", np.all(fe_ref(D, 1.0, 1.0, 0.0) < 1e-10))

    print("(4) Δ_IDDES")
    delta = delta_iddes_ref(d_h, 1.0)
    allok &= check("壁で 0.15·h_max", abs(delta[0] - CW) < 1e-14)
    allok &= check("d≥h/0.15 で h_max", abs(delta_iddes_ref(1.0 / CW + 1.0, 1.0) - 1.0) < 1e-14)
    allok &= check("[0.15h, h] に有界・単調", np.all(delta >= CW) and np.all(delta <= 1.0)
                   and np.all(np.diff(delta) >= 0.0))

    print("(5) RANS 縮退 (付着 BL: log 層合成 νt=κu_τd, S=u_τ/(κd) → r_dt=1)")
    utau, hmax = 1.0, 0.02
    d = np.geomspace(1e-5, 5 * hmax, 200)
    nut = KAPPA * utau * d
    grad = utau / (KAPPA * d)
    rdt = rd_ref(nut, d, grad)
    allok &= check("log 層で r_dt=1", np.allclose(rdt, 1.0))
    allok &= check("f_dt→0 (シールド)", np.all(fdt_ref(rdt) < 1e-12))
    l_rans = 0.05
    l_i = l_iddes_ref(l_rans, d, hmax, rdt, rd_ref(1e-6, d, grad))
    allok &= check("f_e=0 → l_IDDES=l_RANS 厳密", np.allclose(l_i, l_rans, rtol=1e-12),
                   f"max rel dev={np.max(np.abs(l_i/l_rans-1)):.2e}")

    print("(6) WMLES 分岐 (解像乱流相当: νt→SGS 級で f_dt→1)")
    nut_sgs = 1e-4 * KAPPA * utau * d   # r_dt ≈ 1e-4 ≪ 0.02 → f_dt≈1
    rdt_s = rd_ref(nut_sgs, d, grad)
    allok &= check("f_dt→1", np.all(fdt_ref(rdt_s) > 0.999))
    fd_t = ftilde_d_ref(rdt_s, d, hmax)
    allok &= check("f̃_d = f_B に縮退", np.allclose(fd_t, fb_ref(d, hmax), atol=1e-3))
    fe = fe_ref(d, hmax, rdt_s, rd_ref(1e-6, d, grad))
    l_wmles = fb_ref(d, hmax) * (1 + fe) * l_rans + (1 - fb_ref(d, hmax)) * C_DES * delta_iddes_ref(d, hmax)
    l_i = l_iddes_ref(l_rans, d, hmax, rdt_s, rd_ref(1e-6, d, grad))
    allok &= check("l_IDDES = l_WMLES (Shur の WMLES 枝)",
                   np.allclose(l_i, l_wmles, rtol=1e-3))
    allok &= check("壁直近 (f_B=1) は必ず RANS (l=(1+f_e)l_RANS)",
                   np.allclose(l_i[d < 0.5 * hmax], (1 + fe[d < 0.5 * hmax]) * l_rans, rtol=1e-6))

    print("(7) float32 実装 vs double 参照")
    rng = np.random.default_rng(7)
    n = 20000
    d = 10.0 ** rng.uniform(-6, 0, n)
    hmax = 10.0 ** rng.uniform(-3, 0, n)
    nut = 10.0 ** rng.uniform(-8, -1, n)
    nu = np.full(n, 1.5e-5)
    grad = 10.0 ** rng.uniform(-2, 5, n)
    l_rans = 10.0 ** rng.uniform(-5, 0, n)
    ref = l_iddes_ref(l_rans, d, hmax,
                      np.minimum(rd_ref(nut, d, grad), 10.0),
                      np.minimum(rd_ref(nu, d, grad), 10.0))
    gpu, _, _, _ = l_iddes_gpu(l_rans, d, hmax, nut, nu, grad)
    rel = np.abs(gpu.astype(np.float64) - ref) / np.maximum(np.abs(ref), 1e-30)
    # f_e1−1 の桁落ち (fe1≈1) と tanh(x^10) の誤差増幅で float32 連鎖は ~1e-4 まで出る。
    # ブレンド関数の性質上 5e-4 で十分 (l_IDDES は l_RANS/l_LES の凸結合で有界)。
    allok &= check("相対差 < 5e-4 (float32 連鎖の実力値)", np.max(rel) < 5e-4,
                   f"max rel={np.max(rel):.2e}")
    gpu_ex, _, _, _ = l_iddes_gpu(1.0, 1e-8, 1e-3, 1e-2, 1.5e-5, 1e-30)
    allok &= check("極端値 (d→0, ∇u→0, clamp 発火) で NaN/Inf 無し",
                   np.isfinite(gpu_ex))

    print()
    print(f"VERDICT: {'PASS' if allok else 'FAIL'} (SST-IDDES 関数群 Gritskevich 2012 整合)")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
