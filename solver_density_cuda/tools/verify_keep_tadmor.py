#!/usr/bin/env python3
"""KEEP 中心流束の Tadmor 条件 (entropy conservation) 数値検証。

`convectiveFlux_keep_d.inc.cuh` の CPG 枝の中心流束 F* について
    E = Δwᵀ F* − Δψ    (w: エントロピー変数, ψ = ρ u_n: エントロピーポテンシャル流束)
を乱数状態対で評価する。E==0 なら entropy conservative (EC)、E≠0 なら非 EC。

検証内容 (すべて自動判定):
  1. ハーネス妥当性: 厳密 EC な Chandrashekar (2013) KEPEC 流束 (log-mean) で |E| が
     機械精度ゼロになること。
  2. KEEP は厳密 EC ではない: O(1) ランダムジャンプで |E| が有意に非ゼロ・両符号
     (E>0 = エントロピー生成側も出る) であること。
  3. 準保存性: 微小ジャンプ ε に対し |E| = O(ε³) (log-log 勾配 ≈ 3) であること。
     → 滑らかな場では設計次数でエントロピー整合、大ジャンプ (衝撃等) では非保存。

エントロピー対: η = −ρs/(γ−1), s = ln p − γ ln ρ, ψ_n = ρ u_n (Chandrashekar 2013)。

結論の使い方: KEEP は「KE 保存・エントロピー O(Δ³) 準保存」であり、コード/文書で
「エントロピー厳密保存」「散逸込みで entropy stable」と主張してはならない。
厳密 EC が必要になった場合は中心流束の log-mean 化 (KEPEC) が必要。

実行: python3 verify_keep_tadmor.py   (依存: numpy のみ)
"""
import numpy as np

ga = 1.4


def entropy_vars(ro, u, p):
    """w = ∂η/∂U = [(γ−s)/(γ−1) − β|u|², 2βu, −2β], β = ρ/(2p)"""
    q2 = np.sum(u*u, axis=-1)
    b = ro/(2.0*p)
    s = np.log(p) - ga*np.log(ro)
    w0 = (ga - s)/(ga - 1.0) - b*q2
    return np.stack([w0, 2*b*u[..., 0], 2*b*u[..., 1], 2*b*u[..., 2], -2*b], axis=-1)


def psi(ro, u, p, n):
    return ro*np.sum(u*n, axis=-1)


def keep_flux(ro0, u0, p0, ro1, u1, p1, n):
    """forge KEEP_d の !advGauge CPG 枝そのまま (pRef=0, 散逸なし, 面積 s=1)"""
    ub = 0.5*(u0+u1)
    C = 0.5*(ro0+ro1)*np.sum(ub*n, axis=-1)
    M = C[..., None]*ub
    K = C*0.5*np.sum(u0*u1, axis=-1)
    I = C*0.5*(p0/ro0 + p1/ro1)/(ga-1.0)
    G = 0.5*(p0+p1)[..., None]*n
    P = 0.5*np.sum((u0*p1[..., None] + u1*p0[..., None])*n, axis=-1)
    return np.concatenate([C[..., None], M+G, (K+I+P)[..., None]], axis=-1)


def logmean(a, b):
    z = a/b
    f = (z-1)/(z+1)
    u = f*f
    small = u < 1e-4
    F = np.where(small, 1 + u/3 + u*u/5 + u**3/7,
                 np.where(u > 0, np.log(z)/(2*f + 1e-300), 1.0))
    return (a+b)/(2*F)


def kepec_flux(ro0, u0, p0, ro1, u1, p1, n):
    """Chandrashekar (2013) KEPEC (厳密 EC) — ハーネス検証用リファレンス"""
    b0, b1 = ro0/(2*p0), ro1/(2*p1)
    rol = logmean(ro0, ro1)
    bl = logmean(b0, b1)
    rob = 0.5*(ro0+ro1)
    bb = 0.5*(b0+b1)
    ub = 0.5*(u0+u1)
    un = np.sum(ub*n, axis=-1)
    Fro = rol*un
    ptil = rob/(2*bb)
    Fm = Fro[..., None]*ub + ptil[..., None]*n
    q2b = 0.5*(np.sum(u0*u0, axis=-1) + np.sum(u1*u1, axis=-1))
    FE = (1/(2*(ga-1)*bl) - 0.5*q2b)*Fro + np.sum(ub*Fm, axis=-1)
    return np.concatenate([Fro[..., None], Fm, FE[..., None]], axis=-1)


def tadmor_E(flux, ro0, u0, p0, ro1, u1, p1, n):
    F = flux(ro0, u0, p0, ro1, u1, p1, n)
    dw = entropy_vars(ro1, u1, p1) - entropy_vars(ro0, u0, p0)
    dpsi = psi(ro1, u1, p1, n) - psi(ro0, u0, p0, n)
    return np.sum(dw*F, axis=-1) - dpsi


def main():
    rng = np.random.default_rng(0)
    N = 200000
    ro0 = rng.uniform(0.5, 2.0, N); ro1 = rng.uniform(0.5, 2.0, N)
    p0 = rng.uniform(0.5, 2.0, N);  p1 = rng.uniform(0.5, 2.0, N)
    u0 = rng.uniform(-1, 1, (N, 3)); u1 = rng.uniform(-1, 1, (N, 3))
    n = rng.normal(size=(N, 3)); n /= np.linalg.norm(n, axis=-1, keepdims=True)

    ok = True

    # 1. ハーネス妥当性 (KEPEC は機械精度で EC)
    E_kepec = tadmor_E(kepec_flux, ro0, u0, p0, ro1, u1, p1, n)
    e1 = np.max(np.abs(E_kepec))
    r1 = e1 < 1e-11
    ok &= r1
    print(f"[1] harness (KEPEC log-mean, 厳密 EC): max|E| = {e1:.3e} "
          f"(< 1e-11) ... {'PASS' if r1 else 'FAIL'}")

    # 2. KEEP は厳密 EC ではない (O(1) ジャンプで非ゼロ・両符号)
    E_keep = tadmor_E(keep_flux, ro0, u0, p0, ro1, u1, p1, n)
    e2 = np.max(np.abs(E_keep))
    fpos = np.mean(E_keep > 0)
    r2 = (e2 > 1e-3) and (0.05 < fpos < 0.95)
    ok &= r2
    print(f"[2] KEEP 非 EC 確認: max|E| = {e2:.3e} (≫0), E>0 割合 = {fpos*100:.1f}% "
          f"(両符号) ... {'PASS' if r2 else 'FAIL'}")

    # 3. 微小ジャンプ ε で |E| = O(ε³)
    roL, pL = 1.0, 1.0
    uL = np.array([0.3, -0.2, 0.1])
    d = np.array([0.7, -0.4, 0.5, 0.3, -0.6])
    nn = np.array([0.6, 0.48, 0.64]); nn /= np.linalg.norm(nn)
    eps_list = [1e-1, 1e-2, 1e-3, 1e-4]
    Es = []
    for eps in eps_list:
        roR = roL + eps*d[0]; uR = uL + eps*d[1:4]; pR = pL + eps*d[4]
        E = tadmor_E(keep_flux, np.array([roL]), uL[None, :], np.array([pL]),
                     np.array([roR]), uR[None, :], np.array([pR]), nn[None, :])
        Es.append(abs(E[0]))
    slope = np.polyfit(np.log(eps_list), np.log(Es), 1)[0]
    r3 = 2.7 < slope < 3.3
    ok &= r3
    print(f"[3] ε スケーリング: |E| = " +
          ", ".join(f"{e:.2e}@ε={x:.0e}" for e, x in zip(Es, eps_list)))
    print(f"    log-log 勾配 = {slope:.3f} (≈3 → O(ε³) 準保存) ... {'PASS' if r3 else 'FAIL'}")

    print()
    print(f"VERDICT: {'PASS' if ok else 'FAIL'} — KEEP 中心流束は厳密 EC ではなく "
          f"エントロピー O(Δ³) 準保存 (KE 保存は別途成立)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
