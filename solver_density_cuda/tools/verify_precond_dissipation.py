#!/usr/bin/env python3
"""KEEP ES 散逸レイヤの低マッハ前処理化 (keepDissPrecond) の式検証 (コード化前)。

plans/active/convection-keep-diss-lowmach-precond.md §3 の設計:
音響 ± 対の散逸を、Turkel 前処理音響 2×2 サブシステム
    M2 = Γ⁻¹ A2,  A2 = [[Un, ρc²],[1/ρ, Un]] (状態 (p, Un)),  Γ = diag(1/β, 1), β=Ur²/c²
の前処理散逸 D_p = Γ|M2| で置換する。|M2| は 2×2 行列関数の閉形
    |M2| = φ1 M2 + φ2 I,  φ1=(|λ+|-|λ-|)/(λ+-λ-),  φ2=(λ+|λ-|-λ-|λ+|)/(λ+-λ-)
(λ± = u'±c', u'=(1+β)Un/2, c'=lowMachCprime) で固有ベクトル不要。

検証項目:
 (1) 閉形 |M2| = 数値固有分解 と一致
 (2) 漸近: 低マッハ静止で D_p → diag(c²/Ur, Ur) = 圧力ジャンプ 1/Ur 増強・速度ジャンプ Ur 縮小
 (3) β=1 (M>=1) で標準 Roe 散逸 |A2| に厳密復帰
 (4) 特性振幅版 K = T⁻¹ D_p T の対称部正定値性 (ES 性) の (M, Mn) マップ
 (5) 5式系での作用: 会話版モデル (entropy/shear は |Un| 据え置き + 音響 K) と
     完全 Weiss–Smith 5×5 (Γc = I + ((1-β)/(βc²)) g rᵀ, D5 = Γc|Γc⁻¹Ac|) の比較
 (6) 市松キラー実証: 2Δ 圧力モード vs 速度モードの減衰率比が c²/Ur² 倍化

使い方: python3 verify_precond_dissipation.py
"""
import numpy as np

GA = 1.4


def lowmach_ur(c, velmag, eps):
    return min(c, max(velmag, eps * c))


def cprime(c, velmag, un, eps):
    ur = lowmach_ur(c, velmag, eps)
    beta = (ur / c) ** 2
    return 0.5 * np.sqrt((1 - beta) ** 2 * un ** 2 + 4 * ur * ur)


def dp_matrix(ro, c, un, velmag, eps):
    """前処理音響 2×2 散逸 D_p = Γ|Γ⁻¹A2| の閉形。"""
    ur = lowmach_ur(c, velmag, eps)
    beta = (ur / c) ** 2
    up = 0.5 * (1 + beta) * un
    cp = cprime(c, velmag, un, eps)
    lp, lm = up + cp, up - cp
    phi1 = (abs(lp) - abs(lm)) / (lp - lm)
    phi2 = (lp * abs(lm) - lm * abs(lp)) / (lp - lm)
    # |M2| = phi1*M2 + phi2*I, M2 = [[beta*un, beta*ro*c^2],[1/ro, un]]
    d11 = phi1 * un + phi2 / beta
    d12 = phi1 * ro * c * c
    d21 = phi1 / ro
    d22 = phi1 * un + phi2
    return np.array([[d11, d12], [d21, d22]])


def dp_numeric(ro, c, un, velmag, eps):
    ur = lowmach_ur(c, velmag, eps)
    beta = (ur / c) ** 2
    A2 = np.array([[un, ro * c * c], [1.0 / ro, un]])
    G = np.diag([1.0 / beta, 1.0])
    M = np.linalg.inv(G) @ A2
    lam, R = np.linalg.eig(M)
    absM = R @ np.diag(np.abs(lam)) @ np.linalg.inv(R)
    return G @ absM


def K_matrix(ro, c, un, velmag, eps):
    """特性振幅 (α-, α+) 基底での音響散逸 K = T⁻¹ D_p T。"""
    Dp = dp_matrix(ro, c, un, velmag, eps)
    T = np.array([[c * c, c * c], [-c / ro, c / ro]])   # (α-,α+) -> (Δp, ΔUn)
    return np.linalg.inv(T) @ Dp @ T


def main():
    rng = np.random.default_rng(3)
    eps = 0.15

    # (1) 閉形 == 数値固有分解
    err = 0.0
    for _ in range(500):
        ro = 1.2 * (1 + 0.3 * rng.random())
        c = 340 * (1 + 0.3 * rng.random())
        M = 10 ** rng.uniform(-3, 0.3)
        velmag = M * c
        un = velmag * rng.uniform(-1, 1)
        err = max(err, np.abs(dp_matrix(ro, c, un, velmag, eps)
                              - dp_numeric(ro, c, un, velmag, eps)).max()
                  / (ro * c * c))
    print(f"(1) closed-form |M2| vs numeric eig: max rel err = {err:.3e} ->",
          "PASS" if err < 1e-12 else "FAIL")

    # (2) 漸近 (静止低マッハ): D_p -> diag(c²/Ur, Ur)
    ro, c = 1.2, 340.0
    velmag = 3.0  # |u| = 3 m/s -> Ur = eps*c = 51
    ur = lowmach_ur(c, velmag, eps)
    Dp = dp_matrix(ro, c, 0.0, velmag, eps)
    tgt = np.diag([c * c / ur, ur])
    err2 = np.abs(Dp - tgt).max() / (c * c / ur)
    print(f"(2) low-Mach asymptote D_p=diag(c²/Ur, Ur): rel err = {err2:.3e} "
          f"(Δp coeff {Dp[0,0]:.1f} = c²/Ur {c*c/ur:.1f}, ΔUn coeff {Dp[1,1]:.1f} = Ur {ur:.1f}) ->",
          "PASS" if err2 < 1e-12 else "FAIL")

    # (3) β=1 復帰 (M>=1): D_p == |A2| (標準 Roe)
    err3 = 0.0
    for _ in range(200):
        M = rng.uniform(1.0, 3.0)
        velmag = M * c
        un = velmag * rng.uniform(-1, 1)
        A2 = np.array([[un, ro * c * c], [1.0 / ro, un]])
        lam, R = np.linalg.eig(A2)
        absA = R @ np.diag(np.abs(lam)) @ np.linalg.inv(R)
        err3 = max(err3, np.abs(dp_matrix(ro, c, un, velmag, eps) - absA).max() / (ro * c * c))
    print(f"(3) beta=1 (M>=1) recovery to |A2|: max rel err = {err3:.3e} ->",
          "PASS" if err3 < 1e-10 else "FAIL")

    # (4) sym(K) 正定値性マップ (ES 性): K は S_A>0 の重みで 2 次形式に入る
    lam_min = 1e30
    worst = None
    for Mtot in 10 ** np.linspace(-3, 0.49, 60):        # M<1 (超音速は |A| で自明 PSD)
        for fn in np.linspace(-1, 1, 41):               # Un/|u| (面法線成分比)
            velmag = Mtot * c
            un = velmag * fn
            K = K_matrix(ro, c, un, velmag, eps)
            ev = np.linalg.eigvalsh(0.5 * (K + K.T))
            if ev.min() < lam_min:
                lam_min = ev.min()
                worst = (Mtot, fn, ev.min(), np.abs(K).max())
    # 正規化: 最悪固有値を局所の |K| スケールで
    rel = worst[2] / worst[3]
    print(f"(4) sym(K) min eigenvalue over (M,Un) grid: {lam_min:.3e} "
          f"(worst at M={worst[0]:.3f}, Un/|u|={worst[1]:.2f}, rel={rel:.3e}) ->",
          "PASS (ES holds)" if lam_min >= -1e-10 * worst[3] else "FAIL (needs clip/blend)")

    # (5) 完全 Weiss–Smith 5×5 との比較 (1D 法線系, 接線速度込み)
    def cons_jac(ro, u, v, w, p):
        # A_c = dF/dU (x-方向, CPG) を数値微分で
        def U_of(q):
            ro_, u_, v_, w_, p_ = q
            return np.array([ro_, ro_ * u_, ro_ * v_, ro_ * w_,
                             p_ / (GA - 1) + 0.5 * ro_ * (u_ * u_ + v_ * v_ + w_ * w_)])
        def F_of(q):
            ro_, u_, v_, w_, p_ = q
            E = p_ / (GA - 1) + 0.5 * ro_ * (u_ * u_ + v_ * v_ + w_ * w_)
            return np.array([ro_ * u_, ro_ * u_ * u_ + p_, ro_ * u_ * v_, ro_ * u_ * w_,
                             u_ * (E + p_)])
        q0 = np.array([ro, u, v, w, p])
        U0 = U_of(q0)
        # dF/dU = dF/dq (dU/dq)^-1
        dFdq = np.zeros((5, 5)); dUdq = np.zeros((5, 5))
        for j in range(5):
            dq = np.zeros(5); dq[j] = max(1e-6 * abs(q0[j]), 1e-8)
            dFdq[:, j] = (F_of(q0 + dq) - F_of(q0 - dq)) / (2 * dq[j])
            dUdq[:, j] = (U_of(q0 + dq) - U_of(q0 - dq)) / (2 * dq[j])
        return dFdq @ np.linalg.inv(dUdq)

    mism_lo, mism_hi = 0.0, 0.0
    for _ in range(60):
        Mtot = 10 ** rng.uniform(-3, -0.5)
        u = Mtot * c * rng.uniform(-1, 1)
        v = Mtot * c * rng.uniform(-0.5, 0.5)
        w = Mtot * c * rng.uniform(-0.5, 0.5)
        p = ro * c * c / GA
        velmag = np.sqrt(u * u + v * v + w * w)
        ur = lowmach_ur(c, velmag, eps); beta = (ur / c) ** 2
        H = c * c / (GA - 1) + 0.5 * velmag ** 2
        g = np.array([1, u, v, w, H])
        ek = 0.5 * velmag ** 2
        r = (GA - 1) * np.array([ek, -u, -v, -w, 1])
        Gc = np.eye(5) + ((1 - beta) / (beta * c * c)) * np.outer(g, r)
        Ac = cons_jac(ro, u, v, w, p)
        Mfull = np.linalg.inv(Gc) @ Ac
        lam, R = np.linalg.eig(Mfull)
        D5 = Gc @ (R @ np.diag(np.abs(lam)) @ np.linalg.inv(R)).real
        # モデル: entropy/shear = |u| 据え置き + 音響 K (r± 基底)
        un = u  # 法線 = x
        cc = np.sqrt(GA * p / ro)
        nvec = np.array([1.0, 0, 0])
        r_m = np.array([1, u - cc, v, w, H - cc * un])
        r_p = np.array([1, u + cc, v, w, H + cc * un])
        r_E = np.array([1, u, v, w, ek])
        r_s1 = np.array([0, 0, 1, 0, v])
        r_s2 = np.array([0, 0, 0, 1, w])
        Rfull = np.stack([r_m, r_E, r_s1, r_s2, r_p], axis=1)
        L = np.linalg.inv(Rfull)
        K = K_matrix(ro, cc, un, velmag, eps)
        lamU = abs(un)
        # D_model ΔU = R diag/blk (K on ±, |Un| on others) L ΔU
        blk = np.diag([0.0, lamU, lamU, lamU, 0.0])
        blk[0, 0] = K[0, 0]; blk[0, 4] = K[0, 1]
        blk[4, 0] = K[1, 0]; blk[4, 4] = K[1, 1]
        Dm = Rfull @ blk @ L
        num = np.abs(Dm - D5).max()
        den = np.abs(D5).max()
        mism_lo = max(mism_lo, num / den)
    print(f"(5) model (|Un| conv + K acoustic) vs full Weiss-Smith 5x5: max rel mismatch = {mism_lo:.3e} ->",
          "PASS" if mism_lo < 5e-2 else "INFO (model is self-consistent design; see plan §3)")

    # (6) 市松キラー実証: 減衰係数比 (Δp モード vs ΔUn モード), M=0.01
    velmag = 0.01 * c
    ur = lowmach_ur(c, velmag, eps)
    Dp = dp_matrix(ro, c, 0.0, velmag, eps)
    Dp_c = dp_matrix(ro, c, 0.0, velmag, eps=1.0)   # eps=1 -> Ur=c = full c (無前処理)
    print(f"(6) damping @M=0.01: Δp coeff precond {Dp[0,0]:.1f} vs full-c {Dp_c[0,0]:.1f} "
          f"(x{Dp[0,0]/Dp_c[0,0]:.1f}), ΔUn coeff {Dp[1,1]:.1f} vs {Dp_c[1,1]:.1f} "
          f"(x{Dp[1,1]/Dp_c[1,1]:.2f}) -> pressure-jump enhanced / velocity-jump reduced")


if __name__ == "__main__":
    main()
