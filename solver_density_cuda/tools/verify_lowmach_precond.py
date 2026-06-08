#!/usr/bin/env python3
"""Phase 4 (完全 Γ⁻¹A 低マッハ前処理) の閉形を検証するホスト側ユニット検査。

device 実装 (`cuda_forge/lowMachPrecond_d.cuh` の lowMachGammaC / lowMachGammaCinv、
`cuda_forge/timeIntegration_d.cu` の前処理ヤコビアン) が依拠する数式を独立に確認する:

  1. Γ_c = Γ_p M⁻¹ = I + ((1-β)/(βc²)) g rᵀ        (ランク 1 閉形)
  2. Γ_c⁻¹ = I - ((1-β)/c²) g rᵀ                     (Sherman-Morrison)
  3. rᵀg = c²
  4. β=1 (超音速) で Γ_c = Γ_c⁻¹ = I
  5. 前処理ヤコビアン P = Γ_c⁻¹ A_c の固有値 = {Uₙ,Uₙ,Uₙ, ½(1+β)Uₙ ± c'}

理論は docs/time_integration/theory.md「低マッハ前処理固有系 (Weiss-Smith)」節、
計画は .github/plans/time_integration-lowmach-preconditioning.md §5 Phase 4。
依存は numpy のみ。`python3 tools/verify_lowmach_precond.py` で実行。
"""
import numpy as np

GA = 1.4
R = 287.0


def flux_normal(Q, n):
    """Euler 法線フラックス F_n(Q), Q=(ρ,ρu,ρv,ρw,ρE)。"""
    rho, mx, my, mz, rhoE = Q
    u, v, w = mx / rho, my / rho, mz / rho
    Un = u * n[0] + v * n[1] + w * n[2]
    p = (GA - 1.0) * (rhoE - 0.5 * (mx * u + my * v + mz * w))
    return np.array([rho * Un,
                     mx * Un + p * n[0],
                     my * Un + p * n[1],
                     mz * Un + p * n[2],
                     (rhoE + p) * Un])


def jac_fd(Q, n, eps=1e-4):
    """A_c = ∂F_n/∂Q を中心差分で。"""
    A = np.zeros((5, 5))
    for j in range(5):
        dQ = np.zeros(5); dQ[j] = eps * max(abs(Q[j]), 1.0)
        A[:, j] = (flux_normal(Q + dQ, n) - flux_normal(Q - dQ, n)) / (2 * dQ[j])
    return A


def gamma_c(gamma, c, u, v, w, beta, coeff):
    gm1 = gamma - 1.0
    ek = 0.5 * (u * u + v * v + w * w)
    H = c * c / gm1 + ek
    g = np.array([1.0, u, v, w, H])
    r = gm1 * np.array([ek, -u, -v, -w, 1.0])
    return np.eye(5) + coeff * np.outer(g, r), g, r


def main():
    rng = np.random.default_rng(2)
    worst = {"gc": 0.0, "inv": 0.0, "rg": 0.0, "b1": 0.0, "eig": 0.0}
    for _ in range(2000):
        u, v, w = rng.normal(scale=250, size=3)
        T = 250 + rng.random() * 600
        p = 1e5 * (0.3 + rng.random() * 1.5)
        rho = p / (R * T)
        c = np.sqrt(GA * R * T)
        ek = 0.5 * (u * u + v * v + w * w)
        rhoE = p / (GA - 1.0) + rho * ek
        Q = np.array([rho, rho * u, rho * v, rho * w, rhoE])
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        Un = u * n[0] + v * n[1] + w * n[2]

        for beta in (rng.random() * 0.97 + 0.02, 1.0):
            Ur2 = beta * c * c
            # --- 参照 M, Γ_p, Γ_c=Γ_p M⁻¹ ---
            H = c * c / (GA - 1.0) + ek
            rho_p, rho_T, cv = rho / p, -rho / T, R / (GA - 1.0)
            M = np.zeros((5, 5))
            M[:, 0] = [rho_p, rho_p * u, rho_p * v, rho_p * w, rho_p * H - 1.0]
            M[:, 1] = [0, rho, 0, 0, rho * u]
            M[:, 2] = [0, 0, rho, 0, rho * v]
            M[:, 3] = [0, 0, 0, rho, rho * w]
            M[:, 4] = [rho_T, rho_T * u, rho_T * v, rho_T * w,
                       rho_T * H - rho_T * p / rho + rho * cv]
            Theta = 1.0 / Ur2 + (GA - 1.0) / (c * c)
            Gp = M.copy(); Gp[:, 0] = [Theta, Theta * u, Theta * v, Theta * w, Theta * H - 1.0]
            Gc_ref = Gp @ np.linalg.inv(M)
            # --- ランク 1 閉形 (device 実装と同じ) ---
            Gc, g, r = gamma_c(GA, c, u, v, w, beta, (1 - beta) / (beta * c * c))
            Gci, _, _ = gamma_c(GA, c, u, v, w, beta, -(1 - beta) / (c * c))

            worst["gc"] = max(worst["gc"], np.abs(Gc - Gc_ref).max() / np.abs(Gc_ref).max())
            worst["inv"] = max(worst["inv"], np.abs(Gc @ Gci - np.eye(5)).max())
            worst["rg"] = max(worst["rg"], abs(r @ g - c * c) / (c * c))
            if beta == 1.0:
                worst["b1"] = max(worst["b1"], np.abs(Gc - np.eye(5)).max(),
                                  np.abs(Gci - np.eye(5)).max())

            # --- 前処理固有値 ---
            A = jac_fd(Q, n)
            P = Gci @ A
            eig = np.sort(np.real(np.linalg.eigvals(P)))
            cp = 0.5 * np.sqrt((1 - beta) ** 2 * Un * Un + 4 * Ur2)
            ref = np.sort([Un, Un, Un, 0.5 * (1 + beta) * Un + cp, 0.5 * (1 + beta) * Un - cp])
            worst["eig"] = max(worst["eig"], np.abs(eig - ref).max() / (abs(c) + abs(Un)))

    print("Phase 4 低マッハ前処理 閉形検証 (2000 ランダム状態, 相対誤差の最大):")
    print(f"  1. Γ_c == Γ_p M⁻¹ (ランク1)      : {worst['gc']:.2e}")
    print(f"  2. Γ_c · Γ_c⁻¹ == I              : {worst['inv']:.2e}")
    print(f"  3. rᵀg == c²                     : {worst['rg']:.2e}")
    print(f"  4. β=1 で Γ_c=Γ_c⁻¹=I            : {worst['b1']:.2e}")
    print(f"  5. eig(Γ_c⁻¹A) == λ' (前処理固有値): {worst['eig']:.2e}")
    tol = 1e-3  # 固有値は FD ヤコビアン由来で緩め
    ok = all(v < tol for v in worst.values())
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
