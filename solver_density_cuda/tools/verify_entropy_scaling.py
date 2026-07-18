#!/usr/bin/env python3
"""Step 2 実装前の式検証: CPG Euler の entropy-scaled 固有ベクトル。

検証項目:
 1. H = ∂U/∂w (解析形) が w = ∂η/∂U (数値微分) と整合 (η = -ρs/(γ-1), s = ln p - γ ln ρ)
 2. S = R⁻¹ H R⁻ᵀ が対角・正 → RSRᵀ = H (Barth scaling)
 3. S の対角が文献値 (acoustic: ρ/(2γ), entropy: (γ-1)ρ/γ, shear: p) と一致
 4. 小ジャンプで R|Λ|SRᵀ Δw ≈ R|Λ|R⁻¹ ΔU (Roe 型と漸近一致)
 5. Q = R|Λ|SRᵀ が対称半正定値
"""
import numpy as np
rng = np.random.default_rng(7)
ga = 1.4

def prim2U(ro, u, p):
    return np.array([ro, ro*u[0], ro*u[1], ro*u[2], p/(ga-1)+0.5*ro*(u@u)])

def U2prim(U):
    ro = U[0]; u = U[1:4]/ro; p = (ga-1)*(U[4]-0.5*ro*(u@u))
    return ro, u, p

def entropy_eta(U):
    ro, u, p = U2prim(U)
    s = np.log(p) - ga*np.log(ro)
    return -ro*s/(ga-1)

def w_of(U):
    ro, u, p = U2prim(U)
    s = np.log(p) - ga*np.log(ro)
    b = ro/(2*p)
    return np.array([(ga-s)/(ga-1) - b*(u@u), 2*b*u[0], 2*b*u[1], 2*b*u[2], -2*b])

def H_analytic(U):
    # H = ∂U/∂w (対称)。数値微分で作る (解析形のরefではなく検証目的なので数値で十分)
    H = np.zeros((5,5))
    w0 = w_of(U)
    # invert numerically: perturb w -> U via Newton (or use dU/dw = (dw/dU)^-1)
    dwdU = np.zeros((5,5))
    eps = 1e-7
    for j in range(5):
        Up = U.copy(); Um = U.copy()
        h = eps*max(abs(U[j]), 1.0)
        Up[j] += h; Um[j] -= h
        dwdU[:, j] = (w_of(Up) - w_of(Um))/(2*h)
    return np.linalg.inv(dwdU)

def eig_RS(U, n, t1, t2):
    ro, u, p = U2prim(U)
    c = np.sqrt(ga*p/ro)
    un = u@n; Ht = (U[4]+p)/ro
    R = np.zeros((5,5))
    R[:,0] = np.r_[1.0, u - c*n, Ht - c*un]          # un - c
    R[:,1] = np.r_[1.0, u,        0.5*(u@u)]          # un (entropy)
    R[:,2] = np.r_[0.0, t1,       u@t1]               # un (shear1)
    R[:,3] = np.r_[0.0, t2,       u@t2]               # un (shear2)
    R[:,4] = np.r_[1.0, u + c*n,  Ht + c*un]          # un + c
    lam = np.array([un-c, un, un, un, un+c])
    S_lit = np.diag([ro/(2*ga), ro*(ga-1)/ga, p, p, ro/(2*ga)])
    return R, lam, S_lit, c

ok = True
for trial in range(5):
    ro = rng.uniform(0.3, 3.0); p = rng.uniform(0.3, 3.0)
    u = rng.uniform(-1, 1, 3)
    n = rng.normal(size=3); n /= np.linalg.norm(n)
    a = np.array([1.0,0,0]) if abs(n[0]) < 0.9 else np.array([0,1.0,0])
    t1 = np.cross(n, a); t1 /= np.linalg.norm(t1); t2 = np.cross(n, t1)
    U = prim2U(ro, u, p)

    # 1. w = ∂η/∂U 数値微分照合
    grad = np.zeros(5); eps = 1e-6
    for j in range(5):
        Up = U.copy(); Um = U.copy(); h = eps*max(abs(U[j]),1.0)
        Up[j]+=h; Um[j]-=h
        grad[j] = (entropy_eta(Up)-entropy_eta(Um))/(2*h)
    e1 = np.max(np.abs(grad - w_of(U)))

    # 2-3. S = R⁻¹HR⁻ᵀ が対角・文献値一致
    H = H_analytic(U)
    R, lam, S_lit, c = eig_RS(U, n, t1, t2)
    S_num = np.linalg.inv(R) @ H @ np.linalg.inv(R).T
    offdiag = np.max(np.abs(S_num - np.diag(np.diag(S_num))))
    e3 = np.max(np.abs(np.diag(S_num) - np.diag(S_lit)))

    # 4. 小ジャンプ整合: R|Λ|S Rᵀ Δw ≈ R|Λ|R⁻¹ ΔU
    dU = 1e-5*rng.normal(size=5)*np.abs(U)
    U2 = U + dU
    dw = w_of(U2) - w_of(U)
    D_es  = R @ np.diag(np.abs(lam)) @ S_lit @ R.T @ dw
    D_roe = R @ np.diag(np.abs(lam)) @ np.linalg.inv(R) @ dU
    e4 = np.max(np.abs(D_es - D_roe))/max(np.max(np.abs(D_roe)), 1e-30)

    # 5. Q SPD
    Q = R @ np.diag(np.abs(lam)) @ S_lit @ R.T
    ev = np.linalg.eigvalsh(0.5*(Q+Q.T))
    e5 = ev.min()

    print(f"trial{trial}: |w-∂η/∂U|={e1:.2e}  S offdiag={offdiag:.2e}  |S-lit|={e3:.2e}  "
          f"ES↔Roe rel={e4:.2e}  minEig(Q)={e5:+.2e}")
    ok &= e1 < 1e-5 and offdiag < 1e-5 and e3 < 1e-5 and e4 < 1e-3 and e5 > -1e-10

print("VERDICT:", "PASS" if ok else "FAIL")
