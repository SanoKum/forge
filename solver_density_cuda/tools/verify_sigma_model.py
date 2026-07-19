#!/usr/bin/env python3
"""σ-model (Nicoud 2011) のカーネル実装式をコード前に数値検証する。

(1) 3×3 SPD (G=gᵀg) の固有値三角関数閉形式 vs numpy SVD の特異値 (乱数 10⁴ 個)
(2) 設計性質: 純せん断 / 剛体回転 / 2成分流 (TG 初期場含む) / 等方伸長で D_σ=0
(3) 乱流的な乱数勾配で D_σ>0
"""
import numpy as np


def sigma_closed_form(g):
    """カーネルに書く形の閉形式 (double)。返り値: (σ1, σ2, σ3) 降順。"""
    G = g.T @ g
    I1 = np.trace(G)
    I2 = 0.5*(I1*I1 - np.trace(G @ G))
    I3 = np.linalg.det(G)
    a1 = I1*I1/9.0 - I2/3.0
    a2 = I1*I1*I1/27.0 - I1*I2/6.0 + I3/2.0
    a1 = max(a1, 0.0)
    denom = a1**1.5
    arg = a2/denom if denom > 0 else 0.0
    arg = min(1.0, max(-1.0, arg))
    a3 = np.arccos(arg)/3.0
    s = np.sqrt(max(a1, 0.0))
    l1 = I1/3.0 + 2.0*s*np.cos(a3)
    l2 = I1/3.0 - 2.0*s*np.cos(np.pi/3.0 + a3)
    l3 = I1/3.0 - 2.0*s*np.cos(np.pi/3.0 - a3)
    lam = sorted([max(l1, 0.0), max(l2, 0.0), max(l3, 0.0)], reverse=True)
    return tuple(np.sqrt(v) for v in lam)


def D_sigma(g):
    s1, s2, s3 = sigma_closed_form(g)
    if s1 <= 0.0:
        return 0.0
    return s3*(s1-s2)*(s2-s3)/(s1*s1)


def main():
    rng = np.random.default_rng(5)
    # (1) 閉形式 vs SVD
    err = 0.0
    for _ in range(10000):
        g = rng.standard_normal((3, 3))
        ref = np.linalg.svd(g, compute_uv=False)
        got = sigma_closed_form(g)
        err = max(err, max(abs(a-b)/max(ref[0], 1e-30) for a, b in zip(got, ref)))
    # 縮退固有値近傍で三角関数式は O(1e-10) の丸めを持つ (既知)。粘性モデル用途には十分。
    print(f"(1) closed-form vs SVD: max rel err = {err:.3e} ->", "PASS" if err < 1e-8 else "FAIL")

    # (2) 設計性質
    cases = {
        "pure shear (g12=1)": np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0]], float),
        "solid rotation": np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], float),
        "isotropic strain": np.eye(3),
        "axisym strain (1,1,-2)/|..|": np.diag([1.0, 1.0, -2.0]),
    }
    # 2成分流 (TG 初期場: w=0, ∂/∂z 任意でも g の第3行=0)
    x, y, z = 0.7, 1.1, 0.4
    M0 = 0.4
    gtg = np.array([
        [M0*np.cos(x)*np.cos(y)*np.cos(z), -M0*np.sin(x)*np.sin(y)*np.cos(z), -M0*np.sin(x)*np.cos(y)*np.sin(z)],
        [M0*np.sin(x)*np.sin(y)*np.cos(z), -M0*np.cos(x)*np.cos(y)*np.cos(z),  M0*np.cos(x)*np.sin(y)*np.sin(z)],
        [0.0, 0.0, 0.0]])
    cases["TG initial field (2-component)"] = gtg
    ok = True
    for name, g in cases.items():
        d = D_sigma(g)
        good = d < 1e-12*max(1.0, np.abs(g).max()**3)
        ok &= good
        print(f"(2) {name:34s} D_sigma = {d:.3e} ->", "PASS (=0)" if good else "FAIL")

    # (3) 乱数勾配 (3D 乱流的) で正
    dmin = min(D_sigma(rng.standard_normal((3, 3))) for _ in range(1000))
    davg = np.mean([D_sigma(rng.standard_normal((3, 3))) for _ in range(1000)])
    print(f"(3) random 3D gradients: min={dmin:.3e} mean={davg:.3e} ->",
          "PASS (>=0, typically >0)" if dmin >= 0 and davg > 1e-3 else "FAIL")
    print("ALL PASS" if ok and err < 1e-8 else "CHECK FAILURES")


if __name__ == "__main__":
    main()
