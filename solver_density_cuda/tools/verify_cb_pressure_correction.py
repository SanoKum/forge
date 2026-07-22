#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""高周波圧力欠陥駆動 mass-flux 補正 (keepDissCbCoeff) の理論検証。

plan: plans/active/convection-keep-cb-pressure-correction.md §4.5

補正 flux: δF_f = −k_f δp^HF g_f,  g_f=[1,uF,vF,wF,HtF]ᵀ, k_f=½C_cb S_f/Ur > 0
entropy 寄与: Δwᵀ δF_f = −k_f δp^HF r_f,  r_f = Δwᵀ g_f
sign gate: δp^HF·r_f > 0 の面のみ有効 → 寄与は −k_f·(正) ≤ 0 (散逸的) を構造保証。

検証項目 (レビュー反映):
  ① gate なしの符号分布 (δp·r_f < 0 の頻度と状態依存性)
  ② gate 後の負寄与 (=反散逸) が厳密ゼロ
  ③ 純圧力市松が gate で削られない (素通し率 100%)
  ④ gate 発火率 (乱流的ランダム状態でどの程度の面が落ちるか)
  ⑤ 2 セル市松の減衰方向: res 符号が高圧側→低圧側の質量輸送になっている
"""
import numpy as np

ga = 1.4
rng = np.random.default_rng(20260723)


def entropy_vars(ro, u, v, w, p):
    """w = [(γ−s)/(γ−1) − β|u|², 2βu, 2βv, 2βw, −2β],  β=ρ/(2p), s=ln p − γ ln ρ"""
    b = ro / (2.0 * p)
    q = u * u + v * v + w * w
    s = np.log(p) - ga * np.log(ro)
    return np.array([(ga - s) / (ga - 1.0) - b * q, 2 * b * u, 2 * b * v, 2 * b * w, -2 * b])


def face_quantities(st0, st1):
    """算術平均の g_f=[1,uF,vF,wF,HtF] と Δw, 生 Δp を返す (カーネルと同じ組み方)"""
    ro0, u0, v0, w0, p0 = st0
    ro1, u1, v1, w1, p1 = st1
    roF = 0.5 * (ro0 + ro1)
    PsF = 0.5 * (p0 + p1)
    uF, vF, wF = 0.5 * (u0 + u1), 0.5 * (v0 + v1), 0.5 * (w0 + w1)
    qF = uF * uF + vF * vF + wF * wF
    Ht = ga * PsF / ((ga - 1.0) * roF) + 0.5 * qF
    g = np.array([1.0, uF, vF, wF, Ht])
    dw = entropy_vars(ro1, u1, v1, w1, p1) - entropy_vars(ro0, u0, v0, w0, p0)
    return g, dw, p1 - p0


def sample_state(mach_scale=0.3):
    ro = 10 ** rng.uniform(-1, 1)
    p = 10 ** rng.uniform(4, 6)
    c = np.sqrt(ga * p / ro)
    vel = rng.normal(0, mach_scale * c / np.sqrt(3), 3)
    return ro, *vel, p


def perturb(st, amp):
    ro, u, v, w, p = st
    c = np.sqrt(ga * p / ro)
    # 正値性ガード: カーネル側も pL,pR>0 の面しか補正しないため負圧サンプルは対象外
    return (ro * max(1 + rng.normal(0, amp), 0.05),
            u + rng.normal(0, amp * c), v + rng.normal(0, amp * c), w + rng.normal(0, amp * c),
            p * max(1 + rng.normal(0, amp), 0.05))


def main():
    ok = True
    N = 200000

    # ① / ④ ランダム状態対での符号分布と gate 発火率
    for amp, label in [(0.01, "smooth-ish (1% jump)"), (0.1, "turbulent (10% jump)"),
                       (0.5, "extreme (50% jump)")]:
        neg = 0
        used = 0
        for _ in range(N // 10):
            st0 = sample_state()
            st1 = perturb(st0, amp)
            g, dw, dp = face_quantities(st0, st1)
            r = dw @ g
            if dp * r < 0:
                neg += 1
            else:
                used += 1
        rate = neg / (neg + used)
        print(f"[①④] {label:24s}: δp·r<0 (gate 落ち) = {rate*100:5.1f}%  素通し = {(1-rate)*100:5.1f}%")

    # ② gate 後の entropy 寄与が非正 (=散逸的) を数値で確認
    worst = -np.inf
    for _ in range(N // 10):
        st0 = sample_state()
        st1 = perturb(st0, rng.uniform(0.01, 0.5))
        g, dw, dp = face_quantities(st0, st1)
        r = dw @ g
        dp_used = dp if dp * r > 0 else 0.0
        contrib = -dp_used * r          # ∝ Δwᵀ δF (k_f>0 は省略)
        worst = max(worst, contrib)
    print(f"[②] gate 後の最大 entropy 寄与 (負=散逸): {worst:.3e}  ", end="")
    if worst <= 0.0:
        print("PASS (厳密非正)")
    else:
        print("FAIL"); ok = False

    # ③ 純圧力市松 (ρ,u 一定・p のみ振動) が素通しされるか
    passed = 0
    total = 0
    for _ in range(2000):
        ro = 10 ** rng.uniform(-1, 1)
        p0 = 10 ** rng.uniform(4, 6)
        dp_amp = p0 * rng.uniform(0.001, 0.3)
        vel = rng.normal(0, 0.3 * np.sqrt(ga * p0 / ro) / np.sqrt(3), 3)
        st0 = (ro, *vel, p0 - 0.5 * dp_amp)
        st1 = (ro, *vel, p0 + 0.5 * dp_amp)
        g, dw, dp = face_quantities(st0, st1)
        r = dw @ g
        total += 1
        if dp * r > 0:
            passed += 1
    print(f"[③] 純圧力市松の素通し率: {passed}/{total}  ", end="")
    if passed == total:
        print("PASS")
    else:
        print("FAIL"); ok = False

    # ⑤ 2 セル市松の減衰方向: flux F_cb = −k δp g の質量成分は δp>0 (ic1 高圧) で負
    #    = ic1→ic0 (高圧→低圧) の質量輸送。ρ,u 同一・p1>p0 の面で符号確認。
    st0 = (1.2, 10.0, 0.0, 0.0, 100000.0)
    st1 = (1.2, 10.0, 0.0, 0.0, 102000.0)
    g, dw, dp = face_quantities(st0, st1)
    r = dw @ g
    F_mass = -dp * g[0]  # k_f>0 省略
    print(f"[⑤] p1>p0 で mass flux 符号 = {np.sign(F_mass):+.0f} (期待 −1: 高圧→低圧)  "
          f"gate: δp·r = {dp*r:.3e} (>0 で素通し)  ", end="")
    if F_mass < 0 and dp * r > 0:
        print("PASS")
    else:
        print("FAIL"); ok = False

    print("\nVERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
