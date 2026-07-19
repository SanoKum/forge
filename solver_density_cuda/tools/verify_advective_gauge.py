#!/usr/bin/env python3
"""KEEP 移流項の基準差分ゲージ (uRef/roRef) の数値検証 (コード化前の式検証)。

plans/active/convection-freestream-preserving-flux.md §8.2 の因数分解が
 (1) 解析的に「元の KEEP 流束 − 参照一様流束 F_inf(s)」と一致すること (float64)
 (2) 一様場 (ρ,u,p)=(ρ∞,u∞,p∞) で float32 でも全項ビット単位ゼロになること
 (3) wavy メッシュのセル和で、非一様場でもゲージ有無の残差差が丸め誤差レベルであること
を確認する。

使い方: python3 verify_advective_gauge.py [wavy.h5]
"""
import sys
import numpy as np

GA = 1.4
F32 = np.float32


def keep_flux(ro0, u0, p0, ro1, u1, p1, n, pRef, dtype):
    """元の KEEP 中心流束 (単位面積, Gtilde は pRef 差分済) を dtype 演算で返す。"""
    f = dtype
    half = f(0.5)
    Ctilde = half*(ro0+ro1) * (half*(u0[0]+u1[0])*n[0]
                               + half*(u0[1]+u1[1])*n[1]
                               + half*(u0[2]+u1[2])*n[2])
    Mt = np.array([Ctilde*half*(u0[i]+u1[i]) for i in range(3)], dtype=f)
    Kt = Ctilde*half*(u0[0]*u1[0]+u0[1]*u1[1]+u0[2]*u1[2])
    It = Ctilde*half*(p0/ro0 + p1/ro1)/f(GA-1.0)
    Gt = np.array([half*((p0-pRef)+(p1-pRef))*n[i] for i in range(3)], dtype=f)
    Pt = half*((u0[0]*p1+u1[0]*p0)*n[0]+(u0[1]*p1+u1[1]*p0)*n[1]+(u0[2]*p1+u1[2]*p0)*n[2])
    return Ctilde, Mt+Gt, Kt+It+Pt


def keep_flux_gauged(ro0, u0, p0, ro1, u1, p1, n, pRef, roRef, uRef, dtype):
    """§8.2 の差分形因数分解 (ゲージ済) 流束。"""
    f = dtype
    half = f(0.5)
    # 差分 (小さい量を先に作る)
    dro = half*((ro0-roRef)+(ro1-roRef))
    du = np.array([half*((u0[i]-uRef[i])+(u1[i]-uRef[i])) for i in range(3)], dtype=f)
    ubar = np.array([half*(u0[i]+u1[i]) for i in range(3)], dtype=f)
    Unbar = ubar[0]*n[0]+ubar[1]*n[1]+ubar[2]*n[2]
    dUn = du[0]*n[0]+du[1]*n[1]+du[2]*n[2]
    UnRef = uRef[0]*n[0]+uRef[1]*n[1]+uRef[2]*n[2]
    Cinf = roRef*UnRef
    dC = dro*Unbar + roRef*dUn
    # 運動量 (Gtilde は従来どおり)
    Mt = np.array([dC*ubar[i] + Cinf*du[i] for i in range(3)], dtype=f)
    Gt = np.array([half*((p0-pRef)+(p1-pRef))*n[i] for i in range(3)], dtype=f)
    # エネルギー
    kbar = half*(u0[0]*u1[0]+u0[1]*u1[1]+u0[2]*u1[2])
    dK = dC*kbar + Cinf*half*((u0[0]-uRef[0])*u1[0]+(u0[1]-uRef[1])*u1[1]+(u0[2]-uRef[2])*u1[2]
                              + uRef[0]*(u1[0]-uRef[0])+uRef[1]*(u1[1]-uRef[1])+uRef[2]*(u1[2]-uRef[2]))
    eRefG = pRef/roRef  # (γ-1)e∞: カーネル内同一 float 除算
    de = half*((p0/ro0-eRefG)+(p1/ro1-eRefG))/f(GA-1.0)
    ebar = half*(p0/ro0+p1/ro1)/f(GA-1.0)
    dI = dC*ebar + Cinf*de
    dP = half*(((u0[0]-uRef[0])*p1+uRef[0]*(p1-pRef) + (u1[0]-uRef[0])*p0+uRef[0]*(p0-pRef))*n[0]
               + ((u0[1]-uRef[1])*p1+uRef[1]*(p1-pRef) + (u1[1]-uRef[1])*p0+uRef[1]*(p0-pRef))*n[1]
               + ((u0[2]-uRef[2])*p1+uRef[2]*(p1-pRef) + (u1[2]-uRef[2])*p0+uRef[2]*(p0-pRef))*n[2])
    return dC, Mt+Gt, dK+dI+dP


def ref_flux(roRef, uRef, pRef, n, dtype):
    """参照一様流束 F_inf(s)/S (単位面積, 圧力項は pRef ゲージで 0)。"""
    f = dtype
    UnRef = uRef[0]*n[0]+uRef[1]*n[1]+uRef[2]*n[2]
    Cinf = roRef*UnRef
    M = np.array([Cinf*uRef[i] for i in range(3)], dtype=f)
    k = f(0.5)*(uRef[0]**2+uRef[1]**2+uRef[2]**2)
    e = pRef/roRef/f(GA-1.0)
    E = Cinf*(k+e) + pRef*UnRef
    return Cinf, M, E


def main():
    rng = np.random.default_rng(7)
    roRef, pRef = 1.177, 101325.0
    c = np.sqrt(GA*pRef/roRef)
    uRef = np.array([0.1*c, 0.0, 0.0])

    # (1) 解析等価性 (float64, ランダム状態 × ランダム法線)
    err = 0.0
    for _ in range(2000):
        ro0, ro1 = roRef*(1+0.1*rng.standard_normal(2))
        p0, p1 = pRef*(1+0.1*rng.standard_normal(2))
        u0 = uRef + 30*rng.standard_normal(3)
        u1 = uRef + 30*rng.standard_normal(3)
        n = rng.standard_normal(3); n /= np.linalg.norm(n)
        a = keep_flux(ro0, u0, p0, ro1, u1, p1, n, pRef, np.float64)
        g = keep_flux_gauged(ro0, u0, p0, ro1, u1, p1, n, pRef, roRef, uRef, np.float64)
        r = ref_flux(roRef, uRef, pRef, n, np.float64)
        for i in range(3):
            d = np.abs(np.asarray(a[i]) - np.asarray(r[i]) - np.asarray(g[i]))
            s = np.abs(np.asarray(a[i])) + np.abs(np.asarray(r[i])) + 1.0
            err = max(err, float(np.max(d/s)))
    print(f"(1) analytic equivalence (float64): max rel err = {err:.3e}  ->",
          "PASS" if err < 1e-13 else "FAIL")

    # (2) 一様場 float32 ビット零
    f = F32
    ro = f(roRef); p = f(pRef); u = uRef.astype(f)
    worst = 0.0
    for _ in range(500):
        n = rng.standard_normal(3).astype(f); n /= f(np.linalg.norm(n))
        g = keep_flux_gauged(ro, u, p, ro, u, p, n, f(pRef), f(roRef), u, f)
        worst = max(worst, abs(float(g[0])), float(np.max(np.abs(g[1]))), abs(float(g[2])))
    print(f"(2) uniform-state float32 gauged flux: max |F| = {worst:.3e}  ->",
          "PASS" if worst == 0.0 else "FAIL")

    # (3) wavy メッシュのセル残差: 一様場ゲージ=厳密零 / 摂動場でゲージ有無が丸め一致
    if len(sys.argv) > 1:
        import h5py
        fh = h5py.File(sys.argv[1], "r")
        strct2 = fh["CELLS/STRUCT"][:]
        sv = fh["PLANES/surfVect"][:].astype(F32).reshape(-1, 3)
        nC = fh["CELLS/volume"].shape[0]
        cells = []
        icc = 0
        for _ic in range(nC):
            nn = strct2[icc]; icc += 1 + nn
            npl = strct2[icc]; icc += 1
            pl = strct2[icc:icc+npl]; icc += npl
            nd = strct2[icc]; icc += 1
            dr = strct2[icc:icc+nd]; icc += nd
            icc += 1
            cells.append((pl, dr))
        # 摂動場 (セル単位の擬似乱数; プレーン両側で同一状態を引けるようセル ID から生成)
        rng2 = np.random.default_rng(11)
        roc = (roRef*(1+1e-3*rng2.standard_normal(nC))).astype(F32)
        pc = (pRef*(1+1e-3*rng2.standard_normal(nC))).astype(F32)
        uc = (uRef[None, :] + 1.0*rng2.standard_normal((nC, 3))).astype(F32)
        # plane -> (ic0,ic1): 内部面のみ検証 (STRUCT の cell 側から復元)
        pc_map = {}
        for ic, (pl, dr) in enumerate(cells):
            for pidx, d in zip(pl, dr):
                pc_map.setdefault(pidx, []).append((ic, d))
        resU = np.zeros((nC, 5))  # 一様場ゲージ残差 (float32 で厳密 0 のはず)
        resA = np.zeros((nC, 5)); resG = np.zeros((nC, 5))
        uu = uRef.astype(F32); rr = F32(roRef); pp = F32(pRef)
        for pidx, adj in pc_map.items():
            if len(adj) != 2:
                continue
            (i0, d0), (i1, _d1) = adj if adj[0][1] > 0 else (adj[1], adj[0])
            s = sv[pidx]; S = F32(np.linalg.norm(s.astype(np.float64)))
            n = (s/S).astype(F32)
            gU = keep_flux_gauged(rr, uu, pp, rr, uu, pp, n, pp, rr, uu, F32)
            a = keep_flux(roc[i0], uc[i0], pc[i0], roc[i1], uc[i1], pc[i1], n, pp, F32)
            g = keep_flux_gauged(roc[i0], uc[i0], pc[i0], roc[i1], uc[i1], pc[i1], n, pp, rr, uu, F32)
            rf = ref_flux(rr, uu, pp, n, F32)
            # CELLS/STRUCT の dir=+1 は surfVect が「そのセルへ向く」規約 (検証済):
            # 流出 = -dir*F·s なので i0 (dir=+1) へ +vec, i1 へ -vec を足すと
            # ゲージ有無差が +F_ref(r_c) に一致する。
            for (dst, src) in ((resU, gU), (resA, a), (resG, g)):
                vec = np.array([src[0], src[1][0], src[1][1], src[1][2], src[2]], dtype=np.float64)*float(S)
                dst[i0] += vec; dst[i1] -= vec
            del rf
        rmsU = np.sqrt((resU**2).mean(axis=0))
        print(f"(3a) wavy uniform-state gauged cell residual (internal faces): max rms = {rmsU.max():.3e}  ->",
              "PASS" if rmsU.max() == 0.0 else "FAIL")
        # ゲージ有無差: (A - G) は Σ F_ref(s) = F_ref(r_c) に一致するはず (丸め内)
        diff = resA - resG
        scale = np.abs(resA).max(axis=0) + 1e-30
        # F_ref(r_c) 予測
        rcs = np.zeros((nC, 3), np.float64)
        for pidx, adj in pc_map.items():
            if len(adj) != 2:
                continue
            for ic, d in adj:
                rcs[ic] += d*sv[pidx].astype(np.float64)
        UnR = rcs @ uRef
        pred = np.stack([roRef*UnR,
                         roRef*UnR*uRef[0],
                         roRef*UnR*uRef[1],
                         roRef*UnR*uRef[2],
                         roRef*UnR*(0.5*uRef@uRef + pRef/roRef/(GA-1.0)) + pRef*UnR], axis=1)
        # pred≠0 の成分 (uRef 方向) は厳密一致、pred=0 の成分は面流束の丸めノイズのみ許容
        relerr = np.abs(diff-pred).max(axis=0)/(np.abs(pred).max(axis=0)+1e-30)
        active = np.abs(pred).max(axis=0) > 0
        noise = np.abs(diff).max(axis=0)
        ok = (relerr[active] < 1e-6).all() and (noise[~active] < 1e-3).all()
        print(f"(3b) (orig - gauged) == F_ref(r_c): rel err (active) = {relerr[active].max():.3e}, "
              f"noise (inactive) = {noise[~active].max():.3e}  ->", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
