#!/usr/bin/env python3
"""forge 結果から推力・境界別力・断面運動量流束を積分する後処理ツール。

使い方 (run ディレクトリで):
    # 壁積分: 境界面ダンプ res_<name>_<physID>_<step>.h5 から力 (圧力+摩擦)
    python3 integrate_forces.py wall <run_dir> [--step 最終] [--pa 1e5]

    # 断面積分: 体積 res_<step>.h5 + forge.h5 の節点座標から x=const 断面の ṁ / 運動量推力
    python3 integrate_forces.py slice <run_dir> --x 0.035 [--rmax 0.00433] [--step N] [--pa 1e5]

    # 時系列: 全 res_*.h5 で wall/slice を回し頭打ち確認 (派生量の準定常ルール)
    python3 integrate_forces.py history <run_dir> [--x 0.035 --rmax 0.00433] [--pa 1e5]

前提・定義:
- 半割モデル (y>=0)。y 成分は対称で相殺するため |Fy| を検証表示し、x/z は 2 倍して報告。
- 力の定義: 構造に働く力 F = ∫(Ps-pa)(+n_out)dA + ∫τw dA。n_out は流体から外向き
  (=構造へ向く)。面ダンプの三角形 winding から法線を作り、既知境界 (inlet の -z) で
  向きを較正して全境界に適用する。
- slice: 節点場を極座標グリッド (r,θ∈[0,π]) に逆距離補間 (KDTree k=8) して積分。
  ṁ=∫ρux dA, F=∫(ρux²+Ps-pa)dA。ParaView 不要の純 Python。
- 摩擦 τw は面ダンプの twall_x/y/z (節点値)。node の壁は nodeWallDirichlet=1 で u=0 pinned。
- Isp = F/ṁ [m/s]、効率 η=Isp/Isp_theo (--isp-theo, 既定 2212 = Pan/Chen/Ye 2017)。
"""
import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from res_h5_to_vtu import parse_conne, _count_cells


# ---------------- 共通 ----------------

def find_steps(run, pattern):
    out = {}
    for p in glob.glob(os.path.join(run, pattern)):
        m = re.search(r"_(\d+)\.h5$", p)
        if m:
            out.setdefault(int(m.group(1)), []).append(p)
    return out


def surface_geometry(f):
    """面ダンプの三角形/四角形ごとの面積ベクトル (winding 由来) と節点->面平均行列を返す。"""
    coord = np.array(f["MESH/COORD"]).reshape(-1, 3)
    conne = np.array(f["MESH/CONNE"])
    ncell = _count_cells(conne)
    conn, offs, types = parse_conne(conne, ncell)
    sv = np.zeros((ncell, 3))
    faces = []
    s = 0
    for i, o in enumerate(offs):
        idx = conn[s:o]
        s = o
        p = coord[idx]
        if len(idx) == 3:
            sv[i] = 0.5 * np.cross(p[1] - p[0], p[2] - p[0])
        else:  # quad: 対角線ベクトル積
            sv[i] = 0.5 * np.cross(p[2] - p[0], p[3] - p[1])
        faces.append(idx)
    return coord, faces, sv


def face_mean(vals, faces):
    return np.array([vals[list(ix)].mean() for ix in faces])


def wall_force(path, pa):
    """面ダンプ 1 つの力 (圧力/摩擦、半割そのまま)。法線は winding 由来 (未較正)。"""
    with h5py.File(path) as f:
        coord, faces, sv = surface_geometry(f)
        Ps = face_mean(np.array(f["VALUE/Ps"]), faces)
        if "VALUE/twall_x" in f:   # 摩擦は壁系ダンプのみ (inlet/outlet 等には無い)
            tw = np.stack([face_mean(np.array(f["VALUE/twall_" + c]), faces)
                           for c in ("x", "y", "z")], axis=1)
        else:
            tw = np.zeros((len(faces), 3))
    area = np.linalg.norm(sv, axis=1)
    Fp = ((Ps - pa)[:, None] * sv).sum(axis=0)   # 圧力 (n の向きは呼び出し側で較正)
    Ff = (tw * area[:, None]).sum(axis=0)        # 摩擦
    return Fp, Ff, area.sum(), sv.sum(axis=0)


def calibrate_sign(run, step, pa):
    """inlet ダンプの面法線合計は -z を向くべき (流体外向き=下向き)。符号を返す。"""
    cands = glob.glob(os.path.join(run, f"res_inlet_*_{step}.h5"))
    if not cands:
        return 1.0, "inlet ダンプ無し (符号未較正: winding のまま)"
    with h5py.File(cands[0]) as f:
        _, _, sv = surface_geometry(f)
    nz = sv.sum(axis=0)[2]
    return (1.0 if nz < 0 else -1.0), f"inlet Σn_z={nz:.3e} -> sign={'+1' if nz<0 else '-1'}"


def run_wall(run, step, pa, isp_theo):
    dumps = find_steps(run, "res_*_*.h5")
    surf_steps = sorted(s for s in dumps
                        if any(re.search(r"res_[A-Za-z_]+_\d+_%d\.h5$" % s, p) for p in dumps[s]))
    if step is None:
        step = max(surf_steps)
    sgn, note = calibrate_sign(run, step, pa)
    print(f"[wall] step={step} pa={pa:.3g} 法線較正: {note}")
    total_p = np.zeros(3)
    total_f = np.zeros(3)
    rows = []
    for p in sorted(glob.glob(os.path.join(run, f"res_*_{step}.h5"))):
        m = re.match(r"res_([A-Za-z_]+)_(\d+)_%d\.h5$" % step, os.path.basename(p))
        if not m:
            continue
        name = m.group(1)
        Fp, Ff, area, _ = wall_force(p, pa)
        Fp = sgn * Fp
        rows.append((name, Fp, Ff, area))
        if name.startswith("wall") or name == "base":
            total_p += Fp
            total_f += Ff
    print(f"{'境界':14s} {'A[m2]':>10s} {'Fp_x[N]':>12s} {'Ff_x[N]':>12s} {'Fp_y':>10s} {'Fp_z':>12s}")
    for name, Fp, Ff, area in rows:
        print(f"{name:14s} {area:10.4g} {2*Fp[0]:12.4f} {2*Ff[0]:12.4f} {Fp[1]:10.2e} {2*Fp[2]:12.4f}")
    Fx = 2 * (total_p[0] + total_f[0])   # 半割 -> 全周
    print(f"\n壁合計 (wall+wall_pintle+base, 全周換算): Fx = {Fx:.3f} N "
          f"(圧力 {2*total_p[0]:.3f} + 摩擦 {2*total_f[0]:.3f})")
    print(f"  参考 Fy(半割) = {total_p[1]+total_f[1]:.3e} N "
          f"(半割モデルでは対称面の圧力と釣り合う量; ゼロにはならない)")
    for name, Fp, Ff, _ in rows:
        if name == "wall_pintle":
            print(f"  ピントル軸力 (全周): {2*(Fp[0]+Ff[0]):.4f} N "
                  f"(圧力 {2*Fp[0]:.4f} + 摩擦 {2*Ff[0]:.4f})")
    return Fx, step


def mdot_from_dump(run, step, prefix="res_inlet"):
    cands = glob.glob(os.path.join(run, f"{prefix}_*_{step}.h5"))
    if not cands:
        return None
    with h5py.File(cands[0]) as f:
        coord, faces, sv = surface_geometry(f)
        ro = face_mean(np.array(f["VALUE/ro"]), faces)
        u = np.stack([face_mean(np.array(f["VALUE/roU" + c]), faces) for c in ("x", "y", "z")],
                     axis=1) / ro[:, None]
    # 流入質量流量 (向きは面法線と U の内積の符号で自動; 半割->2倍)
    md = (ro * np.abs((u * sv).sum(axis=1) / np.maximum(np.linalg.norm(sv, axis=1), 1e-30))
          * np.linalg.norm(sv, axis=1)).sum()
    return 2 * md


def run_slice(run, x0, rmax, step, pa, isp_theo, quiet=False):
    from scipy.spatial import cKDTree
    vols = find_steps(run, "res_*.h5")
    vsteps = sorted(s for s in vols
                    if any(re.fullmatch(r".*/res_\d+\.h5", p) for p in vols[s]))
    if step is None:
        step = max(vsteps)
    vres = [p for p in vols[step] if re.fullmatch(r".*/res_\d+\.h5", p)][0]
    mesh = os.path.join(run, "forge.h5")
    with h5py.File(mesh) as f:
        X = np.array(f["MESH/COORD"]).reshape(-1, 3)
    with h5py.File(vres) as f:
        ro = np.array(f["VALUE/ro"])
        Ux = np.array(f["VALUE/Ux"])
        P = np.array(f["VALUE/P"])
    # 断面近傍の節点だけで KDTree (高速化)
    band = np.abs(X[:, 0] - x0) < max(0.004, rmax * 0.5)
    tree = cKDTree(X[band])
    rob, Uxb, Pb = ro[band], Ux[band], P[band]
    # 極座標サンプリング (y>=0 半円)
    nr, nt = 160, 90
    rc = (np.arange(nr) + 0.5) / nr * rmax
    tc = (np.arange(nt) + 0.5) / nt * np.pi
    R, T = np.meshgrid(rc, tc, indexing="ij")
    pts = np.stack([np.full(R.size, x0), (R * np.sin(T)).ravel(), (R * np.cos(T)).ravel()], axis=1)
    d, idx = tree.query(pts, k=8)
    w = 1.0 / np.maximum(d, 1e-12) ** 2
    w /= w.sum(axis=1, keepdims=True)
    roS = (rob[idx] * w).sum(axis=1)
    UxS = (Uxb[idx] * w).sum(axis=1)
    PS = (Pb[idx] * w).sum(axis=1)
    dA = (R * (rmax / nr) * (np.pi / nt)).ravel()   # r dr dθ
    md = 2 * (roS * UxS * dA).sum()
    Fm = 2 * (roS * UxS ** 2 * dA).sum()
    Fpr = 2 * ((PS - pa) * dA).sum()
    F = Fm + Fpr
    if not quiet:
        print(f"[slice] x={x0*1000:.1f}mm rmax={rmax*1000:.2f}mm step={step} ({os.path.basename(vres)})")
        print(f"  ṁ = {md:.4f} kg/s | F_mom = {Fm:.3f} N | F_pres = {Fpr:.3f} N | F = {F:.3f} N")
        print(f"  Isp = F/ṁ = {F/md:.1f} m/s | η = {F/md/isp_theo*100:.2f}% (Isp_theo {isp_theo})")
    return dict(step=step, mdot=md, F=F, Fm=Fm, Fp=Fpr)


def run_history(run, x0, rmax, pa, isp_theo):
    vols = find_steps(run, "res_*.h5")
    vsteps = sorted(s for s in vols
                    if any(re.fullmatch(r".*/res_\d+\.h5", p) for p in vols[s]))
    print(f"[history] {len(vsteps)} snapshots: F(t) at x={x0*1000:.1f}mm")
    hist = []
    for s in vsteps:
        r = run_slice(run, x0, rmax, s, pa, isp_theo, quiet=True)
        hist.append(r)
        print(f"  step {s:6d}: ṁ={r['mdot']:.4f} kg/s  F={r['F']:.3f} N  Isp={r['F']/r['mdot']:.1f} m/s")
    if len(hist) >= 3:
        F = np.array([h["F"] for h in hist])
        tail = F[-3:]
        drift = abs(tail[-1] - tail[0]) / max(abs(tail).max(), 1e-30)
        verdict = "SETTLED" if drift < 0.01 else ("DRIFTING" if drift < 0.05 else "TRANSIENT")
        print(f"  末尾3点の相対ドリフト {drift*100:.2f}% -> {verdict}")
    return hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["wall", "slice", "history"])
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--pa", type=float, default=1.0e5)
    ap.add_argument("--x", type=float, default=0.035)
    ap.add_argument("--rmax", type=float, default=0.00433)
    ap.add_argument("--isp-theo", type=float, default=2212.0)
    a = ap.parse_args()
    if a.mode == "wall":
        run_wall(a.run, a.step, a.pa, a.isp_theo)
        md = mdot_from_dump(a.run, a.step if a.step else
                            max(find_steps(a.run, "res_inlet_*.h5") or [0]))
        if md:
            print(f"  ṁ (inlet ダンプ) = {md:.4f} kg/s")
    elif a.mode == "slice":
        run_slice(a.run, a.x, a.rmax, a.step, a.pa, a.isp_theo)
    else:
        run_history(a.run, a.x, a.rmax, a.pa, a.isp_theo)


if __name__ == "__main__":
    main()
