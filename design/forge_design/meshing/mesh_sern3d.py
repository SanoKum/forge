"""⑤ SERN の 3D 構造 hex メッシュ (S7: 側壁・有限スパンの確認用)。

2D の 2 バンド構造格子 (mesh_sern) を z 方向に押し出す。z=0 は対称面 (半スパン)、z = W/2 に側壁 (スリット面)、
その外側 z ∈ (W/2, Z_far] は外部流が回り込む空間。ランプ (上境界) は全スパン (機体下面)、カウルは z ≤ W/2 の
有限幅の板 (z > W/2 では中間線ノードを共有 = 板なし)。側壁は x ≤ L_sw の上バンド (カウル〜ランプ) にある。

ノード重複 (スリット):
  D1 カウル: i < i_te かつ z_k ≤ W/2 の中間線ノードを上下 2 重 (2D と同じ)。
  D2 側壁: k = k_sw (z = W/2), i < i_sw (後縁 station は共有), 上バンド内部 j (jm < j < NJ−1) を内外 2 重。j = jm では内側 = cowl_in 側の
     上コピー、外側 = 中間線の元ノード (cowl_out 兼)、j = NJ−1 (ランプ) は共有。
境界面 (quad): inlet_nozzle / inlet_ext / outlet / ramp / top_out / cowl_in / cowl_out / bottom / sym (z=0) /
              side_far (z = Z_far) / sidewall_in / sidewall_out。
hex は gmsh 型 5 (底面 4 点 CCW → 上面 4 点)。座標一致ノードを持つので stage 間 restart は index コピーにすること。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .mesh2d import _radial_fracs
from .mesh_sern import _cluster_stations, _tanh_two_sided

PHYS_SERN3D = {"inlet_nozzle": 1, "inlet_ext": 2, "outlet": 3, "ramp": 4, "cowl_in": 5, "cowl_out": 6, "bottom": 7,
               "top_out": 8, "sym": 9, "side_far": 10, "sidewall_in": 11, "sidewall_out": 12, "fluid": 13}


@dataclass
class SernMesh3DParams:
    ni_up: int = 10
    ni_noz: int = 60
    ni_plume: int = 110
    nj_top: int = 49
    nj_bot: int = 31
    nz_in: int = 25          # z ∈ [0, W/2]
    nz_out: int = 17         # z ∈ (W/2, Z_far]
    W: float = 2.0           # ノズル幅 / H (全幅)
    Z_ext: float = 1.5       # 側壁外側の空間 / H
    L_sw: float | None = None  # 側壁の x 範囲 (None → L_cowl)
    L_up: float = 0.5
    x_out_extra: float = 2.0
    bot_depth: float = 3.0
    first_wall_frac: float = 4.0e-3
    first_z_frac: float = 4.0e-3
    interface_angle: float = 0.0
    top_ext_angle: float = 0.0
    scale: float = 1.0
    x_cluster_w: float = 0.15
    x_cluster_a: float = 3.0


def generate_sern_mesh3d(design, prm: SernMesh3DParams):
    L_cowl = float(design.cowl_xy[-1, 0]); y_te = float(design.cowl_xy[-1, 1])
    tan_c = -y_te / L_cowl if L_cowl > 0 else 0.0
    L_ramp = float(design.L_ramp); y_e = float(design.ramp_xy[-1, 1])
    L_sw = L_cowl if prm.L_sw is None else float(prm.L_sw)
    x_out = L_ramp + prm.x_out_extra
    xs = np.concatenate([
        _cluster_stations(-prm.L_up, 0.0, prm.ni_up, (False, True), prm.x_cluster_w, prm.x_cluster_a),
        _cluster_stations(0.0, L_cowl, prm.ni_noz, (True, True), prm.x_cluster_w, prm.x_cluster_a)[1:],
        _cluster_stations(L_cowl, x_out, prm.ni_plume, (True, False), prm.x_cluster_w, prm.x_cluster_a)[1:],
    ])
    xs[int(np.argmin(np.abs(xs - L_ramp)))] = L_ramp
    i_te = int(np.argmin(np.abs(xs - L_cowl))); assert abs(xs[i_te] - L_cowl) < 1e-12
    i_sw = int(np.argmin(np.abs(xs - L_sw)))
    rx, ry = design.ramp_xy[:, 0], design.ramp_xy[:, 1]

    def y_top(x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0.0, 1.0, np.where(x <= L_ramp, np.interp(x, rx, ry), y_e + (x - L_ramp) * np.tan(prm.top_ext_angle)))

    def y_mid(x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0.0, 0.0, np.where(x <= L_cowl, -x * tan_c, y_te + (x - L_cowl) * np.tan(prm.interface_angle)))
    y_bot = float(y_mid(x_out)) - prm.bot_depth
    yt, ym = y_top(xs), y_mid(xs)
    ni, njt, njb = len(xs), prm.nj_top, prm.nj_bot
    NJ = njb + njt - 1; jm = njb - 1
    # z 分布: [0, W/2] は側壁側 (z=W/2) にクラスタ、(W/2, Z_far] は側壁側にクラスタ
    hw = 0.5 * prm.W
    s_in = _radial_fracs(prm.nz_in, min(prm.first_z_frac / hw, 0.5 / (prm.nz_in - 1)))
    z_in = hw * s_in
    s_out = _radial_fracs(prm.nz_out, min(prm.first_z_frac / prm.Z_ext, 0.5 / (prm.nz_out - 1)))
    z_out = hw + prm.Z_ext * (1.0 - s_out[::-1])   # 側壁側が細かい
    zs = np.concatenate([z_in, z_out[1:]]); nz = len(zs); k_sw = prm.nz_in - 1
    # --- 2D 断面 (各 station の y 列) ---
    Y2 = np.zeros((ni, NJ))
    for i in range(ni):
        h_lo = max(ym[i] - y_bot, 1e-12); h_up = max(yt[i] - ym[i], 1e-12)
        s_bot = _radial_fracs(njb, min(prm.first_wall_frac / h_lo, 0.5 / (njb - 1)))
        s_top = _tanh_two_sided(njt, min(prm.first_wall_frac / h_up, 0.5 / (njt - 1)))
        Y2[i, :njb] = y_bot + s_bot * (ym[i] - y_bot)
        Y2[i, jm:] = ym[i] + s_top * (yt[i] - ym[i])
    # --- ノード番号 ---
    N_base = ni * NJ * nz
    def base(i, j, k): return (i * NJ + j) * nz + k
    dup1 = {}   # (i,k) -> id  (カウル上コピー)
    dup2 = {}   # (i,j) -> id  (側壁外コピー, k=k_sw)
    nid = N_base
    for i in range(i_te):
        for k in range(k_sw + 1):
            dup1[(i, k)] = nid; nid += 1
    for i in range(i_sw):                      # 側壁後縁 (i_sw) は共有 (カウル TE と同じ)
        for j in range(jm + 1, NJ - 1):
            dup2[(i, j)] = nid; nid += 1
    coords = np.zeros((nid, 3))
    for i in range(ni):
        for j in range(NJ):
            b = base(i, j, 0)
            coords[b:b + nz, 0] = xs[i]; coords[b:b + nz, 1] = Y2[i, j]; coords[b:b + nz, 2] = zs
    for (i, k), n in dup1.items():
        coords[n] = (xs[i], ym[i], zs[k])
    for (i, j), n in dup2.items():
        coords[n] = (xs[i], Y2[i, j], zs[k_sw])

    def node(i, j, k, side: str):
        """side: 'lo' (下バンド), 'up_in' (上バンド, 内側 z<=W/2), 'up_out' (上バンド, 外側 z>W/2)。"""
        if j == jm:
            if side == "up_in" and (i, k) in dup1:
                return dup1[(i, k)]
            return base(i, j, k)
        if side == "up_out" and k == k_sw and (i, j) in dup2:
            return dup2[(i, j)]
        return base(i, j, k)

    hexes = []
    for i in range(ni - 1):
        for k in range(nz - 1):
            for j in range(NJ - 1):
                if j < jm:
                    sd = ("lo", "lo")
                elif k < k_sw:
                    sd = ("up_in", "up_in")
                elif k == k_sw:
                    sd = ("up_in", "up_out")    # k=k_sw 側は内側面 (境界 k) と外側 (k+1)
                else:
                    sd = ("up_out", "up_out")
                # k = k_sw の cell (k_sw → k_sw+1) は外側 cell: 面 k_sw のノードは outer コピー
                s0, s1 = ("up_out", "up_out") if (j >= jm and k == k_sw) else sd
                n = lambda ii, jj, kk, ss: node(ii, jj, kk, ss)
                hexes.append((n(i, j, k, s0), n(i + 1, j, k, s0), n(i + 1, j + 1, k, s0), n(i, j + 1, k, s0),
                              n(i, j, k + 1, s1), n(i + 1, j, k + 1, s1), n(i + 1, j + 1, k + 1, s1), n(i, j + 1, k + 1, s1)))
    hexes = np.asarray(hexes, dtype=np.int64)
    # --- 境界 quad ---
    B = {n: [] for n in PHYS_SERN3D if n != "fluid"}
    for k in range(nz - 1):
        side_k = "up_in" if k < k_sw else "up_out"
        side_k1 = "up_in" if k + 1 <= k_sw and k + 1 != k_sw else ("up_in" if k + 1 < k_sw else "up_out")
        for j in range(NJ - 1):
            sd = "lo" if j < jm else ("up_in" if k < k_sw else "up_out")
            sd1 = "lo" if j < jm else ("up_in" if k + 1 < k_sw or (k + 1 == k_sw and k < k_sw) else "up_out")
            # k+1 == k_sw の面ノードは内側 (inner) — 外側 cell は k = k_sw から
            nm = "inlet_ext" if j < jm else "inlet_nozzle"
            B[nm].append((node(0, j, k, sd), node(0, j + 1, k, sd), node(0, j + 1, k + 1, sd1), node(0, j, k + 1, sd1)))
            B["outlet"].append((node(ni - 1, j, k, sd), node(ni - 1, j, k + 1, sd1), node(ni - 1, j + 1, k + 1, sd1), node(ni - 1, j + 1, k, sd)))
    for i in range(ni - 1):
        for k in range(nz - 1):
            sdk = "up_in" if k < k_sw else "up_out"; sdk1 = "up_in" if k + 1 < k_sw or k + 1 == k_sw and k < k_sw else "up_out"
            B["bottom"].append((base(i, 0, k), base(i + 1, 0, k), base(i + 1, 0, k + 1), base(i, 0, k + 1)))
            xm = 0.5 * (xs[i] + xs[i + 1])
            top = (node(i, NJ - 1, k, sdk), node(i, NJ - 1, k + 1, sdk1), node(i + 1, NJ - 1, k + 1, sdk1), node(i + 1, NJ - 1, k, sdk))
            (B["ramp"] if xm <= L_ramp else B["top_out"]).append(top)
            if i + 1 <= i_te and k + 1 <= k_sw:
                B["cowl_in"].append((node(i, jm, k, "up_in"), node(i + 1, jm, k, "up_in"), node(i + 1, jm, k + 1, "up_in"), node(i, jm, k + 1, "up_in")))
                B["cowl_out"].append((base(i, jm, k), base(i, jm, k + 1), base(i + 1, jm, k + 1), base(i + 1, jm, k)))
        for j in range(NJ - 1):
            sd = "lo" if j < jm else "up_in"
            B["sym"].append((node(i, j, 0, sd), node(i, j + 1, 0, sd), node(i + 1, j + 1, 0, sd), node(i + 1, j, 0, sd)))
            sdo = "lo" if j < jm else "up_out"
            B["side_far"].append((node(i, j, nz - 1, sdo), node(i + 1, j, nz - 1, sdo), node(i + 1, j + 1, nz - 1, sdo), node(i, j + 1, nz - 1, sdo)))
            if i + 1 <= i_sw and j >= jm:
                B["sidewall_in"].append((node(i, j, k_sw, "up_in"), node(i + 1, j, k_sw, "up_in"), node(i + 1, j + 1, k_sw, "up_in"), node(i, j + 1, k_sw, "up_in")))
                B["sidewall_out"].append((node(i, j, k_sw, "up_out"), node(i, j + 1, k_sw, "up_out"), node(i + 1, j + 1, k_sw, "up_out"), node(i + 1, j, k_sw, "up_out")))
    coords *= prm.scale
    info = {"ni": ni, "NJ": NJ, "nz": nz, "jm": jm, "k_sw": k_sw, "i_te": i_te, "i_sw": i_sw, "cells": int(hexes.shape[0]),
            "nodes": int(coords.shape[0]), "W": prm.W, "Z_far": float(zs[-1]), "L_sw": L_sw, "x_out": x_out, "y_bot": y_bot,
            "L_cowl": L_cowl, "L_ramp": L_ramp, "n_dup_cowl": len(dup1), "n_dup_side": len(dup2)}
    return coords, hexes, B, info, y_mid


def write_msh41_3d(path, coords, hexes, bquads: dict, phys: dict) -> None:
    names = [n for n in phys if n != "fluid" and len(bquads.get(n, [])) > 0]
    n_nodes = coords.shape[0]
    n_elems = hexes.shape[0] + sum(len(bquads[g]) for g in names)
    mn, mx = coords.min(0), coords.max(0)
    L = []; ap = L.append
    ap("$MeshFormat\n4.1 0 8\n$EndMeshFormat")
    ap(f"$PhysicalNames\n{len(names) + 1}")
    for nm in names:
        ap(f'2 {phys[nm]} "{nm}"')
    ap(f'3 {phys["fluid"]} "fluid"')
    ap("$EndPhysicalNames")
    ap("$Entities"); ap(f"0 0 {len(names)} 1")
    bb = f"{mn[0]:.9g} {mn[1]:.9g} {mn[2]:.9g} {mx[0]:.9g} {mx[1]:.9g} {mx[2]:.9g}"
    for si, nm in enumerate(names, start=1):
        ap(f"{si} {bb} 1 {phys[nm]} 0")
    ap(f"1 {bb} 1 {phys['fluid']} 0")
    ap("$EndEntities")
    ap("$Nodes"); ap(f"1 {n_nodes} 1 {n_nodes}"); ap(f"3 1 0 {n_nodes}")
    ap("\n".join(str(i + 1) for i in range(n_nodes)))
    ap("\n".join(f"{c[0]:.10g} {c[1]:.10g} {c[2]:.10g}" for c in coords))
    ap("$EndNodes")
    ap("$Elements"); ap(f"{1 + len(names)} {n_elems} 1 {n_elems}")
    et = 1
    ap(f"3 1 5 {hexes.shape[0]}")
    ap("\n".join(f"{et + k} " + " ".join(str(v + 1) for v in h) for k, h in enumerate(hexes)))
    et += hexes.shape[0]
    for si, nm in enumerate(names, start=1):
        q = bquads[nm]
        ap(f"2 {si} 3 {len(q)}")
        ap("\n".join(f"{et + k} {a+1} {b+1} {c+1} {d+1}" for k, (a, b, c, d) in enumerate(q)))
        et += len(q)
    ap("$EndElements")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
