"""⑤ SERN の平面 2D 構造 quad メッシュ (2 バンド + スリットカウル、決定的)。

plan: plans/active/tooling-nozzle-sern-chain.md §4.7 / S2。無次元 (H = 1) → `scale` [m]。

配置 (x 右、y 上):
    上境界 = ランプ (x<0 は直管 y=1、0..L_ramp は設計輪郭、以降は角 top_ext_angle の直線 = outlet 扱い)
    中間線 = カウル (x ≤ L_cowl: x<0 は y=0 の直管、0..L_cowl は −x tanθ_c0) → TE 以降は角 interface_angle の
             直線 (せん断層に沿わせる格子線。境界条件なし)
    下境界 = 遠方 (y = y_mid(x_out) − bot_depth)
バンド: 下 (下境界→中間線, nj_bot 点) / 上 (中間線→ランプ, nj_top 点)。x station は共通。
カウル (x ≤ L_cowl) は**スリット**: 中間線ノードを上下 2 重にし、上側 = cowl_in、下側 = cowl_out。
TE 以降の中間線ノードは共有 (内部線)。
物理タグ: PHYS_SERN。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .mesh2d import _radial_fracs

PHYS_SERN = {"inlet_nozzle": 1, "inlet_ext": 2, "outlet": 3, "ramp": 4, "cowl_in": 5,
             "cowl_out": 6, "bottom": 7, "top_out": 8, "fluid": 9}


@dataclass
class SernMeshParams:
    ni_up: int = 16
    ni_noz: int = 120
    ni_plume: int = 220
    nj_top: int = 101
    nj_bot: int = 61
    L_up: float = 0.5
    x_out_extra: float = 2.0
    bot_depth: float = 3.0
    first_wall_frac: float = 2.0e-3
    first_bot_frac: float = 0.02
    interface_angle: float = 0.0
    top_ext_angle: float = 0.0
    scale: float = 1.0
    x_cluster_w: float = 0.15
    x_cluster_a: float = 3.0


def _cluster_stations(x0, x1, n, ends=(True, True), w=0.15, a=3.0):
    xi = np.linspace(0.0, 1.0, 4001)
    dens = np.ones_like(xi)
    L = x1 - x0
    wn = max(w / L, 1e-3)
    if ends[0]:
        dens += (a - 1.0) * np.exp(-(xi / wn) ** 2)
    if ends[1]:
        dens += (a - 1.0) * np.exp(-((1.0 - xi) / wn) ** 2)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xi))])
    cum /= cum[-1]
    return x0 + L * np.interp(np.linspace(0.0, 1.0, n), cum, xi)


def _tanh_two_sided(n, first):
    xi = np.linspace(0.0, 1.0, n)

    def s_of(d):
        return 0.5 * (1.0 + np.tanh(d * (xi - 0.5)) / np.tanh(0.5 * d))
    lo, hi = 1e-3, 40.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if s_of(mid)[1] > first:
            lo = mid
        else:
            hi = mid
    s = s_of(0.5 * (lo + hi))
    s[0], s[-1] = 0.0, 1.0
    return s


def generate_sern_mesh(design, prm: SernMeshParams):
    """design: SernDesign (無次元)。戻り値 (coords [m], quads, bedges{name: edges}, info)。"""
    L_cowl = float(design.cowl_xy[-1, 0])
    y_te = float(design.cowl_xy[-1, 1])
    tan_c = -y_te / L_cowl if L_cowl > 0 else 0.0
    L_ramp = float(design.L_ramp)
    y_e = float(design.ramp_xy[-1, 1])
    x_out = L_ramp + prm.x_out_extra
    xs = np.concatenate([
        _cluster_stations(-prm.L_up, 0.0, prm.ni_up, (False, True), prm.x_cluster_w, prm.x_cluster_a),
        _cluster_stations(0.0, L_cowl, prm.ni_noz, (True, True), prm.x_cluster_w, prm.x_cluster_a)[1:],
        _cluster_stations(L_cowl, x_out, prm.ni_plume, (True, False), prm.x_cluster_w, prm.x_cluster_a)[1:],
    ])
    # ランプ後縁に station を置く (最寄りを置換)
    k = int(np.argmin(np.abs(xs - L_ramp)))
    xs[k] = L_ramp
    i_te = int(np.argmin(np.abs(xs - L_cowl)))
    assert abs(xs[i_te] - L_cowl) < 1e-12
    rx, ry = design.ramp_xy[:, 0], design.ramp_xy[:, 1]

    def y_top(x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0.0, 1.0,
                        np.where(x <= L_ramp, np.interp(x, rx, ry), y_e + (x - L_ramp) * np.tan(prm.top_ext_angle)))

    def y_mid(x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0.0, 0.0,
                        np.where(x <= L_cowl, -x * tan_c, y_te + (x - L_cowl) * np.tan(prm.interface_angle)))

    y_bot = float(y_mid(x_out)) - prm.bot_depth
    yt, ym = y_top(xs), y_mid(xs)
    s_top = _tanh_two_sided(prm.nj_top, prm.first_wall_frac)
    s_bot = _radial_fracs(prm.nj_bot, prm.first_wall_frac)   # 上端 (カウル外面/中間線) にクラスタ
    # 下端 (遠方) は粗いまま
    ni, njt, njb = len(xs), prm.nj_top, prm.nj_bot
    N_low = ni * njb
    N_up = ni * (njt - 1)
    # スリットの 2 重ノードは i < i_te のみ。TE ノード (i = i_te) は上下で共有する
    # (厚さ 0 のカウル後縁は 1 点。i_te も重複させると TE 直後のセル辺が上下で不一致になり
    #  幅 0 の隙間 = 未タグ境界辺が 2 本できる — 2026-09-04 メッシュテストで検出)
    N_dup = i_te
    coords = np.zeros((N_low + N_up + N_dup, 3))

    def low(i, j):
        return i * njb + j

    def up(i, j):
        if j == 0:
            return N_low + N_up + i if i < i_te else low(i, njb - 1)
        return N_low + i * (njt - 1) + (j - 1)

    for i in range(ni):
        yb = y_bot + s_bot * (ym[i] - y_bot)
        coords[low(i, 0):low(i, njb), 0] = xs[i]
        coords[low(i, 0):low(i, njb), 1] = yb
        yu = ym[i] + s_top * (yt[i] - ym[i])
        coords[up(i, 1):up(i, njt - 1) + 1, 0] = xs[i]
        coords[up(i, 1):up(i, njt - 1) + 1, 1] = yu[1:]
        if i < i_te:
            coords[up(i, 0)] = (xs[i], ym[i], 0.0)
    quads = []
    for i in range(ni - 1):
        for j in range(njb - 1):
            quads.append((low(i, j), low(i + 1, j), low(i + 1, j + 1), low(i, j + 1)))
        for j in range(njt - 1):
            quads.append((up(i, j), up(i + 1, j), up(i + 1, j + 1), up(i, j + 1)))
    quads = np.asarray(quads, dtype=np.int64)
    bedges = {n: [] for n in PHYS_SERN if n != "fluid"}
    for i in range(ni - 1):
        bedges["bottom"].append((low(i, 0), low(i + 1, 0)))
        xm = 0.5 * (xs[i] + xs[i + 1])
        (bedges["ramp"] if xm <= L_ramp else bedges["top_out"]).append((up(i, njt - 1), up(i + 1, njt - 1)))
        if i + 1 <= i_te:
            bedges["cowl_in"].append((up(i, 0), up(i + 1, 0)))
            bedges["cowl_out"].append((low(i, njb - 1), low(i + 1, njb - 1)))
    for j in range(njb - 1):
        bedges["inlet_ext"].append((low(0, j), low(0, j + 1)))
        bedges["outlet"].append((low(ni - 1, j), low(ni - 1, j + 1)))
    for j in range(njt - 1):
        bedges["inlet_nozzle"].append((up(0, j), up(0, j + 1)))
        bedges["outlet"].append((up(ni - 1, j), up(ni - 1, j + 1)))
    coords *= prm.scale
    info = {"ni": ni, "nj_top": njt, "nj_bot": njb, "cells": int(quads.shape[0]), "nodes": int(coords.shape[0]),
            "x_out": x_out, "y_bot": y_bot, "i_te": i_te, "L_cowl": L_cowl, "L_ramp": L_ramp}
    return coords, quads, bedges, info, y_mid


def write_msh41_named(path, coords, quads, bedges: dict, phys: dict) -> None:
    """任意の境界名セットで gmsh msh4.1 を書く (mesh2d.write_msh41_2d の一般化)。"""
    names = [n for n in phys if n != "fluid" and len(bedges.get(n, [])) > 0]
    n_nodes = coords.shape[0]
    n_elems = quads.shape[0] + sum(len(bedges[g]) for g in names)
    mn, mx = coords.min(0), coords.max(0)
    L = []
    ap = L.append
    ap("$MeshFormat\n4.1 0 8\n$EndMeshFormat")
    ap(f"$PhysicalNames\n{len(names) + 1}")
    for nm in names:
        ap(f'1 {phys[nm]} "{nm}"')
    ap(f'2 {phys["fluid"]} "fluid"')
    ap("$EndPhysicalNames")
    ap("$Entities")
    ap(f"0 {len(names)} 1 0")
    bb = f"{mn[0]:.9g} {mn[1]:.9g} {mn[2]:.9g} {mx[0]:.9g} {mx[1]:.9g} {mx[2]:.9g}"
    for ci, nm in enumerate(names, start=1):
        ap(f"{ci} {bb} 1 {phys[nm]} 0")
    ap(f"1 {bb} 1 {phys['fluid']} 0")
    ap("$EndEntities")
    ap("$Nodes")
    ap(f"1 {n_nodes} 1 {n_nodes}")
    ap(f"2 1 0 {n_nodes}")
    ap("\n".join(str(i + 1) for i in range(n_nodes)))
    ap("\n".join(f"{c[0]:.10g} {c[1]:.10g} {c[2]:.10g}" for c in coords))
    ap("$EndNodes")
    ap("$Elements")
    ap(f"{1 + len(names)} {n_elems} 1 {n_elems}")
    et = 1
    ap(f"2 1 3 {quads.shape[0]}")
    ap("\n".join(f"{et + k} {q[0]+1} {q[1]+1} {q[2]+1} {q[3]+1}" for k, q in enumerate(quads)))
    et += quads.shape[0]
    for ci, nm in enumerate(names, start=1):
        edges = bedges[nm]
        ap(f"1 {ci} 1 {len(edges)}")
        ap("\n".join(f"{et + k} {a+1} {b+1}" for k, (a, b) in enumerate(edges)))
        et += len(edges)
    ap("$EndElements")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
