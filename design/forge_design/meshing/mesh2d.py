"""軸対称ノズルの構造化 2D quad メッシュ (トポロジ固定・決定的)。

- x 方向: スロート近傍を滑らかに細分する間隔関数の逆積分で station を配置。
- r 方向: 全 station 共通の正規化分布 s_j (壁側幾何級数クラスタリング)。
  r_j(i) = r_w(x_i) * s_j — 同一トポロジで壁だけ動かせる (親計画 §4.3 R1)。
- 出力: gmsh msh4.1 テキスト (物理タグ inlet=1 outlet=2 wall=3 axis=4 fluid=5)
  → 既存 convertGmshToForge で forge h5 化 (wall_dist・幾何量は変換器が計算)。

軸対称の座標規約: x=軸, y=半径, z=0 平面, 上半分のみ (methods/axisymmetric)。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PHYS = {"inlet": 1, "outlet": 2, "wall": 3, "axis": 4, "fluid": 5}


@dataclass
class Mesh2DParams:
    ni: int = 241              # x 方向ノード数
    nj: int = 81               # r 方向ノード数
    wall_first_frac: float = 2.0e-3  # 壁第一セル厚 / 局所半径
    throat_refine: float = 3.0       # スロート近傍の間隔比 (far/throat)
    throat_width: float = 1.5        # 細分の e^-x^2 幅 (r* 単位)
    scale: float = 1.0               # r* [m] (無次元 → m)
    # --- 接合部近傍の局所細分 (B10: x_reach のアンカー安定性・こぶの格子収束用) ---
    local_center: float = 0.0        # 細分中心 [r*] (0 で無効)
    local_refine: float = 1.0        # 中心での間隔比 (1.0 で無効)
    local_width: float = 0.75        # 細分の e^-x^2 幅 [r*]


def _x_stations(x0: float, x1: float, ni: int, refine: float, width: float,
                local_center: float = 0.0, local_refine: float = 1.0,
                local_width: float = 0.75) -> np.ndarray:
    """間隔 h(x) ∝ 1/(密度) の逆積分で station を置く。

    密度 = スロート細分 (中心 x=0) × 局所細分 (中心 `local_center`)。後者は
    **接合部近傍だけ**を細かくして、$x_{reach}$ の局所微分アンカーと接合こぶの
    格子収束を、全域解像度を上げずに調べるためのもの (B10)。
    """
    xs = np.linspace(x0, x1, 8001)
    dens = 1.0 + (refine - 1.0) * np.exp(-((xs / width) ** 2))
    if local_refine > 1.0:
        dens = dens * (1.0 + (local_refine - 1.0)
                       * np.exp(-(((xs - local_center) / local_width) ** 2)))
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xs))])
    cum /= cum[-1]
    return np.interp(np.linspace(0.0, 1.0, ni), cum, xs)


def _radial_fracs(nj: int, first_frac: float) -> np.ndarray:
    """s_0=0 (軸) → s_{nj-1}=1 (壁)。壁側の最終区間 = first_frac。幾何級数。"""
    n = nj - 1
    if first_frac * n >= 1.0:
        return np.linspace(0.0, 1.0, nj)
    # 比 q>1 を解く: first*(q^n - 1)/(q - 1) = 1
    q = 1.5
    for _ in range(100):
        f = first_frac * (q ** n - 1.0) / (q - 1.0) - 1.0
        df = first_frac * ((n * q ** (n - 1)) * (q - 1.0) - (q ** n - 1.0)) / (q - 1.0) ** 2
        qn = q - f / df
        if abs(qn - q) < 1e-14:
            q = qn
            break
        q = qn
    gaps = first_frac * q ** np.arange(n)          # 壁側から軸側へ拡大
    s = np.concatenate([[0.0], np.cumsum(gaps[::-1])])
    return s / s[-1]


def generate_axisym_mesh(wall, prm: Mesh2DParams):
    """wall: NozzleWall。戻り値 (coords (N,3) [m], quads (M,4), 境界辺 dict)。"""
    xs = _x_stations(wall.x_in, wall.x_e, prm.ni, prm.throat_refine, prm.throat_width,
                     prm.local_center, prm.local_refine, prm.local_width)
    rw = wall.r(xs)
    s = _radial_fracs(prm.nj, prm.wall_first_frac)
    ni, nj = prm.ni, prm.nj
    X = np.repeat(xs[:, None], nj, axis=1)
    R = rw[:, None] * s[None, :]
    coords = np.zeros((ni * nj, 3))
    coords[:, 0] = X.ravel() * prm.scale
    coords[:, 1] = R.ravel() * prm.scale

    def nid(i, j):
        return i * nj + j

    quads = np.empty(((ni - 1) * (nj - 1), 4), dtype=np.int64)
    k = 0
    for i in range(ni - 1):
        for j in range(nj - 1):
            quads[k] = (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
            k += 1
    bedges = {
        "axis": [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)],
        "wall": [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)],
        "inlet": [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)],
        "outlet": [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)],
    }
    return coords, quads, bedges


def write_msh41_2d(path, coords, quads, bedges) -> None:
    """2D 平面メッシュを gmsh msh4.1 テキストで書く (決定的)。"""
    curve_order = ["inlet", "outlet", "wall", "axis"]
    n_nodes = coords.shape[0]
    n_elems = quads.shape[0] + sum(len(bedges[g]) for g in curve_order)
    mn = coords.min(0)
    mx = coords.max(0)
    L = []
    ap = L.append
    ap("$MeshFormat\n4.1 0 8\n$EndMeshFormat")
    ap("$PhysicalNames\n5")
    for nm in curve_order:
        ap(f'1 {PHYS[nm]} "{nm}"')
    ap(f'2 {PHYS["fluid"]} "fluid"')
    ap("$EndPhysicalNames")
    ap("$Entities")
    ap(f"0 {len(curve_order)} 1 0")
    bb = f"{mn[0]:.9g} {mn[1]:.9g} {mn[2]:.9g} {mx[0]:.9g} {mx[1]:.9g} {mx[2]:.9g}"
    for ci, nm in enumerate(curve_order, start=1):
        ap(f"{ci} {bb} 1 {PHYS[nm]} 0")
    ap(f"1 {bb} 1 {PHYS['fluid']} 0")
    ap("$EndEntities")
    ap("$Nodes")
    ap(f"1 {n_nodes} 1 {n_nodes}")
    ap(f"2 1 0 {n_nodes}")
    ap("\n".join(str(i + 1) for i in range(n_nodes)))
    ap("\n".join(f"{c[0]:.10g} {c[1]:.10g} {c[2]:.10g}" for c in coords))
    ap("$EndNodes")
    ap("$Elements")
    ap(f"{1 + len(curve_order)} {n_elems} 1 {n_elems}")
    et = 1
    ap(f"2 1 3 {quads.shape[0]}")  # dim2, surface ent1, type3=quad
    buf = []
    for q in quads:
        buf.append(f"{et} {q[0]+1} {q[1]+1} {q[2]+1} {q[3]+1}")
        et += 1
    ap("\n".join(buf))
    for ci, nm in enumerate(curve_order, start=1):
        edges = bedges[nm]
        ap(f"1 {ci} 1 {len(edges)}")  # dim1, curve ent, type1=line
        buf = []
        for a, b in edges:
            buf.append(f"{et} {a+1} {b+1}")
            et += 1
        ap("\n".join(buf))
    ap("$EndElements")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
