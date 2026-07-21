#!/usr/bin/env python3
"""線形場に対する勾配・面再構成ジャンプの厳密性検証 (GG vs LSQ)。

背景 (`keepDissJump>=1`, plans/accepted/convection-keep-diss-recon-jump.md):
再構成後ジャンプ Δq_rec = (q_j + g_j·(x_f−x_j)) − (q_i + g_i·(x_f−x_i)) が滑らかな場で
O(h³) に落ちるのは **勾配が線形場を厳密に再現する場合のみ**。勾配に O(1) 相対誤差が
あれば Δq_rec は O(h) に留まり、勾配誤差が「格子品質依存の不均一な散逸」として流束に
戻る。本ツールは任意線形場 q = a + b·x に対し
  - 勾配誤差 max|g − b|/|b|
  - 内部面の再構成ジャンプ max|Δq_rec|/(|b|·h)   (h = セル中心間距離; 線形場の真値は 0)
を、forge の式を忠実に写像した GG / LSQ で測る。

写像した式 (いずれも solver 実装と同一):
  - GG: 面値 φf = fx·φ_ic0 + (1−fx)·φ_ic1、fx = clamp(((cc1−pc)·(cc1−cc0))/|cc1−cc0|², 0, 1)
    (`setStructualVariables.cpp`)、勾配 = Σ φf S / V (`calcGradient_cellgather_d`)。
    境界面は厳密境界値 φ(pc) を与える (bvar 相当の理想化 = GG に有利な設定)。
  - LSQ: 逆距離二乗重み w=1/d² の正規方程式 M g = b (`lsqGrad_accumInternal_d` /
    `lsqGrad_accumBoundary_d`)。境界は面重心の厳密値を LSQ 点に追加。
    M/b の蓄積は double だが格納配列は flow_float (float32) のため、float32 格納で
    solve する変異版も測る (plans/active/discretization-lsq-gradient.md §9 の近壁退化の再現)。

メッシュ: 2D 四角形 (一様 / ノードジッタ / 三角形分割 / 高アスペクト比+ジッタ)。

期待 (= 判定):
  - GG 一様直交: 機械精度ゼロ (線形場厳密)。
  - GG ジッタ/三角形: 非ゼロ (O(1) 相対勾配誤差 → Δq_rec は O(h))。**これが keepDissJump の
    「高波数選択性」が非一様メッシュで劣化する定量根拠** (характер確認であり FAIL ではない)。
  - LSQ (double): 全メッシュで機械精度ゼロ (線形場厳密 — 重みによらず成立)。
  - LSQ (float32 格納): double 比で 6 桁以上劣化 (~1e-6)。plan §9 の 1e9 級スパイクの
    完全再現には近特異幾何 (近壁の歪み+近傍不足) が必要で、本合成メッシュでは劣化傾向のみを示す。

実行: python3 verify_linear_recon.py   (依存: numpy のみ)
"""
import numpy as np


def build_quad_mesh(nx, ny, jitter=0.0, aspect=1.0, tri=False, seed=1):
    """2D 四角形 (または対角分割三角形) メッシュを構築。

    戻り値: cells (ノード index list), nodes (N,2), 内部面 [(ic0,ic1,(n0,n1))], 境界面 [(ic,(n0,n1))]
    """
    rng = np.random.default_rng(seed)
    xs = np.linspace(0.0, 1.0, nx+1)
    ys = np.linspace(0.0, 1.0/aspect, ny+1)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nodes = np.stack([X.ravel(), Y.ravel()], axis=-1)
    if jitter > 0:
        hx, hy = xs[1]-xs[0], ys[1]-ys[0]
        pert = rng.uniform(-jitter, jitter, nodes.shape) * np.array([hx, hy])
        interior = ((nodes[:, 0] > xs[0]+1e-12) & (nodes[:, 0] < xs[-1]-1e-12) &
                    (nodes[:, 1] > ys[0]+1e-12) & (nodes[:, 1] < ys[-1]-1e-12))
        nodes[interior] += pert[interior]

    def nid(i, j):
        return i*(ny+1) + j

    cells = []
    for i in range(nx):
        for j in range(ny):
            q = [nid(i, j), nid(i+1, j), nid(i+1, j+1), nid(i, j+1)]
            if tri:
                if (i+j) % 2 == 0:
                    cells.append([q[0], q[1], q[2]]); cells.append([q[0], q[2], q[3]])
                else:
                    cells.append([q[0], q[1], q[3]]); cells.append([q[1], q[2], q[3]])
            else:
                cells.append(q)

    # エッジ → セル対
    edge_cells = {}
    for ic, cell in enumerate(cells):
        m = len(cell)
        for k in range(m):
            e = (cell[k], cell[(k+1) % m])
            key = (min(e), max(e))
            edge_cells.setdefault(key, []).append(ic)
    internal = [(cs[0], cs[1], key) for key, cs in edge_cells.items() if len(cs) == 2]
    boundary = [(cs[0], key) for key, cs in edge_cells.items() if len(cs) == 1]
    return cells, nodes, internal, boundary


def cell_geom(cells, nodes):
    """多角形の面積と面積重心 (forge のセル重心相当)"""
    cc = np.zeros((len(cells), 2))
    area = np.zeros(len(cells))
    for ic, cell in enumerate(cells):
        p = nodes[cell]
        x, y = p[:, 0], p[:, 1]
        xn, yn = np.roll(x, -1), np.roll(y, -1)
        cr = x*yn - xn*y
        A = 0.5*np.sum(cr)
        cx = np.sum((x+xn)*cr)/(6*A)
        cy = np.sum((y+yn)*cr)/(6*A)
        area[ic] = A
        cc[ic] = (cx, cy)
    return cc, area


def face_geom(nodes, e, cc0=None):
    """エッジの重心・法線ベクトル (長さ込み)。cc0 指定時は外向きに正規化。"""
    p0, p1 = nodes[e[0]], nodes[e[1]]
    pc = 0.5*(p0+p1)
    t = p1 - p0
    S = np.array([t[1], -t[0]])          # 90° 回転
    if cc0 is not None and np.dot(S, pc - cc0) < 0:
        S = -S
    return pc, S


def gg_gradients(cells, nodes, cc, area, internal, boundary, qfun):
    """forge GG (cellgather + fx 射影補間 + 厳密境界値) の写像"""
    q = np.array([qfun(c) for c in cc])
    g = np.zeros((len(cells), 2))
    for ic0, ic1, e in internal:
        pc, S = face_geom(nodes, e, cc[ic0])          # ic0 の外向き
        dcc = cc[ic1] - cc[ic0]
        fx = np.clip(np.dot(cc[ic1]-pc, dcc)/np.dot(dcc, dcc), 0.0, 1.0)
        qf = fx*q[ic0] + (1-fx)*q[ic1]
        g[ic0] += S*qf
        g[ic1] -= S*qf
    for ic, e in boundary:
        pc, S = face_geom(nodes, e, cc[ic])
        g[ic] += S*qfun(pc)                            # 理想化 bvar (厳密境界値)
    return g/area[:, None], q


def lsq_gradients(cells, nodes, cc, internal, boundary, qfun, storage=np.float64):
    """forge LSQ (逆距離二乗重み正規方程式) の写像。storage=float32 で flow_float 格納を模擬。"""
    q = np.array([qfun(c) for c in cc])
    nC = len(cells)
    M = np.zeros((nC, 2, 2))
    b = np.zeros((nC, 2))
    neigh = [[] for _ in range(nC)]
    for ic0, ic1, _ in internal:
        neigh[ic0].append(ic1); neigh[ic1].append(ic0)
    for ic in range(nC):
        for jc in neigh[ic]:
            d = cc[jc] - cc[ic]
            w = 1.0/max(np.dot(d, d), 1e-300)
            M[ic] += w*np.outer(d, d)
            b[ic] += w*d*(q[jc]-q[ic])
    for ic, e in boundary:
        pc, _ = face_geom(nodes, e)
        d = pc - cc[ic]
        w = 1.0/max(np.dot(d, d), 1e-300)
        M[ic] += w*np.outer(d, d)
        b[ic] += w*d*(qfun(pc)-q[ic])
    M = M.astype(storage).astype(np.float64)   # 格納精度で丸め、solve は double (kernel と同順)
    b = b.astype(storage).astype(np.float64)
    g = np.zeros((nC, 2))
    for ic in range(nC):
        g[ic] = np.linalg.solve(M[ic] + 1e-300*np.eye(2), b[ic])
    return g, q


def recon_jump(cc, internal, nodes, g, q, qfun):
    """内部面の再構成ジャンプ |Δq_rec| を h=|cc1-cc0| で正規化した最大値 (線形場の真値 0)"""
    worst = 0.0
    for ic0, ic1, e in internal:
        pc, _ = face_geom(nodes, e)
        qL = q[ic0] + np.dot(g[ic0], pc - cc[ic0])
        qR = q[ic1] + np.dot(g[ic1], pc - cc[ic1])
        h = np.linalg.norm(cc[ic1] - cc[ic0])
        worst = max(worst, abs(qR - qL)/h)
    return worst


def run_case(name, nx, ny, jitter, aspect, tri, checks):
    bvec = np.array([0.7, -0.4])
    a0 = 0.3

    def qfun(x):
        return a0 + np.dot(bvec, x)

    cells, nodes, internal, boundary = build_quad_mesh(nx, ny, jitter, aspect, tri)
    cc, area = cell_geom(cells, nodes)
    nb = np.linalg.norm(bvec)

    gg, q = gg_gradients(cells, nodes, cc, area, internal, boundary, qfun)
    e_gg = np.max(np.linalg.norm(gg - bvec, axis=1))/nb
    j_gg = recon_jump(cc, internal, nodes, gg, q, qfun)/nb

    lsq, q2 = lsq_gradients(cells, nodes, cc, internal, boundary, qfun)
    e_lsq = np.max(np.linalg.norm(lsq - bvec, axis=1))/nb
    j_lsq = recon_jump(cc, internal, nodes, lsq, q2, qfun)/nb

    lsqf, q3 = lsq_gradients(cells, nodes, cc, internal, boundary, qfun, storage=np.float32)
    e_lsqf = np.max(np.linalg.norm(lsqf - bvec, axis=1))/nb

    print(f"--- {name} ({len(cells)} cells) ---")
    print(f"  GG  : 勾配誤差 max|g−b|/|b| = {e_gg:.3e}   再構成ジャンプ max|Δq_rec|/(|b|h) = {j_gg:.3e}")
    print(f"  LSQ : 勾配誤差             = {e_lsq:.3e}   再構成ジャンプ             = {j_lsq:.3e}")
    print(f"  LSQ (float32 格納 M/b)     = {e_lsqf:.3e}")

    ok = True
    for label, val, lo, hi in checks(e_gg, j_gg, e_lsq, j_lsq, e_lsqf):
        r = lo <= val <= hi
        ok &= r
        print(f"    [{label}] {val:.3e} ∈ [{lo:.0e}, {hi:.0e}] ... {'PASS' if r else 'FAIL'}")
    return ok


def main():
    ok = True
    # 一様直交: GG も線形場厳密のはず
    ok &= run_case("uniform quad 40x40", 40, 40, 0.0, 1.0, False,
                   lambda eg, jg, el, jl, ef: [
                       ("GG 線形厳密", eg, 0, 1e-10),
                       ("LSQ 線形厳密", el, 0, 1e-10)])
    # ノードジッタ: GG は非厳密 (O(1) 相対誤差)、LSQ(double) は厳密のまま
    ok &= run_case("jittered quad 40x40 (30% jitter)", 40, 40, 0.3, 1.0, False,
                   lambda eg, jg, el, jl, ef: [
                       ("GG 非厳密 (勾配誤差が残る)", eg, 1e-4, 1.0),
                       ("GG ジャンプ O(h) 残留", jg, 1e-4, 1.0),
                       ("LSQ 線形厳密", el, 0, 1e-10),
                       ("LSQ ジャンプ機械ゼロ", jl, 0, 1e-10)])
    # 三角形分割: 同上 (非構造の代表)
    ok &= run_case("triangulated 40x40", 40, 40, 0.15, 1.0, True,
                   lambda eg, jg, el, jl, ef: [
                       ("GG 非厳密", eg, 1e-4, 1.0),
                       ("LSQ 線形厳密", el, 0, 1e-10)])
    # 高アスペクト比 (AR=100) + ジッタ: LSQ float32 格納の退化を式レベルで再現
    ok &= run_case("high-AR 100 jittered 40x40", 40, 40, 0.2, 100.0, False,
                   lambda eg, jg, el, jl, ef: [
                       ("LSQ(double) は高 AR でも厳密", el, 0, 1e-8),
                       ("LSQ(float32 格納) は double 比で桁劣化", ef, 1e-7, np.inf)])
    print()
    print(f"VERDICT: {'PASS' if ok else 'FAIL'} — GG は非一様メッシュで線形場非厳密 "
          f"(keepDissJump の O(h³) 選択性は Cartesian 限定)、LSQ は double では線形場厳密だが "
          f"float32 格納で高 AR 退化 (lsq plan §9 と整合)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
