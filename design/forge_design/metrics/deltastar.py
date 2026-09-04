r"""RANS 場からの排除厚 $\delta^*(x)$ 抽出 (A12 = axis-Mach 粘性補正)。

計画: plans/active/tooling-nozzle-axismach-viscous-deltastar.md §4。

定義 (圧縮性・**軸対称の質量収支** — 2026-09-01 に平面近似から置換):

$$\delta^*\left(1-\frac{\delta^*}{2r_w}\right) = \int_0^{y_e}
  \left(1 - \frac{\rho u}{(\rho u)_e}\right)\left(1-\frac{y}{r_w}\right) dy
  \;\Rightarrow\; \delta^* = r_w\left(1-\sqrt{1-2I/r_w}\right)$$

壁から $y$ 入った欠損は半径 $r_w-y$ の円環を流れるため $(1-y/r_w)$ の重みが付く。
旧・平面近似 ($\delta^*=\int(1-\rho u/(\rho u)_e)dy$) は欠損重心が壁近くに集中する
ため過大評価は小さい (M6 出口 $\delta/r_w\approx0.14$ で **+1.4 %**、出口 M 換算
+0.05 % — case/45 run_0007 実測)。帳簿の正確さのため軸対称形を正とする。

**教訓 (run_0030 で実測、繰り返し禁止)**:
- 探索窓はフリーストリームに届くまで取る (窓 0.35 r* は下流で BL 縁の手前で切れ、
  「δ* が下流ほど減る」という非物理な結論を出した — 撤回済み)。既定 1.5 r*。
- スロート近傍 (x ≲ 8 r_t) はコア流が半径方向に未一様で「ρu 極大 = BL 縁」の
  判定が破綻する。呼び出し側で相関へフォールバックすること。
"""
from __future__ import annotations

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def deltastar_from_run(mesh_h5, res_h5, wall_xy, scale: float,
                       x_stations, window: float = 1.5, n_prof: int = 400,
                       edge_frac: float = 0.995) -> dict:
    """node RANS run から各 x ステーションの δ* を測る。

    wall_xy: (n,2) 物理壁 [m 単位 or r_t 単位 — scale と整合していれば可]。
    x_stations: r_t 単位。戻り: dict(x, dstar, ue_rel, y_edge, ok) 全て r_t 単位。
    """
    with h5py.File(mesh_h5) as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(res_h5) as f:
        ro = f["/VALUE/ro"][:]
        Ux = f["/VALUE/Ux"][:]
    if len(ro) != len(nc):
        raise ValueError("node run でない (VALUE 長 != 節点数)")
    x_n, r_n = nc[:, 0] / scale, nc[:, 1] / scale
    itp = LinearNDInterpolator(np.c_[x_n, r_n], ro * Ux)
    wx, wr = np.asarray(wall_xy)[:, 0], np.asarray(wall_xy)[:, 1]
    out = {"x": [], "dstar": [], "y_edge": [], "ok": []}
    for xs in np.atleast_1d(x_stations):
        rw = float(np.interp(xs, wx, wr))
        y = np.linspace(0.0, min(window, 0.95 * rw), n_prof)   # 壁 → 内側
        q = itp(np.full(n_prof, xs), rw - y)
        good = np.isfinite(q)
        if good.sum() < n_prof // 2:
            out["x"].append(float(xs)); out["dstar"].append(np.nan)
            out["y_edge"].append(np.nan); out["ok"].append(False)
            continue
        q = np.where(good, q, 0.0)
        qe = float(np.max(q))
        # 縁 = 壁から見て最初に edge_frac·(ρu)_max へ達する y
        idx = np.argmax(q >= edge_frac * qe)
        y_e = float(y[idx])
        m = y <= y_e
        # 軸対称の質量収支: 欠損に円環重み (1 - y/r_w) を掛け、
        # δ*(1 - δ*/(2 r_w)) = I を δ* について解く (docstring 参照)
        d = 1.0 - q[m] / max(qe, 1e-30)
        I = float(np.trapezoid(d * (1.0 - y[m] / rw), y[m]))
        integ = float(rw * (1.0 - np.sqrt(max(1.0 - 2.0 * I / rw, 0.0))))
        out["x"].append(float(xs)); out["dstar"].append(integ)
        out["y_edge"].append(y_e); out["ok"].append(bool(idx > 3))
    return {k: np.asarray(v) for k, v in out.items()}


# --- 固定 Euler 基準のコア整合質量欠損抽出 (plans/active/tooling-nozzle-deltastar-core-matched-euler.md) ----
def _load_structured(run_dir):
    """node 構造格子 run (index = i*nj + j) を (info, x[ni,nj], r[ni,nj], q=ρUx[ni,nj]) で返す (r_t 単位)。"""
    import json
    from pathlib import Path
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S = float(info["scale_m"])
    with h5py.File(run_dir / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1]
    with h5py.File(res) as f:
        q = f["/VALUE/ro"][:] * f["/VALUE/Ux"][:]
    if len(q) != ni * nj or len(nc) != ni * nj:
        raise ValueError(f"{run_dir}: node 構造格子 run でない (VALUE {len(q)}, COORD {len(nc)}, ni*nj {ni*nj})")
    return dict(info=info, S=S, res=res.name,
                x=nc[:, 0].reshape(ni, nj) / S, r=nc[:, 1].reshape(ni, nj) / S, q=q.reshape(ni, nj))


def massflow_ratio(ns_run, euler_run) -> dict:
    """NS/Euler 質量流量比 → 有効スロート半径・質量流量由来の等価スロート補正量 (r_t 単位)。
    列ごとに 2π∫ρUx r dr (台形則)、x∈(−2.5, x_E) の中央値を採用。"""
    out = {}
    for tag, rd in (("ns", ns_run), ("euler", euler_run)):
        d = _load_structured(rd)
        S = d["S"]
        m = 2 * np.pi * np.trapezoid(d["q"] * d["r"] * S * S, d["r"] * S, axis=1)
        xs = d["x"][:, 0]
        sel = (xs > -2.5) & (xs < d["info"]["x_E"])
        out[f"mdot_{tag}"] = float(np.median(m[sel]))
        out[f"mdot_{tag}_spread"] = float(m[sel].std() / np.median(m[sel]))
        out[f"res_{tag}"] = d["res"]
        if tag == "ns":
            r_tW = float(d["info"].get("throat_physical", {}).get("r", np.nan)) if d["info"].get("throat_physical") else float(np.min(d["r"][:, -1]))
    ratio = out["mdot_ns"] / out["mdot_euler"]
    out["mdot_ratio"] = float(ratio)
    out["r_throat_eff"] = float(np.sqrt(ratio))            # 設計 r_t 単位
    out["r_throat_phys"] = r_tW
    out["cd_rel"] = float(ratio / r_tW ** 2)
    out["delta_r_throat_eff"] = float(r_tW - np.sqrt(ratio))
    out["delta_r_throat_given"] = float(r_tW - 1.0)
    return out



def core_matched_deficit(r, q_ns, q_e, rw_e: float, core_frac: float = 0.30,
                         n_axis_skip: int = 1, outer_frac: float = 0.25, refine: int = 4):
    r"""1 断面の**コア整合質量欠損 → 半径方向等価排除厚** (純関数、単体試験対象)。

    r: 軸→壁の半径 (NS 格子, 最後が NS 壁)、q_ns: NS の ρUx、q_e: 同じ r に写した Euler の ρUx
    (Euler 壁 rw_e より外は壁値) — 配列 (r と同長) か callable(r)。
    refine>1 で NS 区間を細分 (q_ns は区分線形、q_e は callable なら細分点で評価、配列なら線形補間)。
    壁集中メッシュはコアに節点が ~5 個しか無いので、積分と α を安定させるために使う。
    戻り: dict または None (コア点数不足)。"""
    r0 = np.asarray(r, dtype=float); q0 = np.asarray(q_ns, dtype=float)
    if refine > 1:
        t = np.linspace(0.0, 1.0, refine, endpoint=False)
        r = np.concatenate([(r0[:-1, None] + (r0[1:] - r0[:-1])[:, None] * t[None, :]).ravel(), [r0[-1]]])
        qN = np.interp(r, r0, q0)
        qE = np.asarray(q_e(r), dtype=float) if callable(q_e) else np.interp(r, r0, np.asarray(q_e, dtype=float))
    else:
        r, qN = r0, q0
        qE = np.asarray(q_e(r), dtype=float) if callable(q_e) else np.asarray(q_e, dtype=float)
    rwN = float(r[-1])
    mc = r <= core_frac * rw_e
    if mc.sum() < n_axis_skip + 3:
        return None
    num = np.trapezoid(qN[mc] * r[mc], r[mc]); den = np.trapezoid(qE[mc] * r[mc], r[mc])
    alpha = float(num / den) if den > 0 else np.nan
    qref = alpha * qE
    rel = qN[mc] / np.maximum(qref[mc], 1e-30) - 1.0
    core_rms = float(np.sqrt(np.mean(rel ** 2)))
    core_rms_noaxis = float(np.sqrt(np.mean(rel[n_axis_skip:] ** 2)))
    core_maxdev = float(np.max(np.abs(rel[n_axis_skip:])))
    integrand = (qref - qN) * r
    D = float(2 * np.pi * np.trapezoid(integrand, r))
    # 壁からの累積 F(y) = 2π∫_{rw-y}^{rw} q_ref r' dr' (y = 壁からの距離, 単調増加)
    seg = 0.5 * (qref[1:] * r[1:] + qref[:-1] * r[:-1]) * np.diff(r)
    F_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    r_desc = r[::-1]
    reason = []
    if D < 0:
        q_scale = max(float(np.mean(qref[mc])), 1e-30)
        delta_r = D / (2 * np.pi * rwN * q_scale)                          # 線形化 (負値, 診断用)
        reason.append("negative_deficit")
    elif D > F_from_wall[-1]:
        delta_r = float("nan"); reason.append("no_root")
    else:
        r_eff = float(np.interp(D, F_from_wall, r_desc))
        delta_r = rwN - r_eff
    segD = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)
    G_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(segD[::-1])])
    y_desc = rwN - r_desc
    G_outer = float(np.interp(outer_frac * rwN, y_desc, G_from_wall))
    outer_share = G_outer / D if D > 0 else float("nan")
    return dict(r_wall_euler=float(rw_e), r_wall_ns=rwN, delta_in=rwN - float(rw_e), alpha=alpha,
                core_rms=core_rms, core_rms_noaxis=core_rms_noaxis, core_maxdev=core_maxdev,
                mass_deficit=D, delta_r=float(delta_r), outer_share=outer_share, reason=reason)


def deltastar_from_core_matched_euler(ns_run, euler_run, core_frac: float = 0.30,
                                      core_frac_sens=(0.25, 0.35), n_axis_skip: int = 1,
                                      outer_frac: float = 0.25, smooth_lam: float = 1e-3,
                                      gate_core_rms: float = 0.01, gate_sens: float = 0.05,
                                      gate_sens_floor: float = 1e-3, refine: int = 4, out_dir=None) -> dict:
    r"""**固定 Euler 基準・コア整合**の半径方向等価排除厚 $\delta_r(x)$ を NS 全列で抽出する。

    定義 (plan §4.2–4.4):

    - 対応は**同じ物理 $r$** (設計 $r_t$ 単位)。NS の有効壁 $r_{ref}=r_{w,NS}-\delta_{in}$ は、
      $\delta_{in}$ を「NS 壁 − Euler 壁」で定義するので $r_{ref}=r_{w,E}$ となり、η 写像は恒等になる。
      Euler 壁の外側 ($r>r_{w,E}$) は Euler 壁値で定数外挿。
    - コア $r \le \eta_c r_{w,E}$ の面積重み流量を一致させる倍率 $\alpha$、$q_{ref}=\alpha q_E$。
    - 符号付き欠損 $D = 2\pi\int_0^{r_{w,NS}}(q_{ref}-q_{NS})\,r\,dr$ (クリップなし)。
      $2\pi\int_{r_{eff}}^{r_{w,NS}} q_{ref}\,r\,dr = D$ の根から $\delta_r = r_{w,NS}-r_{eff}$。
    - ゲート (診断; 値は捨てない): コア相対 RMS、コア範囲感度、$D<0$、根なし、欠損のコア分布。

    戻り値: dict of arrays (x, r_wall_euler, r_wall_ns, delta_in, alpha, core_rms, core_rms_noaxis,
    core_maxdev, mass_deficit, delta_r_raw, delta_r_smooth, delta_r_sens (n,2), ok, reason) と mdot 帳簿。
    `out_dir` を与えると `delta_r_equiv.csv` / `delta_r_equiv_diag.json` を書く。"""
    import json
    from pathlib import Path
    from scipy.interpolate import make_smoothing_spline
    N = _load_structured(ns_run)
    E = _load_structured(euler_run)
    xE = E["x"][:, 0]
    rwE_col = E["r"][:, -1]

    def euler_q_at(x, r):
        """Euler 場を (x, r) へ: 列間線形 × 列内線形 (r>壁は壁値)。"""
        k = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
        w = (x - xE[k]) / max(xE[k + 1] - xE[k], 1e-30)
        w = float(np.clip(w, 0.0, 1.0))
        q0 = np.interp(r, E["r"][k], E["q"][k])          # np.interp は範囲外を端値で埋める
        q1 = np.interp(r, E["r"][k + 1], E["q"][k + 1])
        rw = (1 - w) * rwE_col[k] + w * rwE_col[k + 1]
        return (1 - w) * q0 + w * q1, rw

    def extract_one(i, fc):
        x = float(N["x"][i, 0]); r0 = N["r"][i]; q0 = N["q"][i]
        if x < xE[0] - 1e-9 or x > xE[-1] + 1e-9:
            return None
        _, rwE = euler_q_at(x, r0[:1])
        res = core_matched_deficit(r0, q0, lambda rr: euler_q_at(x, rr)[0], float(rwE), core_frac=fc,
                                   n_axis_skip=n_axis_skip, outer_frac=outer_frac, refine=refine)
        if res is None:
            return None
        res["x"] = x
        return res

    rows = []
    for i in range(N["x"].shape[0]):
        base = extract_one(i, core_frac)
        if base is None:
            continue
        sens = [extract_one(i, fc) for fc in core_frac_sens]
        base["delta_r_sens"] = [s["delta_r"] if s else np.nan for s in sens]
        rows.append(base)
    x = np.array([r_["x"] for r_ in rows])
    draw = np.array([r_["delta_r"] for r_ in rows])
    sens = np.array([r_["delta_r_sens"] for r_ in rows])
    HARD = ("negative_deficit", "no_root", "nan")
    ok = np.ones(len(rows), bool); hard_ok = np.ones(len(rows), bool); reasons = []
    for k, r_ in enumerate(rows):
        rs = list(r_["reason"])
        if r_["core_rms_noaxis"] > gate_core_rms: rs.append("core_shape")
        if np.isfinite(draw[k]) and draw[k] > 0 and np.all(np.isfinite(sens[k])):
            if np.max(np.abs(sens[k] - draw[k])) > max(gate_sens * draw[k], gate_sens_floor): rs.append("core_frac_sens")
        if np.isfinite(r_["outer_share"]) and r_["outer_share"] < 0.8: rs.append("deficit_not_near_wall")
        if not np.isfinite(draw[k]): rs.append("nan")
        ok[k] = len(rs) == 0
        hard_ok[k] = not any(h in rs for h in HARD)
        reasons.append(",".join(rs))
    # 平滑化: hard-ok の生値だけを使い、範囲外は端値 (外挿しない)。soft 不合格の値は使う (plan §4.4)。
    dsm = np.full_like(draw, np.nan)
    if hard_ok.sum() > 10:
        xs_, ds_ = x[hard_ok], draw[hard_ok]
        spl = make_smoothing_spline(xs_, ds_, lam=smooth_lam)
        dsm = spl(np.clip(x, xs_[0], xs_[-1]))
    # 壁更新に渡す値: hard-ok なら平滑化値、hard 不合格なら NaN (呼び出し側が前回値を保持する)
    duse = np.where(hard_ok, dsm, np.nan)
    out = dict(x=x, r_wall_euler=np.array([r_["r_wall_euler"] for r_ in rows]),
               r_wall_ns=np.array([r_["r_wall_ns"] for r_ in rows]),
               delta_in=np.array([r_["delta_in"] for r_ in rows]),
               alpha=np.array([r_["alpha"] for r_ in rows]),
               core_rms=np.array([r_["core_rms"] for r_ in rows]),
               core_rms_noaxis=np.array([r_["core_rms_noaxis"] for r_ in rows]),
               core_maxdev=np.array([r_["core_maxdev"] for r_ in rows]),
               mass_deficit=np.array([r_["mass_deficit"] for r_ in rows]),
               outer_share=np.array([r_["outer_share"] for r_ in rows]),
               delta_r_raw=draw, delta_r_smooth=dsm, delta_r_use=duse, delta_r_sens=sens,
               ok=ok, hard_ok=hard_ok, reason=np.array(reasons),
               settings=dict(core_frac=core_frac, core_frac_sens=list(core_frac_sens), n_axis_skip=n_axis_skip,
                             outer_frac=outer_frac, smooth_lam=smooth_lam, gate_core_rms=gate_core_rms,
                             gate_sens=gate_sens, gate_sens_floor=gate_sens_floor, refine=refine,
                             ns_res=N["res"], euler_res=E["res"],
                             ns_run=str(ns_run), euler_run=str(euler_run)))
    out["massflow"] = massflow_ratio(ns_run, euler_run)
    if out_dir is not None:
        od = Path(out_dir); od.mkdir(parents=True, exist_ok=True)
        np.savetxt(od / "delta_r_equiv.csv",
                   np.c_[x, draw, dsm, duse, out["delta_in"], out["alpha"], out["core_rms_noaxis"], sens,
                         ok.astype(int), hard_ok.astype(int)],
                   delimiter=",", comments="",
                   header="x_rt,delta_r_raw,delta_r_smooth,delta_r_use,delta_in,alpha,core_rms_noaxis,delta_r_c25,delta_r_c35,ok,hard_ok")
        diag = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in out.items()}
        (od / "delta_r_equiv_diag.json").write_text(json.dumps(diag))
    return out
