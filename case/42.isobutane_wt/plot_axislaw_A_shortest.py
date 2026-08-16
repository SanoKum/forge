"""最短ロバスト軸則 A の軸中心Mach線図と逆MOC Machコンターを生成する。

再生成:
  design/.venv-opt/bin/python case/42.isobutane_wt/plot_axislaw_A_shortest.py
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


CASE = Path(__file__).resolve().parent
RESULT = CASE / "optimize_axislaw_A_shortest.json"
AXIS_PNG = CASE / "best_axislaw_A_shortest_axis_mach.png"
CONTOUR_PNG = CASE / "best_axislaw_A_shortest_mach_contour.png"
FIELD_NPZ = CASE / "best_axislaw_A_shortest_moc_field.npz"

_CMP_PATH = CASE / "compare_axislaw_D.py"
_SPEC = importlib.util.spec_from_file_location("compare_axislaw_D_for_shortest_plot", _CMP_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"比較モジュールを読み込めない: {_CMP_PATH}")
_CMP = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CMP)


def _best_law():
    data = json.loads(RESULT.read_text())
    winner = data["winner"]
    law = _CMP.KnotQuinticAxisLaw(
        _CMP.X_A, winner["L_c"], _CMP.M_A, _CMP.MP_A, _CMP.MPP_A,
        _CMP.MD, winner["M_K"],
    )
    return winner, law


def _calculate_field(law, n_axis: int = 1200):
    x_end = law.x_E + 2.3 * _CMP.RF * np.sqrt(_CMP.MD * _CMP.MD - 1.0)
    return _CMP.inverse_design(
        _CMP.ht, lambda x: float(law(max(x, _CMP.X_A))),
        x_axis_end=x_end, n_axis=n_axis, n_start=41, gamma=_CMP.gas,
        M_start=1.05, exit_mode="characteristic", x_E=law.x_E, M_d=_CMP.MD,
        start_line="throat_char", wall_mode="cplus", axis_dx0=_CMP.AXIS_DX0,
    )


def _axis_plot(law, x_f: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(_CMP.X_A, x_f, 2400)
    mach = np.asarray(law(x), dtype=float)
    fig, ax = plt.subplots(figsize=(10.5, 4.8), constrained_layout=True)
    ax.plot(x, mach, color="#1f6f8b", lw=2.5, label="axis target M(x)")
    ax.axvline(law.x_K, color="#d66b3d", lw=1.4, ls="--")
    ax.axvspan(law.x_E, x_f, color="#dff3ea", alpha=0.55, label="uniform M=6 region")
    ax.scatter([_CMP.X_A, law.x_K, law.x_E], [_CMP.M_A, law.M_K, _CMP.MD],
               c=["#243b53", "#d66b3d", "#25765a"], s=48, zorder=5)
    for label, xx, mm in (("A", _CMP.X_A, _CMP.M_A), ("K", law.x_K, law.M_K),
                          ("E", law.x_E, _CMP.MD)):
        offset = (8, -24) if label == "E" else (8, 9)
        ax.annotate(f"{label}  ({xx:.2f}, M={mm:.2f})", (xx, mm),
                    xytext=offset, textcoords="offset points", fontsize=9)
    ax.scatter([x_f], [_CMP.MD], facecolors="white", edgecolors="#25765a", s=48, zorder=5)
    ax.annotate(f"F  ({x_f:.2f}, M={_CMP.MD:.2f})", (x_f, _CMP.MD),
                xytext=(-122, -24), textcoords="offset points", fontsize=9)
    ax.set(xlabel=r"$x/r_t$", ylabel="Mach number",
           title=rf"Axis-center Mach — $L_c={law.L_c:.2f}$, $M_K={law.M_K:.2f}$")
    ax.set_xlim(_CMP.X_A - 1.0, x_f + 1.5)
    ax.set_ylim(0.9, 6.25)
    ax.grid(True, color="#dce2ea", lw=0.8)
    ax.legend(frameon=False, loc="lower right")
    fig.savefig(AXIS_PNG, dpi=180, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return x, mach


def _contour_plot(res) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    levels = np.asarray(res["levels"], dtype=float)
    wall = np.asarray(res["wall"], dtype=float)
    # 1200²点をそのままDelaunay化せず、規則的に1/4抽出する。図示解像度には十分で、
    # 元のMOC値を補間前に保持する。壁外へ橋渡しする三角形は後段でmaskする。
    sampled = levels[::4, ::4, :]
    x = sampled[..., 0].ravel()
    r = sampled[..., 1].ravel()
    mach = sampled[..., 4].ravel()
    ok = np.isfinite(x) & np.isfinite(r) & np.isfinite(mach)
    x, r, mach = x[ok], r[ok], mach[ok]
    x_f = float(wall[-1, 0])
    wall_at_point = np.interp(x, wall[:, 0], wall[:, 1], left=np.nan, right=np.nan)
    inside = ((x >= wall[0, 0] - 1e-8) & (x <= x_f + 1e-8) & (r >= -1e-10)
              & np.isfinite(wall_at_point) & (r <= wall_at_point + 1e-8))
    x, r, mach = x[inside], r[inside], mach[inside]

    tri = mtri.Triangulation(x, r)
    tx = x[tri.triangles].mean(axis=1)
    tr = r[tri.triangles].mean(axis=1)
    wall_r = np.interp(tx, wall[:, 0], wall[:, 1], left=np.nan, right=np.nan)
    # 点間隔に対して異常に長いDelaunay辺も、領域外を横切る橋渡しとして除外する。
    vertices = tri.triangles
    edge_max = np.maximum.reduce([
        np.hypot(x[vertices[:, 0]] - x[vertices[:, 1]], r[vertices[:, 0]] - r[vertices[:, 1]]),
        np.hypot(x[vertices[:, 1]] - x[vertices[:, 2]], r[vertices[:, 1]] - r[vertices[:, 2]]),
        np.hypot(x[vertices[:, 2]] - x[vertices[:, 0]], r[vertices[:, 2]] - r[vertices[:, 0]]),
    ])
    tri.set_mask(~np.isfinite(wall_r) | (tr > wall_r + 0.03) | (edge_max > 2.0))

    fig, ax = plt.subplots(figsize=(12.5, 4.2), constrained_layout=True)
    contour_levels = np.linspace(1.0, 6.0, 26)
    cf = ax.tricontourf(tri, mach, levels=contour_levels, cmap="turbo", extend="both")
    ax.plot(wall[:, 0], wall[:, 1], color="#111827", lw=1.7, label="designed wall")
    ax.plot([wall[0, 0], x_f], [0.0, 0.0], color="#111827", lw=0.9)
    ax.axvline(float(res["exit"]["x_F"]), color="white", lw=1.0, ls="--", alpha=0.9)
    cb = fig.colorbar(cf, ax=ax, pad=0.015)
    cb.set_label("Mach number")
    ax.set(xlabel=r"$x/r_t$", ylabel=r"$r/r_t$",
           title="Inverse-MOC Mach contour — best robust-shortest design")
    ax.set_xlim(wall[0, 0], x_f)
    ax.set_ylim(0.0, 1.05 * float(wall[-1, 1]))
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=True, loc="upper left")
    fig.savefig(CONTOUR_PNG, dpi=180, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return x, r, mach, wall


def main() -> None:
    winner, law = _best_law()
    res = _calculate_field(law, n_axis=int(winner["n_axis"]))
    # 結果JSONのx_Fは、診断・順位付けに使った壁テーブル最終点。レポート図も同じ定義へ合わせる。
    axis_x, axis_mach = _axis_plot(law, float(winner["x_F"]))
    field_x, field_r, field_mach, wall = _contour_plot(res)
    np.savez_compressed(
        FIELD_NPZ,
        axis_x=axis_x, axis_mach=axis_mach,
        field_x=field_x, field_r=field_r, field_mach=field_mach,
        wall=wall,
        L_c=float(winner["L_c"]), M_K=float(winner["M_K"]),
        n_axis=int(winner["n_axis"]),
    )
    print(f"saved {AXIS_PNG}")
    print(f"saved {CONTOUR_PNG}")
    print(f"saved {FIELD_NPZ}")


if __name__ == "__main__":
    main()
