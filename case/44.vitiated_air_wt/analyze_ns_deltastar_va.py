#!/usr/bin/env python3
"""A12: NS run から δ*(x) を抽出し相関と比較、v2 用 CSV を出力する。

使い方: design/.venv-opt/bin/python analyze_ns_deltastar.py <run_dir>
出力: <run_dir>/dstar_extracted.csv (x_rt, dstar_rt — v2 の prepare_ns --dstar-csv 入力)
      <run_dir>/dstar_compare.png
"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from matplotlib import font_manager as _fm  # noqa: E402

for _f in ("/home/sano/.fonts/NotoSansCJKjp-Regular.otf",):
    try:
        _fm.fontManager.addfont(_f)
    except Exception:  # noqa: BLE001
        pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "design"))

from forge_design.feedback.deltastar import deltastar_offset  # noqa: E402
from forge_design.metrics.deltastar import deltastar_from_run  # noqa: E402

SCALE = None  # prepare_info.json の scale_m から取る


def main(run_dir):
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    info = json.loads((rd / "prepare_info.json").read_text())
    global SCALE; SCALE = float(info["scale_m"])
    wall_phys = np.loadtxt(rd / "wall_physical.csv", delimiter=",", skiprows=1) / SCALE
    wall_inv = np.loadtxt(rd / "wall_design.csv", delimiter=",", skiprows=1)
    wall_inv = wall_inv / [SCALE, SCALE, 1.0, 1.0]
    x_E = float(info["x_E"])
    # 測定ステーション: x∈[8, x_F−0.5] (x<8 はコア未一様でエッジ検出が破綻 — run_0030)
    xs = np.arange(8.0, wall_phys[-1, 0] - 0.5, 0.25)
    d = deltastar_from_run(rd / "nozzle.h5", res, wall_phys, SCALE, xs)
    ok = d["ok"] & np.isfinite(d["dstar"])
    # 相関 (v1 で使った壁オフセット量) との比較: phys − inv の法線距離
    ds_corr = np.interp(xs, wall_inv[:, 0],
                        (np.interp(wall_inv[:, 0], wall_phys[:, 0], wall_phys[:, 1])
                         - wall_inv[:, 1]) * np.cos(wall_inv[:, 2]))
    ratio = d["dstar"][ok] / np.maximum(ds_corr[ok], 1e-12)
    print(f"measured {ok.sum()}/{len(xs)} stations")
    print(f"δ*_CFD/δ*_corr: 中央値 {np.median(ratio):.3f}  範囲 "
          f"[{np.min(ratio):.3f}, {np.max(ratio):.3f}]")
    for i in range(0, ok.sum(), max(ok.sum() // 8, 1)):
        j = np.where(ok)[0][i]
        print(f"  x={xs[j]:6.2f}  δ*_CFD={d['dstar'][j]:.5f}  corr={ds_corr[j]:.5f}  "
              f"比 {d['dstar'][j]/max(ds_corr[j],1e-12):.3f}")
    # v2 用 CSV: 測定域は CFD、x<8 は相関へ滑らかにブレンド。
    # 測定域の外は端値クリップでなく**端勾配の線形外挿** (クリップは折れ目を作り、
    # 壁 B-spline が下流端で非単調になる — run_0071 初回で実測)。最後に全域を
    # もう一度平滑化してブレンド点の C1 折れも均す。
    from scipy.interpolate import make_smoothing_spline
    spl = make_smoothing_spline(xs[ok], d["dstar"][ok], lam=1e-2)
    x_full = wall_inv[:, 0]
    ds_c_full = (np.interp(x_full, wall_phys[:, 0], wall_phys[:, 1])
                 - wall_inv[:, 1]) * np.cos(wall_inv[:, 2])
    x_lo, x_hi = xs[ok][0], xs[ok][-1]
    ds_cfd = np.clip(spl(np.clip(x_full, x_lo, x_hi)), 0.0, None)
    slope_hi = float(spl.derivative()(x_hi))
    m_hi = x_full > x_hi
    ds_cfd[m_hi] = np.clip(float(spl(x_hi)) + slope_hi * (x_full[m_hi] - x_hi), 0.0, None)
    w = np.clip((x_full - 6.0) / (9.0 - 6.0), 0.0, 1.0)
    w = w * w * (3 - 2 * w)
    ds_v2 = (1 - w) * ds_c_full + w * ds_cfd
    ds_v2 = make_smoothing_spline(x_full, ds_v2, lam=1e-3)(x_full)
    ds_v2 = np.clip(ds_v2, 0.0, None)
    np.savetxt(rd / "dstar_extracted.csv", np.c_[x_full, ds_v2], delimiter=",",
               header="x_rt,dstar_rt (CFD x>9 / corr x<6 / blend)", comments="")
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(xs[ok], d["dstar"][ok], "o", ms=3, label="δ* CFD 抽出")
    ax.plot(x_full, ds_c_full, "-", label="δ* 相関 (v1 で使用)")
    ax.plot(x_full, ds_v2, "--", label="v2 壁へ渡す δ*")
    ax.axvline(x_E, color="gray", lw=0.7, ls=":")
    ax.set_xlabel("x / r_t"); ax.set_ylabel("δ* / r_t")
    ax.legend(); ax.set_title(f"{rd.name}: 排除厚 (中央値比 {np.median(ratio):.3f})")
    fig.tight_layout(); fig.savefig(rd / "dstar_compare.png", dpi=110)
    print(f"saved {rd/'dstar_extracted.csv'}")


if __name__ == "__main__":
    main(sys.argv[1])
