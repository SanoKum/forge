"""最短ロバスト軸長探索の自己完結型HTML結果ページを生成する。

入力:
  optimize_axislaw_A_shortest.json
出力:
  report_axislaw_A_shortest.html

再生成:
  python3 case/42.isobutane_wt/report_axislaw_A_shortest.py
"""
from __future__ import annotations

import html
import json
from collections import Counter
from pathlib import Path


CASE = Path(__file__).resolve().parent
INPUT = CASE / "optimize_axislaw_A_shortest.json"
OUTPUT = CASE / "report_axislaw_A_shortest.html"
R_T_M = 0.05375


def _line_chart(
    series: list[tuple[str, str, list[tuple[float, float]]]],
    *,
    title: str,
    x_label: str,
    y_label: str,
    threshold: float | None = None,
    width: int = 760,
    height: int = 330,
) -> str:
    """依存ライブラリなしの小さなインラインSVG折れ線図。"""
    left, right, top, bottom = 70, 22, 42, 55
    plot_w, plot_h = width - left - right, height - top - bottom
    points = [p for _, _, pp in series for p in pp]
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    if threshold is not None:
        ys.append(threshold)
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    xpad = max((x1 - x0) * 0.04, 1e-9)
    ypad = max((y1 - y0) * 0.16, 1e-6)
    x0, x1 = x0 - xpad, x1 + xpad
    y0, y1 = y0 - ypad, y1 + ypad

    def sx(x: float) -> float:
        return left + (x - x0) / (x1 - x0) * plot_w

    def sy(y: float) -> float:
        return top + (y1 - y) / (y1 - y0) * plot_h

    out = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title)}">',
        f'<text x="{left}" y="24" class="chart-title">{html.escape(title)}</text>',
    ]
    for i in range(6):
        y = y0 + (y1 - y0) * i / 5
        py = sy(y)
        out.append(f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" class="grid"/>')
        out.append(f'<text x="{left-9}" y="{py+4:.1f}" text-anchor="end" class="tick">{y:.4f}</text>')
    for i in range(6):
        x = x0 + (x1 - x0) * i / 5
        px = sx(x)
        out.append(f'<text x="{px:.1f}" y="{height-bottom+24}" text-anchor="middle" class="tick">{x:.3f}</text>')
    if threshold is not None:
        py = sy(threshold)
        out.append(f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" class="threshold"/>')
        out.append(f'<text x="{width-right-4}" y="{py-7:.1f}" text-anchor="end" class="threshold-label">基準 {threshold:g}</text>')
    for label, color, pp in series:
        coords = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in pp)
        out.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2.7"/>')
        for x, y in pp:
            out.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.7" fill="{color}"/>')
    legend_x = left
    for label, color, _ in series:
        out.append(f'<circle cx="{legend_x+5}" cy="{height-11}" r="4" fill="{color}"/>')
        out.append(f'<text x="{legend_x+15}" y="{height-7}" class="legend">{html.escape(label)}</text>')
        legend_x += 155
    out.append(f'<text x="{left+plot_w/2:.1f}" y="{height-bottom+46}" text-anchor="middle" class="axis-label">{html.escape(x_label)}</text>')
    out.append(f'<text transform="translate(17 {top+plot_h/2:.1f}) rotate(-90)" text-anchor="middle" class="axis-label">{html.escape(y_label)}</text>')
    out.append("</svg>")
    return "".join(out)


def _fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def main() -> None:
    data = json.loads(INPUT.read_text())
    records = [r for r in data["records"] if "error" not in r]
    winner = data["winner"]
    win_rec = next(
        r for r in records
        if r.get("search_stage") == "final_MK"
        and r["L_c"] == winner["L_c"] and r["M_K"] == winner["M_K"]
    )

    micro = [r for r in records if r.get("search_stage") == "micro"]
    probes = sorted(
        (r for r in records if r.get("search_stage") == "boundary_probe"),
        key=lambda r: r["L_c"],
    )
    final_mk = sorted(
        (r for r in records if r.get("search_stage") == "final_MK"),
        key=lambda r: r["M_K"],
    )
    stage_counts = Counter(r.get("search_stage") for r in records)

    micro_best = []
    for lc in sorted({r["L_c"] for r in micro}):
        rr = [r for r in micro if r["L_c"] == lc and r["hard_gate_pass"]
              and r["axis_gates"]["mpp_quality_ok"]
              and r["topology"]["min_sin_angle_interior"] >= 0.02]
        if rr:
            best = max(rr, key=lambda r: r["margin"]["min_mu_minus_theta_deg"])
            micro_best.append((lc, best["margin"]["min_mu_minus_theta_deg"]))

    margin_chart = _line_chart(
        [
            ("600点・各長さの最大余裕", "#74829a", micro_best),
            ("1200点・境界確認", "#e05252", [
                (r["L_c"], r["margin"]["min_mu_minus_theta_deg"]) for r in probes
            ]),
        ],
        title="最短境界の解像度確認",
        x_label="L_c / r_t",
        y_label="min(μ−θ) [deg]",
        threshold=1.0,
    )
    mk_margin_chart = _line_chart(
        [("1200点", "#2a7f9e", [
            (r["M_K"], r["margin"]["min_mu_minus_theta_deg"]) for r in final_mk
        ])],
        title="L_c=39.05 における M_K 感度",
        x_label="M_K",
        y_label="min(μ−θ) [deg]",
        threshold=1.0,
    )
    mk_theta_chart = _line_chart(
        [("1200点", "#7a5bb6", [(r["M_K"], r["theta_max_deg"]) for r in final_mk])],
        title="壁最大角の M_K 感度",
        x_label="M_K",
        y_label="θ_max [deg]",
    )

    boundary_rows = "".join(
        f"<tr class=\"{'pass' if r['margin']['min_mu_minus_theta_deg'] >= 1 else 'fail'}\">"
        f"<td>{r['L_c']:.2f}</td><td>{r['M_K']:.2f}</td>"
        f"<td>{r['margin']['min_mu_minus_theta_deg']:.6f}°</td>"
        f"<td>{r['topology']['min_sin_angle_interior']:.5f}</td>"
        f"<td>{'PASS' if r['margin']['min_mu_minus_theta_deg'] >= 1 else 'FAIL'}</td></tr>"
        for r in probes[:8]
    )
    final_rows = "".join(
        f"<tr class=\"{'pass' if r['margin']['min_mu_minus_theta_deg'] >= 1 else 'fail'}\">"
        f"<td>{r['M_K']:.2f}</td><td>{r['theta_max_deg']:.5f}°</td>"
        f"<td>{r['margin']['min_mu_minus_theta_deg']:.6f}°</td>"
        f"<td>{r['topology']['min_sin_angle_interior']:.5f}</td>"
        f"<td>{'PASS' if r['margin']['min_mu_minus_theta_deg'] >= 1 else 'FAIL'}</td></tr>"
        for r in final_mk if 2.73 <= r["M_K"] <= 2.81
    )

    lc_m = winner["L_c"] * R_T_M
    xf_m = winner["x_F"] * R_T_M
    full_m = (winner["x_F"] + 12.5) * R_T_M
    generated = f"""<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>軸則 A — 最短ロバスト軸長探索</title>
<style>
:root{{--ink:#17202b;--muted:#687385;--line:#dce2ea;--paper:#f4f6f9;--card:#fff;--navy:#19344d;--blue:#2a7f9e;--green:#25765a;--red:#c74343;--amber:#a66b16;--violet:#7253ad}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--paper);color:var(--ink);font-family:Inter,"Noto Sans JP","Yu Gothic",sans-serif;line-height:1.72}}
.top{{background:linear-gradient(125deg,#10283c,#214c68 62%,#2d6c7e);color:#fff;padding:52px 24px 46px}}
.wrap{{max-width:1120px;margin:auto}} .eyebrow{{font-size:12px;letter-spacing:.16em;text-transform:uppercase;opacity:.75;font-weight:700}}
h1{{font-size:clamp(30px,5vw,53px);line-height:1.18;margin:10px 0 14px}} .lead{{max-width:850px;font-size:18px;color:#dbe8ee}}
.status{{display:flex;gap:10px;flex-wrap:wrap;margin-top:22px}} .badge{{border-radius:999px;padding:6px 12px;font-size:12px;font-weight:800;letter-spacing:.05em}}
.ok{{background:#dff6eb;color:#145f43}} .pending{{background:#fff0cf;color:#80530d}} main{{padding:28px 20px 70px}}
.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-top:-48px;position:relative}}
.metric,.card{{background:var(--card);border:1px solid var(--line);border-radius:15px;box-shadow:0 5px 20px rgba(25,52,77,.06)}}
.metric{{padding:20px}} .metric .label{{font-size:12px;color:var(--muted);font-weight:700}} .metric .value{{font-size:26px;font-weight:800;color:var(--navy);margin-top:4px}}
.metric .sub{{font-size:12px;color:var(--muted)}} section{{margin-top:34px}} h2{{font-size:24px;color:var(--navy);margin:0 0 14px}} h3{{font-size:16px;margin:0 0 8px;color:var(--navy)}}
.card{{padding:24px}} .decision{{border-left:5px solid var(--green)}} .warning{{border-left:5px solid var(--amber);background:#fffdf7}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:16px}} .grid3{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}}
.mini{{padding:18px;border:1px solid var(--line);border-radius:12px;background:#fbfcfe}} .mini b{{display:block;font-size:18px;color:var(--navy)}}
.formula{{background:#112638;color:#e8f0f4;border-radius:12px;padding:19px 22px;font-family:"SFMono-Regular",Consolas,monospace;overflow:auto}}
.flow{{display:flex;align-items:stretch;gap:9px;overflow:auto}} .step{{min-width:150px;flex:1;background:#eef3f7;border-radius:11px;padding:14px;text-align:center}}
.step b{{font-size:22px;color:var(--navy);display:block}} .arrow{{display:flex;align-items:center;color:#91a0af;font-size:22px}}
figure{{margin:0}} svg{{width:100%;height:auto;background:#fff;border-radius:12px}} .chart-title{{font:700 15px sans-serif;fill:#19344d}} .grid{{stroke:#e6ebf1;stroke-width:1}} .tick,.legend,.axis-label{{font:11px sans-serif;fill:#687385}} .threshold{{stroke:#c74343;stroke-width:1.5;stroke-dasharray:6 5}} .threshold-label{{font:700 11px sans-serif;fill:#c74343}}
table{{width:100%;border-collapse:collapse;font-size:13px}} th{{text-align:left;background:#edf2f6;color:#34495e}} th,td{{padding:9px 10px;border-bottom:1px solid var(--line)}} tr.pass td:last-child{{color:var(--green);font-weight:800}} tr.fail td:last-child{{color:var(--red);font-weight:800}}
.callout{{font-size:20px;font-weight:800;color:var(--navy)}} .muted{{color:var(--muted)}} code{{background:#edf1f5;padding:2px 5px;border-radius:5px}} ul{{padding-left:20px}}
.footer{{border-top:1px solid var(--line);margin-top:45px;padding-top:18px;color:var(--muted);font-size:12px}}
@media(max-width:800px){{.metrics,.grid2,.grid3{{grid-template-columns:1fr 1fr}}}} @media(max-width:540px){{.metrics,.grid2,.grid3{{grid-template-columns:1fr}} .metric{{padding:16px}}}}
</style>
</head>
<body>
<header class="top"><div class="wrap">
<div class="eyebrow">Axis-Mach nozzle design · M6 / R=3 · 2026-08-17</div>
<h1>軸則 A — 最短ロバスト軸長探索</h1>
<p class="lead">壁角を最小にするのではなく、品質ゲートと工学的余裕を満たす範囲で、閉包後端 <i>x</i><sub>F</sub> を最短化した。</p>
<div class="status"><span class="badge ok">DESIGN GATE PASS</span><span class="badge pending">CFD PENDING</span><span class="badge pending">DESIGN LIMIT — MARGIN ACTIVE</span></div>
</div></header>
<main><div class="wrap">
<div class="metrics">
  <div class="metric"><div class="label">最短軸長 L<sub>c</sub></div><div class="value">39.05 r<sub>t</sub></div><div class="sub">{lc_m:.3f} m</div></div>
  <div class="metric"><div class="label">Knot Mach</div><div class="value">M<sub>K</sub> = 2.76</div><div class="sub">合格帯 2.76–2.78</div></div>
  <div class="metric"><div class="label">閉包後端 x<sub>F</sub></div><div class="value">94.349 r<sub>t</sub></div><div class="sub">スロート基準 {xf_m:.3f} m</div></div>
  <div class="metric"><div class="label">支配余裕 min(μ−θ)</div><div class="value">1.00095°</div><div class="sub">基準 1.00000°</div></div>
</div>

<section><div class="card decision">
<h2>結論</h2><div class="callout">設計上の最短境界は L<sub>c</sub>=39.05, M<sub>K</sub>=2.76。</div>
<p>1200点で hard gate、単峰条件、μ−θ余裕、特性線トポロジ余裕をすべて満たす最短候補である。入口上流長を含む概算端間長は <b>{full_m:.3f} m</b>。ただし余裕制約が能動で、CFD未実施のため、生産既定変更ではなく「設計境界」として採択する。</p>
</div></section>

<section><h2>1. 何をもって「最良」としたか</h2>
<div class="grid2"><div class="card"><h3>必須制約</h3><ul>
<li>軸Mach単調、M″反転1回以下</li><li>壁QA・CFD壁構築ゲート合格</li><li>内部flip 0、非隣接交差 0、内部点堆積 0</li><li>min(μ−θ) ≥ 1°、最小sin角 ≥ 0.02</li>
</ul></div><div class="card"><h3>辞書順の目的</h3><ol>
<li>合格集合で x<sub>F</sub> を最小化</li><li>同じ長さなら θ<sub>max</sub> を最小化</li><li>さらに同等なら余裕・曲率を比較</li></ol>
<p class="muted">1°と0.02は数学的限界ではなく工学的余裕。</p></div></div></section>

<section><h2>2. A の数学的構成</h2><div class="card">
<p>Aはknot前後の2区間。第1区間はA端のHallアンカーとK端の値・勾配・M″=0を結ぶ5次Hermite、第2区間は終端でM′=M″=0になる4次式。</p>
<div class="formula">区分1 [x_A, x_K] : quintic Hermite<br>
区分2 [x_K, x_E] : M = M_K + ΔM₂ (2u − 2u³ + u⁴)<br>
s_K = 2ΔM₂/L₂, &nbsp; L₁ = 2ΔM₁/(M′_A + s_K)<br>
DOF = (L_c, M_K), &nbsp; x_E = x_A + L_c</div>
<div class="grid3" style="margin-top:14px">
<div class="mini"><span class="muted">knot位置</span><b>x<sub>K</sub> = {_fmt(win_rec['x_K'],4)} r<sub>t</sub></b>ξ<sub>K</sub> = {_fmt(win_rec['xi_K'],5)}</div>
<div class="mini"><span class="muted">区分長</span><b>L₁ / L₂</b>{_fmt(win_rec['L_1'],4)} / {_fmt(win_rec['L_2'],4)}</div>
<div class="mini"><span class="muted">knot勾配</span><b>s<sub>K</sub> = {_fmt(win_rec['s_K'],5)}</b>M″反転 = 1</div>
</div></div></section>

<section><h2>3. 探索の絞り込み</h2><div class="card">
<div class="flow">
<div class="step"><span>粗探索</span><b>{stage_counts['coarse']}</b>600点</div><div class="arrow">→</div>
<div class="step"><span>細探索</span><b>{stage_counts['fine']}</b>600点</div><div class="arrow">→</div>
<div class="step"><span>0.01格子</span><b>{stage_counts['micro']}</b>600点</div><div class="arrow">→</div>
<div class="step"><span>境界確認</span><b>{stage_counts['boundary_probe']}</b>1200点</div><div class="arrow">→</div>
<div class="step"><span>M<sub>K</sub>確定</span><b>{stage_counts['final_MK']}</b>1200点</div>
</div><p class="muted">全探索87.76秒、最大RSS 308 MiB。2400点は使用していない。</p>
</div></section>

<section><h2>4. 最短境界 — 600点では39.03を誤採択</h2>
<div class="grid2"><figure class="card">{margin_chart}</figure><div class="card"><table>
<thead><tr><th>L<sub>c</sub></th><th>M<sub>K</sub></th><th>μ−θ</th><th>min sin</th><th>判定</th></tr></thead><tbody>{boundary_rows}</tbody></table>
<p class="muted">600点の39.03は1.00027°だったが、1200点で0.99748°へ低下。39.05が最初の合格点になった。</p></div></div></section>

<section><h2>5. M<sub>K</sub>感度 — 合格帯は狭い</h2>
<div class="grid2"><figure class="card">{mk_margin_chart}</figure><figure class="card">{mk_theta_chart}</figure></div>
<div class="card" style="margin-top:16px"><table><thead><tr><th>M<sub>K</sub></th><th>θ<sub>max</sub></th><th>μ−θ</th><th>min sin</th><th>判定</th></tr></thead><tbody>{final_rows}</tbody></table>
<p class="muted">合格は2.76–2.78。二次目的 θ<sub>max</sub> が最小の2.76を選択。</p></div></section>

<section><h2>6. Characteristic topology と壁品質</h2><div class="grid3">
<div class="mini"><span class="muted">内部flip / 交差</span><b>0 / 0</b>hard gate PASS</div>
<div class="mini"><span class="muted">最小sin角</span><b>{win_rec['topology']['min_sin_angle_interior']:.5f}</b>基準0.02の2.09倍</div>
<div class="mini"><span class="muted">最小signed area</span><b>{win_rec['topology']['min_signed_area_rel_interior']:.4f}</b>内部wrong sign 0</div>
<div class="mini"><span class="muted">θ<sub>max</sub></span><b>{win_rec['theta_max_deg']:.4f}°</b>出口角 ≈ 0°</div>
<div class="mini"><span class="muted">κ₀R</span><b>{win_rec['kappa0_R_ratio']:.4f}</b>Hallアンカーと整合</div>
<div class="mini"><span class="muted">出口半径誤差</span><b>{100*win_rec['r_F_err_rel']:.4f}%</b>1D予測比</div>
</div></section>

<section><h2>7. 採否と残る不確実性</h2><div class="grid2">
<div class="card decision"><h3>採択</h3><p><b>L<sub>c</sub>=39.05, M<sub>K</sub>=2.76</b>を設計上の最短境界として記録する。従来の40より軸Mach指定長を0.95 r<sub>t</sub>（約51 mm）短縮。</p></div>
<div class="card warning"><h3>生産採用は保留</h3><ul><li>μ−θ余裕は閾値の0.00095°上で、制約が能動</li><li>CFD性能と粘性補正後の余裕は未検証</li><li>1°という閾値自体は工学的選択</li><li>0.01 r<sub>t</sub>格子より細かい連続最適値は主張しない</li></ul></div>
</div></section>

<section><h2>8. 再現性</h2><div class="card"><ul>
<li>探索: <code>optimize_axislaw_A_shortest.py</code></li>
<li>全記録: <code>optimize_axislaw_A_shortest.json</code></li>
<li>表: <code>optimize_axislaw_A_shortest_summary.csv</code></li>
<li>解像度: 探索600、最終1200、2400未使用</li>
<li>回帰: moc_diagnostics / axislaw / onepoint / bspline / inverse / axismach_wall — ALL PASS</li>
</ul></div></section>

<div class="footer">Forge nozzle design report · generated from optimize_axislaw_A_shortest.json · design-only, no CFD run created.</div>
</div></main></body></html>"""
    OUTPUT.write_text(generated)
    print(f"saved {OUTPUT}")


if __name__ == "__main__":
    main()
