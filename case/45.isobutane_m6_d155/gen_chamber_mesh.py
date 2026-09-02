"""チャンバー付きノズルメッシュ生成 (出口コーナー不安定の A/B 用, Codex レビュー準拠)。

3 ブロック transfinite (gmsh):
  N : ノズル内部 (壁スプライン=run_0011 の CPG 物理壁, 軸まで)
  C1: 出口後方の噴流域 (x_lip→x_end, r 0→r_lip)  — N と出口面で共形
  C2: 外周 (x_lip→x_end, r_lip→R_out)。x=x_lip の縦面 (r_lip→R_out) はフランジ壁
チャンバー: 長さ 4 D_e・外径 3 D_e。BL は粗 (wall_first_frac 5e-4 = y+~50 相当,
coarse メッシュで不安定が再現することを確認済み [run_0015_cflsweep/coarse_cfl8])。

usage: design/.venv-opt/bin/python gen_chamber_mesh.py <outdir>
出力: <outdir>/chamber.geo, chamber.msh (gmsh -2 実行込み)
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else CASE / "mesh_chamber"
OUT.mkdir(parents=True, exist_ok=True)

# 完全壁輪郭 (入口管+収縮+スロート+拡大)。wall_physical.csv はスロート以降のみなので
# 実メッシュの wall bcond ノードから抽出した wall_full.csv を使う (罠: 2026-09-02)。
wp = np.loadtxt(CASE / "mesh_chamber" / "wall_full.csv", delimiter=",", skiprows=1)
x_lip, r_lip = wp[-1, 0], wp[-1, 1]
D_e = 2 * r_lip
X_END = x_lip + 4.0 * D_e
R_OUT = 3.0 * D_e
NI_CONV, NI_DIV, NJ = 180, 750, 65   # ノズル: スロートで分割し両側からスロートへ細分
NX_C, NJ2 = 110, 60           # チャンバー
WALL_FRAC = 5.0e-4            # 半径方向第一セル割合 (coarse 相当)


def prog_for_first_frac(frac, n_intervals):
    """transfinite Progression q: 第一区間割合 = (q-1)/(q^N-1) = frac を解く。
    frac が一様割合以上なら細分不要 → q=1 (一様)。"""
    if frac >= 1.0 / n_intervals:
        return 1.0
    from scipy.optimize import brentq
    f = lambda q: (q - 1.0) / (q ** n_intervals - 1.0) - frac
    return brentq(f, 1.0000001, 2.0)


q_rad = prog_for_first_frac(WALL_FRAC, NJ - 1)
# スロート位置 (最小半径)
it = int(np.argmin(wp[:, 1]))
x_t = wp[it, 0]
DX_T = 0.003   # スロート近傍の目標 Δx [m] (coarse 実績相当)
L_conv = x_t - wp[0, 0]
L_div = x_lip - x_t
q_conv = prog_for_first_frac(DX_T / L_conv, NI_CONV - 1)   # スロート側が細 (向き反転で使う)
q_div  = prog_for_first_frac(DX_T / L_div, NI_DIV - 1)
# チャンバー x: リップ側第一区間 ~ ノズル出口 Δx に合わせる
dx_noz = L_div * (q_div - 1.0) / (q_div ** (NI_DIV - 1) - 1.0) * q_div ** (NI_DIV - 2)
q_cx = prog_for_first_frac(min(dx_noz, 0.1) / (X_END - x_lip), NX_C - 1)
# 外周 r: リップ側第一区間 ~ 2 mm
q_c2 = prog_for_first_frac(0.002 / (R_OUT - r_lip), NJ2 - 1)

# 壁スプライン点 (間引き, スロート点は必ず含める)
sk = max(1, len(wp) // 350)
keep = sorted(set(list(range(0, len(wp), sk)) + [it, len(wp) - 1]))
pts = wp[keep]
i_throat_pt = keep.index(it)

g = []
g.append("Mesh.MshFileVersion = 4.1;")
g.append("lc = 1.0;")
pid = 1
for xw, rw in pts:
    g.append(f"Point({pid}) = {{{xw:.9f}, {rw:.9f}, 0, lc}};")
    pid += 1
n_wall = pid - 1
P_IN_W, P_LIP = 1, n_wall            # 壁: 入口端点 → リップ
x_in = pts[0, 0]
P_IN_A = pid; g.append(f"Point({P_IN_A}) = {{{x_in:.9f}, 0, 0, lc}};"); pid += 1   # 入口軸
P_EXIT_A = pid; g.append(f"Point({P_EXIT_A}) = {{{x_lip:.9f}, 0, 0, lc}};"); pid += 1  # 出口軸
P_END_A = pid; g.append(f"Point({P_END_A}) = {{{X_END:.9f}, 0, 0, lc}};"); pid += 1
P_END_LIP = pid; g.append(f"Point({P_END_LIP}) = {{{X_END:.9f}, {r_lip:.9f}, 0, lc}};"); pid += 1
P_END_OUT = pid; g.append(f"Point({P_END_OUT}) = {{{X_END:.9f}, {R_OUT:.9f}, 0, lc}};"); pid += 1
P_LIP_OUT = pid; g.append(f"Point({P_LIP_OUT}) = {{{x_lip:.9f}, {R_OUT:.9f}, 0, lc}};"); pid += 1

P_T_W = 1 + i_throat_pt                               # スロート壁点
P_T_A = pid; g.append(f"Point({P_T_A}) = {{{x_t:.9f}, 0, 0, lc}};"); pid += 1
g.append(f"Spline(1) = {{{', '.join(str(i) for i in range(1, P_T_W + 1))}}};")           # 収縮壁
g.append(f"Spline(11) = {{{', '.join(str(i) for i in range(P_T_W, n_wall + 1))}}};")     # 拡大壁
g.append(f"Line(2) = {{{P_IN_W}, {P_IN_A}}};")        # 入口 (壁→軸)
g.append(f"Line(3) = {{{P_IN_A}, {P_T_A}}};")         # 軸 (収縮側)
g.append(f"Line(12) = {{{P_T_A}, {P_EXIT_A}}};")      # 軸 (拡大側)
g.append(f"Line(13) = {{{P_T_A}, {P_T_W}}};")         # スロート面 (軸→壁)
g.append(f"Line(4) = {{{P_EXIT_A}, {P_LIP}}};")       # 出口面 (軸→リップ) = N/C1 共有
g.append(f"Line(5) = {{{P_EXIT_A}, {P_END_A}}};")     # チャンバー軸
g.append(f"Line(6) = {{{P_END_A}, {P_END_LIP}}};")    # 末端面 下半
g.append(f"Line(7) = {{{P_LIP}, {P_END_LIP}}};")      # r=r_lip 水平線 = C1/C2 共有
g.append(f"Line(8) = {{{P_END_LIP}, {P_END_OUT}}};")  # 末端面 上半
g.append(f"Line(9) = {{{P_LIP_OUT}, {P_END_OUT}}};")  # 外周 r=R_out
g.append(f"Line(10) = {{{P_LIP}, {P_LIP_OUT}}};")     # フランジ (x=x_lip, r_lip→R_out)

g.append("Curve Loop(1) = {2, 3, 13, -1};  Plane Surface(1) = {1};")     # N1 (収縮)
g.append("Curve Loop(4) = {12, 4, -11, -13}; Plane Surface(4) = {4};")   # N2 (拡大)
g.append("Curve Loop(2) = {5, 6, -7, -4}; Plane Surface(2) = {2};")     # C1
g.append("Curve Loop(3) = {7, 8, -9, -10}; Plane Surface(3) = {3};")    # C2

g.append(f"Transfinite Curve {{-1, -3}} = {NI_CONV} Using Progression {q_conv:.7f};")  # スロートへ細分
g.append(f"Transfinite Curve {{11, 12}} = {NI_DIV} Using Progression {q_div:.7f};")
g.append(f"Transfinite Curve {{-2, 4, 6, 13}} = {NJ} Using Progression {1.0/q_rad:.7f};")  # 壁側細分
if q_cx > 1.0:
    g.append(f"Transfinite Curve {{5, 7, 9}} = {NX_C} Using Progression {q_cx:.7f};")
else:
    g.append(f"Transfinite Curve {{5, 7, 9}} = {NX_C} Using Progression 1.0;")
g.append(f"Transfinite Curve {{10, 8}} = {NJ2} Using Progression {q_c2:.7f};")
g.append("Transfinite Surface {1} = {%d, %d, %d, %d};" % (P_IN_A, P_T_A, P_T_W, P_IN_W))
g.append("Transfinite Surface {4} = {%d, %d, %d, %d};" % (P_T_A, P_EXIT_A, P_LIP, P_T_W))
g.append("Transfinite Surface {2} = {%d, %d, %d, %d};" % (P_EXIT_A, P_END_A, P_END_LIP, P_LIP))
g.append("Transfinite Surface {3} = {%d, %d, %d, %d};" % (P_LIP, P_END_LIP, P_END_OUT, P_LIP_OUT))
g.append("Recombine Surface {1, 2, 3, 4};")

g.append('Physical Curve("inlet", 1) = {2};')
g.append('Physical Curve("outlet", 2) = {6, 8, 9};')
g.append('Physical Curve("wall", 3) = {1, 11, 10};')     # ノズル壁 + フランジ
g.append('Physical Curve("axis", 4) = {3, 12, 5};')
g.append('Physical Surface("fluid", 5) = {1, 2, 3, 4};')

(OUT / "chamber.geo").write_text("\n".join(g) + "\n")
print(f"geo: throat=({x_t:.4f}) lip=({x_lip:.3f},{r_lip:.3f}) X_END={X_END:.2f} R_OUT={R_OUT:.2f} "
      f"q_rad={q_rad:.4f} q_cx={q_cx:.4f} q_c2={q_c2:.4f}")
r = subprocess.run(["gmsh", "-2", str(OUT / "chamber.geo"), "-o", str(OUT / "chamber.msh"),
                    "-format", "msh41"], capture_output=True, text=True)
print(r.stdout[-400:] if r.returncode == 0 else r.stderr[-800:])
print("rc", r.returncode)
