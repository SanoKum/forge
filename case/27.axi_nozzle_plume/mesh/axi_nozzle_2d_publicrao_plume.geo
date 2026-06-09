// Rao nozzle — axisymmetric 2D mesh with far-field plume extension
//
// outer領域をx方向・r方向それぞれ2分割した構造格子版
//
// Domain layout (mm; ScalingFactor = 0.001):
//
//   D_exit  = 2 * r_exit = 65.479 mm
//   r_far   = r_exit + 5*D_exit = 360.135 mm
//   x_far   = div_end_x + 10*D_exit = 722.682 mm
//   r_mid   = r_exit + frac_r * (r_far - r_exit)
//   x_mid   = div_end_x + frac_x * (x_far - div_end_x)
//
//  x=div_end_x   x=x_mid        x=x_far
//  r=r_far  p_outer1_top ─── p_outer2_mid ─── p_exit_top     ← outer2_1 / outer2_2
//            │     [E1]       │     [E2]        │
//  r=r_mid  p_outer1_mid ─── p_inner_corner ── p_exit_r_mid  ← l_r_mid_1 / l_r_mid_2
//            │     [D1]       │     [D2]        │
//  r=r_exit p_wall_div_end ─ p_interface_mid ─ p_exit_mid    ← interface_1 / interface_2
//            │     [C1]       │     [C2]        │
//  r=0      p_axis_div_end ─ p_x_mid ─────────  p_axis_exit  ← axis_far_1 / axis_far_2
//
//           outer1_bot/top   mid_v_c/d/e      outer3_bot/lo/hi
//
// Physical boundaries:
//   inlet  (1) — nozzle inlet
//   wall   (3) — nozzle wall (conv + div)
//   axis   (4) — symmetry axis (full length)
//   outer1 (5) — radial face at nozzle exit (r_exit → r_far)
//   outer2 (6) — outer boundary r=r_far (downstream)
//   outer3 (7) — exit plane x=x_far (r=0 → r_far)

Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers   = 0;
Geometry.LineNumbers    = 0;
Mesh.ScalingFactor      = 0.001;

lc = 2.0;

// ---- nozzle resolution ----
nx_conv = 80;
nx_div  = 80;
ny      = 60;


prog_r_noz = 0.93;

// ---- plume x-direction: 2-split (adjust nx_p1/nx_p2 and prog_x1/prog_x2 independently) ----
nx_p1   = 60;     // cells in x, 左半 (nozzle exit side)
nx_p2   = 40;     // cells in x, 右半 (exit plane side)
prog_x1   = 1.05;   // x progression, 左半
prog_x1_2 = 1.035;   // x progression, 左半
prog_x1_3 = 1.00;   // x progression, 左半
prog_x2 = 1.03;   // x progression, 右半 (downstream へ拡大)

// ---- outer r-direction: 2-split ----
//   lower (r_exit → r_mid): ny_out_bot, prog_r_bot (fine near r_exit)
//   upper (r_mid  → r_far): ny_out_top, prog_r_top (fine near r_mid)
ny_out_bot = 40;
ny_out_top = 20;
prog_r_bot = 1.18;
prog_r_bot2 = 1.09;
prog_r_bot3 = 1.06;
prog_r_top = 1.08;

// ---- split fractions ----
frac_r = 0.3;    // r_mid = r_exit + frac_r * (r_far - r_exit)
frac_x = 0.3;    // x_mid = div_end_x + frac_x * (x_far - div_end_x)

Include "axi_nozzle_points_public_rao.geo";

// Derived dimensions
D_exit = 2.0 * r_exit;
r_far  = r_exit + 5.0 * D_exit;
x_far  = div_end_x + 10.0 * D_exit;
r_mid  = r_exit + frac_r * (r_far - r_exit);
x_mid  = div_end_x + frac_x * (x_far - div_end_x);

// ==============================================================
// Nozzle wall splines
// ==============================================================
l_conv_top = newl; Spline(l_conv_top) = conv_pts[];
l_div_top  = newl; Spline(l_div_top)  = div_pts[];
Transfinite Curve {l_conv_top} = nx_conv + 1;
Transfinite Curve {l_div_top}  = nx_div  + 1;

// Axis points for nozzle region
p_axis_inlet   = newp; Point(p_axis_inlet)   = {xmin,      0.0, 0.0, lc};
p_axis_throat  = newp; Point(p_axis_throat)  = {throat_x,  0.0, 0.0, lc};
p_axis_div_end = newp; Point(p_axis_div_end) = {div_end_x, 0.0, 0.0, lc};

l_conv_bot = newl; Line(l_conv_bot) = {p_axis_inlet,  p_axis_throat};
l_div_bot  = newl; Line(l_div_bot)  = {p_axis_throat, p_axis_div_end};
Transfinite Curve {l_conv_bot} = nx_conv + 1;
Transfinite Curve {l_div_bot}  = nx_div  + 1;

// Radial dividers (inlet / throat / nozzle-exit)
l_inlet    = newl; Line(l_inlet)    = {p_axis_inlet,  conv_pts[0]};
n_conv     = #conv_pts[];
p_wall_throat  = conv_pts[n_conv - 1];
l_throat_v = newl; Line(l_throat_v) = {p_axis_throat, p_wall_throat};

n_div          = #div_pts[];
p_wall_div_end = div_pts[n_div - 1];
l_divend_v = newl; Line(l_divend_v) = {p_axis_div_end, p_wall_div_end};

Transfinite Curve {l_inlet}    = ny + 1 Using Progression prog_r_noz;
Transfinite Curve {l_throat_v} = ny + 1 Using Progression prog_r_noz;
Transfinite Curve {l_divend_v} = ny + 1 Using Progression prog_r_noz;

// Convergent surface
ll_conv = newll; Curve Loop(ll_conv) = {l_conv_bot, l_throat_v, -l_conv_top, -l_inlet};
s_conv  = news;  Plane Surface(s_conv) = {ll_conv};
Transfinite Surface {s_conv}; Recombine Surface {s_conv};

// Divergent surface
ll_div = newll; Curve Loop(ll_div) = {l_div_bot, l_divend_v, -l_div_top, -l_throat_v};
s_div  = news;  Plane Surface(s_div) = {ll_div};
Transfinite Surface {s_div}; Recombine Surface {s_div};

// ==============================================================
// Far-field plume extension — 6 blocks (C1/C2/D1/D2/E1/E2)
// ==============================================================

// ---- corner / boundary points ----
p_outer1_top    = newp; Point(p_outer1_top)    = {div_end_x, r_far,   0.0, lc};
p_outer2_mid    = newp; Point(p_outer2_mid)    = {x_mid,     r_far,   0.0, lc};
p_exit_top      = newp; Point(p_exit_top)      = {x_far,     r_far,   0.0, lc};
p_exit_mid      = newp; Point(p_exit_mid)      = {x_far,     r_exit,  0.0, lc};
p_axis_exit     = newp; Point(p_axis_exit)     = {x_far,     0.0,     0.0, lc};

// ---- midpoint points (split handles) ----
p_outer1_mid    = newp; Point(p_outer1_mid)    = {div_end_x, r_mid,   0.0, lc};
p_x_mid         = newp; Point(p_x_mid)         = {x_mid,     0.0,     0.0, lc};
p_interface_mid = newp; Point(p_interface_mid) = {x_mid,     r_exit,  0.0, lc};
p_inner_corner  = newp; Point(p_inner_corner)  = {x_mid,     r_mid,   0.0, lc};
p_exit_r_mid    = newp; Point(p_exit_r_mid)    = {x_far,     r_mid,   0.0, lc};

// ---- outer boundary curves ----
// outer1: x=div_end_x, r_exit→r_mid→r_far
l_outer1_bot = newl; Line(l_outer1_bot) = {p_wall_div_end, p_outer1_mid};
l_outer1_top = newl; Line(l_outer1_top) = {p_outer1_mid,   p_outer1_top};

// outer2: r=r_far, div_end_x→x_mid→x_far
l_outer2_1 = newl; Line(l_outer2_1) = {p_outer1_top, p_outer2_mid};
l_outer2_2 = newl; Line(l_outer2_2) = {p_outer2_mid, p_exit_top};

// outer3: x=x_far, r_far→r_mid→r_exit→0
l_outer3_hi  = newl; Line(l_outer3_hi)  = {p_exit_top,   p_exit_r_mid};
l_outer3_lo  = newl; Line(l_outer3_lo)  = {p_exit_r_mid, p_exit_mid};
l_outer3_bot = newl; Line(l_outer3_bot) = {p_exit_mid,   p_axis_exit};

// ---- internal interfaces ----
// extended axis: r=0, div_end_x→x_mid→x_far
l_axis_far_1 = newl; Line(l_axis_far_1) = {p_axis_div_end, p_x_mid};
l_axis_far_2 = newl; Line(l_axis_far_2) = {p_x_mid,        p_axis_exit};

// r=r_exit interface: div_end_x→x_mid→x_far
l_interface_1 = newl; Line(l_interface_1) = {p_wall_div_end,  p_interface_mid};
l_interface_2 = newl; Line(l_interface_2) = {p_interface_mid, p_exit_mid};

// x=x_mid 縦仕切り (r: 0→r_exit→r_mid→r_far)
l_mid_v_c = newl; Line(l_mid_v_c) = {p_x_mid,        p_interface_mid};
l_mid_v_d = newl; Line(l_mid_v_d) = {p_interface_mid, p_inner_corner};
l_mid_v_e = newl; Line(l_mid_v_e) = {p_inner_corner,  p_outer2_mid};

// r=r_mid 横仕切り (x: div_end_x→x_mid→x_far)
l_r_mid_1 = newl; Line(l_r_mid_1) = {p_outer1_mid,  p_inner_corner};
l_r_mid_2 = newl; Line(l_r_mid_2) = {p_inner_corner, p_exit_r_mid};

// ---- transfinite counts ----
// x方向 (outer左半)
Transfinite Curve {l_axis_far_1}  = nx_p1 + 1 Using Progression prog_x1;
Transfinite Curve {l_interface_1} = nx_p1 + 1 Using Progression prog_x1;
Transfinite Curve {l_outer2_1}    = nx_p1 + 1 Using Progression prog_x1_3;
Transfinite Curve {l_r_mid_1}     = nx_p1 + 1 Using Progression prog_x1_2;

// x方向 (outer右半)
Transfinite Curve {l_axis_far_2}  = nx_p2 + 1 Using Progression prog_x2;
Transfinite Curve {l_interface_2} = nx_p2 + 1 Using Progression prog_x2;
Transfinite Curve {l_outer2_2}    = nx_p2 + 1 Using Progression prog_x2;
Transfinite Curve {l_r_mid_2}     = nx_p2 + 1 Using Progression prog_x2;

// r方向 (inner, r=0→r_exit)
Transfinite Curve {l_divend_v}    = ny + 1 Using Progression prog_r_noz;          // 既設定と同じ
Transfinite Curve {l_mid_v_c}     = ny + 1 Using Progression 1.00;
Transfinite Curve {l_outer3_bot}  = ny + 1 Using Progression (1.0/1.0);

// r方向 (lower outer, r_exit→r_mid)
Transfinite Curve {l_outer1_bot}  = ny_out_bot + 1 Using Progression prog_r_bot;
Transfinite Curve {l_mid_v_d}     = ny_out_bot + 1 Using Progression prog_r_bot2;
Transfinite Curve {l_outer3_lo}   = ny_out_bot + 1 Using Progression (1.0/prog_r_bot3);

// r方向 (upper outer, r_mid→r_far)
Transfinite Curve {l_outer1_top}  = ny_out_top + 1 Using Progression prog_r_top;
Transfinite Curve {l_mid_v_e}     = ny_out_top + 1 Using Progression prog_r_top;
Transfinite Curve {l_outer3_hi}   = ny_out_top + 1 Using Progression (1.0/prog_r_top);

// ---- Block C1: inner left (r=0→r_exit, x=div_end_x→x_mid) ----
ll_c1 = newll;
Curve Loop(ll_c1) = {l_axis_far_1, l_mid_v_c, -l_interface_1, -l_divend_v};
s_c1 = news; Plane Surface(s_c1) = {ll_c1};
Transfinite Surface {s_c1} = {p_axis_div_end, p_x_mid, p_interface_mid, p_wall_div_end};
Recombine Surface {s_c1};

// ---- Block C2: inner right (r=0→r_exit, x=x_mid→x_far) ----
ll_c2 = newll;
Curve Loop(ll_c2) = {l_axis_far_2, -l_outer3_bot, -l_interface_2, -l_mid_v_c};
s_c2 = news; Plane Surface(s_c2) = {ll_c2};
Transfinite Surface {s_c2} = {p_x_mid, p_axis_exit, p_exit_mid, p_interface_mid};
Recombine Surface {s_c2};

// ---- Block D1: lower outer left (r=r_exit→r_mid, x=div_end_x→x_mid) ----
ll_d1 = newll;
Curve Loop(ll_d1) = {l_interface_1, l_mid_v_d, -l_r_mid_1, -l_outer1_bot};
s_d1 = news; Plane Surface(s_d1) = {ll_d1};
Transfinite Surface {s_d1} = {p_wall_div_end, p_interface_mid, p_inner_corner, p_outer1_mid};
Recombine Surface {s_d1};

// ---- Block D2: lower outer right (r=r_exit→r_mid, x=x_mid→x_far) ----
ll_d2 = newll;
Curve Loop(ll_d2) = {l_interface_2, -l_outer3_lo, -l_r_mid_2, -l_mid_v_d};
s_d2 = news; Plane Surface(s_d2) = {ll_d2};
Transfinite Surface {s_d2} = {p_interface_mid, p_exit_mid, p_exit_r_mid, p_inner_corner};
Recombine Surface {s_d2};

// ---- Block E1: upper outer left (r=r_mid→r_far, x=div_end_x→x_mid) ----
ll_e1 = newll;
Curve Loop(ll_e1) = {l_r_mid_1, l_mid_v_e, -l_outer2_1, -l_outer1_top};
s_e1 = news; Plane Surface(s_e1) = {ll_e1};
Transfinite Surface {s_e1} = {p_outer1_mid, p_inner_corner, p_outer2_mid, p_outer1_top};
Recombine Surface {s_e1};

// ---- Block E2: upper outer right (r=r_mid→r_far, x=x_mid→x_far) ----
ll_e2 = newll;
Curve Loop(ll_e2) = {l_r_mid_2, -l_outer3_hi, -l_outer2_2, -l_mid_v_e};
s_e2 = news; Plane Surface(s_e2) = {ll_e2};
Transfinite Surface {s_e2} = {p_inner_corner, p_exit_r_mid, p_exit_top, p_outer2_mid};
Recombine Surface {s_e2};

// ==============================================================
// Outer 延伸領域 — L字形三角格子
// ==============================================================
// x方向: x_far  → x_far2 = div_end_x + 20*D_exit  (プルーム領域 2倍)
// r方向: r_far  → r_far2 = r_exit   + 10*D_exit  (半径 2倍)
//
//  r=r_far2  p_tl ──────────────── p_tm ──────────────── p_tr
//             │   [triangle 上ストリップ]  │  [triangle 右ストリップ上]  │
//  r=r_far   p_outer1_top ─ 既存quad ─ p_exit_top         │
//             │          [既存構造格子]            │         │
//  r=0       axis_div_end ─────────── p_axis_exit ── p_br
//
// 既存の outer2 (r=r_far) / outer3 (x=x_far) は内部面になるため
// Physical Curve から除外し、新しい外部境界に差し替える。

r_far2 = r_exit + 10.0 * D_exit;     // 半径方向 2倍 (= r_far + 5*D_exit)
x_far2 = div_end_x + 20.0 * D_exit;  // 流れ方向 2倍 (= x_far + 10*D_exit)

// 格子サイズ目安:
//   既存 outer3 の r=r_far 側セル ≈ 22 mm
//   既存 outer2 の x=x_far 側セル ≈ 19 mm
//   → 界面付近は lc_near = 25 mm に設定して滑らかにつなぐ
lc_near = 25.0;
lc_far  = 70.0;

// 新頂点
p_tl = newp; Point(p_tl) = {div_end_x, r_far2, 0.0, lc_near};
p_tm = newp; Point(p_tm) = {x_far,     r_far2, 0.0, lc_near};
p_tr = newp; Point(p_tr) = {x_far2,    r_far2, 0.0, lc_far};
p_br = newp; Point(p_br) = {x_far2,    0.0,    0.0, lc_far};

// 新境界線
l_new_left   = newl; Line(l_new_left)   = {p_outer1_top, p_tl};         // x=div_end_x, r_far→r_far2
l_new_top_L  = newl; Line(l_new_top_L)  = {p_tl,         p_tm};         // r=r_far2, 左半
l_new_top_R  = newl; Line(l_new_top_R)  = {p_tm,         p_tr};         // r=r_far2, 右半
l_new_right  = newl; Line(l_new_right)  = {p_tr,         p_br};         // x=x_far2, 全高
l_new_bottom = newl; Line(l_new_bottom) = {p_br,         p_axis_exit};  // r=0, 軸延伸

// L字形曲線ループ (CCW)
//   外周: l_new_left → l_new_top_L → l_new_top_R → l_new_right → l_new_bottom
//   内周 (既存quad との界面):
//     x=x_far  ↑ : -l_outer3_bot, -l_outer3_lo, -l_outer3_hi
//     r=r_far  ← : -l_outer2_2, -l_outer2_1
ll_new = newll;
Curve Loop(ll_new) = {
  l_new_left, l_new_top_L, l_new_top_R, l_new_right, l_new_bottom,
  -l_outer3_bot, -l_outer3_lo, -l_outer3_hi,
  -l_outer2_2, -l_outer2_1
};
s_new = news; Plane Surface(s_new) = {ll_new};
// Transfinite/Recombine なし → 三角格子

// ==============================================================
// Physical groups
// ==============================================================
// outer2 (l_outer2_1/2) と outer3 (l_outer3_*) は quad と triangle の
// 内部接続面になるため Physical Curve から除外。
Physical Curve  ("inlet",  1) = {l_inlet};
Physical Curve  ("wall",   3) = {l_conv_top, l_div_top};
Physical Curve  ("axis",   4) = {l_conv_bot, l_div_bot, l_axis_far_1, l_axis_far_2, l_new_bottom};
Physical Curve  ("outer1", 5) = {l_outer1_bot, l_outer1_top, l_new_left};
Physical Curve  ("outer2", 6) = {l_new_top_L, l_new_top_R};
Physical Curve  ("outer3", 7) = {l_new_right};
Physical Surface("fluid",  8) = {s_conv, s_div, s_c1, s_c2, s_d1, s_d2, s_e1, s_e2, s_new};
