// =====================================================================
// 乱流平板 (SST 壁法則検証) メッシュ
//   - 領域: x in [-0.1, 1.0], y in [0, 0.2]  ** 平面 2D (押し出しなし) **
//   - 底面: x in [-0.1, 0] は slip (対称), x in [0, 1.0] は no-slip wall (平板)
//   - 入口 x=-0.1 (全圧全温), 出口 x=1.0 (背圧), 上面 y=0.2 は slip
//   - 壁垂直方向: ny=45, r_y=1.09660 -> mid-plate で y+~30 (壁関数帯)
//   - 流体: 空気 mu=1.8e-5, U~34 m/s (M=0.1), Re/m ~ 2.3e6
// =====================================================================

Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
General.Color.Text = White;
Mesh.ScalingFactor = 1.0;

lc = 0.05;

// --- 幾何パラメータ ---
x0 = -0.1;   // 入口 (上流 slip 区間の始点)
x1 =  0.0;   // 平板前縁 (slip -> wall の境界)
x2 =  1.0;   // 平板後縁 = 出口
H  =  0.2;   // 領域高さ
// (z 厚さなし: 平面 2D)

// --- 分割数 / 等比 ---
nx_up    = 40;     // 上流 slip 区間 (streamwise)
r_up     = 1.08;   // 上流も前縁 (x=0) 側を細かく寄せる (前縁でセルサイズ連続)
nx_plate = 200;    // 平板区間 (streamwise)
r_plate  = 1.02;   // 前縁側を細かく
ny       = 45;
r_y      = 1.09660;

// --- 点 (底面: 前縁で分割) ---
Point(1) = {x0, 0.0, 0.0, lc};   // 入口下
Point(2) = {x1, 0.0, 0.0, lc};   // 前縁
Point(3) = {x2, 0.0, 0.0, lc};   // 後縁下 (出口下)
Point(4) = {x2, H,   0.0, lc};   // 後縁上 (出口上)
Point(5) = {x1, H,   0.0, lc};   // 前縁上
Point(6) = {x0, H,   0.0, lc};   // 入口上

// --- 線 ---
Line(1) = {1, 2};   // bottom slip (上流)
Line(2) = {2, 3};   // bottom wall (平板)
Line(3) = {3, 4};   // outlet (右)
Line(4) = {4, 5};   // top plate
Line(5) = {5, 2};   // 前縁の垂直線 (共有, 壁->上)  ※下=点2
Line(6) = {5, 6};   // top up
Line(7) = {6, 1};   // inlet (左)

// --- 壁垂直方向 transfinite (壁側=底で細かく) ---
// 線7 は 6(上)->1(下) なので始点が上。壁(下)を細かくするため Progression を負側に。
Transfinite Line {7} = ny Using Progression 1.0/r_y;   // 6(上)->1(底=壁): 終点(壁)で細かく
Transfinite Line {5} = ny Using Progression 1.0/r_y;   // 5(上)->2(底=壁): 終点(壁)で細かく
Transfinite Line {3} = ny Using Progression r_y;       // 3(底=壁)->4(上): 始点(壁)で細かく

// --- streamwise transfinite ---
// 上流も前縁 (点2 / 点5 = x=0) 側へ寄せる。前縁で plate 側とセルサイズを連続させる。
Transfinite Line {1} = nx_up Using Progression 1.0/r_up;  // 1(入口)->2(前縁): 終点(前縁)で細かく
Transfinite Line {6} = nx_up Using Progression r_up;      // 5(前縁)->6(入口): 始点(前縁)で細かく
Transfinite Line {2} = nx_plate Using Progression r_plate;       // 2(前縁)->3(後縁): 前縁側細かく
Transfinite Line {4} = nx_plate Using Progression 1.0/r_plate;   // 4(後縁)->5(前縁): 終点(前縁)で細かく

// --- 面 ---
Curve Loop(1) = {7, 1, -5, 6};       // 上流ブロック: inlet,bottom_slip,(2->5),top_up
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface(1);

Curve Loop(2) = {5, 2, 3, 4};        // 平板ブロック: (5->2),bottom_wall,outlet,top_plate
Plane Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface(2);

// =====================================================================
// 平面 2D (node-centered 用): 押し出しなし, Physical Curve でタグ付け
//   node-centered (median-dual) は 1 層押し出しメッシュだと spanwise に
//   ノードが 2 枚できる。2 点しかない方向では 2 次 MUSCL の左右再構成が
//   厳密に一致して迎角上流差分の散逸が完全に消え、spanwise 市松モードが
//   無減衰になる (case/26 run_0032-0037 の 2 次発散の真因)。2D 問題は
//   必ずこの平面メッシュを使うこと。
// =====================================================================
Physical Curve("inlet",  1) = {7};       // x=x0 入口
Physical Curve("outlet", 2) = {3};       // x=x2 出口
Physical Curve("top",    3) = {4, 6};    // y=H 上面 slip
Physical Curve("wall",   4) = {2};       // y=0, x in [0,1] 平板 no-slip
Physical Curve("sym",    5) = {1};       // y=0, x in [-0.1,0] slip
Physical Surface("fluid", 8) = {1, 2};
