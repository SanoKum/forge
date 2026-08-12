// =====================================================================
// 乱流平板 (SST 壁法則検証) メッシュ
//   - 領域: x in [-0.1, 1.0], y in [0, 0.2], z 厚さ 0.01 (1 層押し出し)
//   - 底面: x in [-0.1, 0] は slip (対称), x in [0, 1.0] は no-slip wall (平板)
//   - 入口 x=-0.1 (全圧全温), 出口 x=1.0 (背圧), 上面 y=0.2 は slip
//   - 壁垂直方向: 第一セル高さ ~4e-6 m, 等比 1.10 -> mid-plate で y+~0.4 (<1)
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
t  =  0.01;  // z 厚さ

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

// --- z 方向に 1 層押し出し ---
Extrude {0, 0, t} {
  Surface{1}; Surface{2}; Layers{4}; Recombine;
}

// =====================================================================
// Physical Surface を bounding box で頑健にタグ付け
// =====================================================================
eps = 1e-6;

// inlet: x = x0 平面
Physical Surface("inlet", 1) = { Surface In BoundingBox{
  x0-eps, -eps, -eps,   x0+eps, H+eps, t+eps} };

// outlet: x = x2 平面
Physical Surface("outlet", 2) = { Surface In BoundingBox{
  x2-eps, -eps, -eps,   x2+eps, H+eps, t+eps} };

// top: y = H 平面 (上流 + 平板の 2 面)
Physical Surface("top", 3) = { Surface In BoundingBox{
  x0-eps, H-eps, -eps,   x2+eps, H+eps, t+eps} };

// wall: y = 0, x in [x1, x2] (平板のみ。上流 slip は除外)
Physical Surface("wall", 4) = { Surface In BoundingBox{
  x1-eps, -eps, -eps,   x2+eps, eps, t+eps} };

// sym: y = 0, x in [x0, x1] (上流 slip 区間)
Physical Surface("sym", 5) = { Surface In BoundingBox{
  x0-eps, -eps, -eps,   x1+eps, eps, t+eps} };

// side1: z = 0 平面 (元の 2 面)
Physical Surface("side1", 6) = { Surface In BoundingBox{
  x0-eps, -eps, -eps,   x2+eps, H+eps, eps} };

// side2: z = t 平面
Physical Surface("side2", 7) = { Surface In BoundingBox{
  x0-eps, -eps, t-eps,   x2+eps, H+eps, t+eps} };

Physical Volume("fluid", 8) = { Volume{:} };
