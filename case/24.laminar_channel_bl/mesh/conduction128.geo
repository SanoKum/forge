// 等温壁検証用の純伝導ボックス (静止流体・厳密解あり)
//
//   y=H  +--------------+   <- top wall (wall_isothermal, T_top,  physID 4)
//        |              |
//  slip  |   no flow    | slip
//   (1)  |              |  (2)
//   y=0  +--------------+   <- bottom wall (wall_isothermal, T_bot, physID 3)
//        x=0          x=L
//
// 定常厳密解: T(y) = T_bot + (T_top - T_bot) * y/H (線形),
//             q_w = lambda * (T_bot - T_top) / H  (thermCond 定数 = viscMethod:0)
// channel.geo と同じ壁法線クラスタリング (Bump 0.4) を残し、
// 非一様間隔でも熱流束評価が正しいことを併せて確認する。

L  = 0.02;   // 幅 (解は y のみの 1D)    [m]
H  = 0.01;   // 壁間距離                 [m]
nx = 16;
ny = 128;

Point(1) = {0, 0, 0, 1};
Point(2) = {L, 0, 0, 1};
Point(3) = {L, H, 0, 1};
Point(4) = {0, H, 0, 1};

Line(1) = {1, 2};   // bottom wall
Line(2) = {2, 3};   // right side (slip)
Line(3) = {3, 4};   // top wall
Line(4) = {4, 1};   // left side (slip)

Transfinite Curve {1, 3} = nx + 1;
Transfinite Curve {2, 4} = ny + 1 Using Bump 0.4;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface {1};

Physical Curve  ("side_left",  1) = {4};
Physical Curve  ("side_right", 2) = {2};
Physical Curve  ("wall_bot",   3) = {1};
Physical Curve  ("wall_top",   4) = {3};
Physical Surface("fluid",      5) = {1};
