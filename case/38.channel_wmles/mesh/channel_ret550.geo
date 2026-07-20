// WMLES 周期チャネル Reτ=550 (Hoyas & Jiménez 較正用の最小構成)
//   δ=0.005, 箱 2πδ × 2δ × πδ, 一様格子 48×24×32 (Δx+=72, Δy+=46, Δz+=54 @ u_τ=3.85)
//   x/z 周期、y=0/2δ が等温壁 (wallModelLES)。体積力 f_x=ρu_τ²/δ で駆動。
Lx = 0.031415927;
Ly = 0.01;
Lz = 0.015707963;
nx = 48; ny = 24; nz = 32;

Point(1) = {0, 0, 0, 1};
Point(2) = {Lx, 0, 0, 1};
Point(3) = {Lx, Ly, 0, 1};
Point(4) = {0, Ly, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Curve {1, 3} = nx + 1;
Transfinite Curve {2, 4} = ny + 1;
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface {1};

ext[] = Extrude {0, 0, Lz} { Surface{1}; Layers{nz}; Recombine; };

// ext[]: [0]=z=Lz 面, [1]=volume, [2..5]=側面 (Line1..4 の掃引順)
Physical Surface("xlo",  1) = {ext[5]};   // x=0  (Line4 掃引)
Physical Surface("xhi",  2) = {ext[3]};   // x=Lx (Line2 掃引)
Physical Surface("ylo",  3) = {ext[2]};   // y=0  壁 (Line1 掃引)
Physical Surface("yhi",  4) = {ext[4]};   // y=Ly 壁 (Line3 掃引)
Physical Surface("zlo",  5) = {1};        // z=0
Physical Surface("zhi",  6) = {ext[0]};   // z=Lz
Physical Volume ("fluid", 7) = {ext[1]};
