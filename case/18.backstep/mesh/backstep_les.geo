Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Mesh.ScalingFactor = 1.0 ;

lc=0.1; // characteristic mesh size (transfinite なので実質未使用)

// -------------------
// *** backstep LES mesh: 230(x) x 50(y) x 80(span z) ***
// 段差 x=2H, H=1。x[0,30] y[0,3] z[0,4]。
// x: nx_1=16(15cell, x0-2) + nx_2=216(215cell, x2-30) = 230 cells
// y: ny_1=18(17cell, y0-1) + ny_2=34(33cell, y1-3)    = 50 cells
// z: Layers{80}                                        = 80 cells (Δz=0.05H)
// -------------------
nx_1 = 16;
nx_2 = 216;

ny_1 = 18;
ny_2 = 34;

H = 1.0 ;

Point(1) = {0   ,   H, 0, lc};
Point(2) = {2*H ,   H, 0, lc};
Point(3) = {2*H ,   0, 0, lc};
Point(4) = {30*H, 0.0, 0, lc};
Point(5) = {30*H,   H, 0, lc};
Point(6) = {30*H, 3*H, 0, lc};
Point(7) = {2*H , 3*H, 0, lc};
Point(8) = {0   , 3*H, 0, lc};

l1= newl; Line(l1) =  {1, 2}; Transfinite Line {l1} = nx_1;
l2= newl; Line(l2) =  {2, 3}; Transfinite Line {l2} = ny_1;
l3= newl; Line(l3) =  {3, 4}; Transfinite Line {l3} = nx_2;
l4= newl; Line(l4) =  {4, 5}; Transfinite Line {l4} = ny_1;
l5= newl; Line(l5) =  {5, 6}; Transfinite Line {l5} = ny_2;
l6= newl; Line(l6) =  {6, 7}; Transfinite Line {l6} = nx_2;
l7= newl; Line(l7) =  {7, 8}; Transfinite Line {l7} = nx_1;
l8= newl; Line(l8) =  {7, 2}; Transfinite Line {l8} = ny_2;
l9= newl; Line(l9) =  {8, 1}; Transfinite Line {l9} = ny_2;
l10= newl; Line(l10) = {2, 5}; Transfinite Line {l10} = nx_2;

Curve Loop(1) = {9, 1, -8, 7};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface(1);

Curve Loop(2) = {8, 10, 5, 6};
Plane Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface(2);

Curve Loop(3) = {2, 3, 4, -10};
Plane Surface(3) = {3};
Transfinite Surface {3};
Recombine Surface(3);

Extrude {0, 0, 4*H} {
  Surface{1}; Surface{2}; Surface{3}; Layers {80}; Recombine;
}

Physical Surface("inlet", 1) = {19};
Physical Surface("outlet", 2) = {49, 71};
Physical Surface("top", 3) = {53, 31};
Physical Surface("bot", 4) = {63, 23, 67};
Physical Surface("side2", 5) = {54, 32, 76};
Physical Surface("side1", 6) = {2, 3, 1};
Physical Volume("fluid", 7) = {1, 2, 3};
