// 2D 版 bump (node-centered 検証用)。bump.geo の Extrude を除き、
// 2D 平面メッシュ + Physical Curve 境界 (inlet/top/bottom/outlet) を定義する。
// 構造化 quad (Transfinite + Recombine)。物理 ID は bcondConfig (loM/hiM) に一致:
//   1:inlet  2:top  3:bottom  4:outlet  (3D の front/back symmetry physID5 は 2D では無し)
Geometry.LineNumbers = 1;
Mesh.ScalingFactor = 1.0;

lc=0.025;
nx=129; ny=49;

Point(1) = {0  ,0.0,0,lc};
Point(2) = {1.0,0.0,0,lc};
Point(3) = {1.5,0.1,0,lc};
Point(4) = {2.0,0.0,0,lc};
Point(5) = {3.0,0.0,0,lc};
Point(6) = {3.0,1.0,0,lc};
Point(7) = {2.0,1.0,0,lc};
Point(8) = {1.0,1.0,0,lc};
Point(9) = {0.0,1.0,0,lc};

Line(1) = {1,2}; 
Point(10) = {1.5,-1.2,0,lc};
Circle(2) = {2, 10, 4}; 
Line(3) = {4,5}; 
Line(4) = {5,6}; 
Line(5) = {6,7}; 
Line(6) = {7,8}; 
Line(7) = {8,9}; 
Line(8) = {9,1}; 

Line(9) = {2,8} ; 
Line(10) = {4,7}; 

Line Loop(1) = {1,9,7,8};
Line Loop(2) = {2,10,6,-9};
Line Loop(3) = {3,4,5,-10};

Plane Surface(1) = {1}; 
Plane Surface(2) = {2}; 
Plane Surface(3) = {3}; 

// 2D 境界 (Physical Curve)
Physical Curve("inlet", 1)  = {8};        // x=0 左端
Physical Curve("top", 2)    = {5,6,7};    // y=1 上壁
Physical Curve("bottom", 3) = {1,2,3};    // y=0 下壁 (bump 含む)
Physical Curve("outlet", 4) = {4};        // x=3 右端
Physical Surface("fluid", 6) = {1,2,3};
