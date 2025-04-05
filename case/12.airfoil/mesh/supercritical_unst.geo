
Include "airfoil.geo";

Mesh.Voronoi=1;

//+

ymax = 10;
xmax = 15;
n_inlet = 130;
n_vertical = 90;
r_vertical = 1/0.92;
n_airfoil = 80;
n_wake = 200;
r_wake1 = 1/0.985;
r_wake2 = 1.0;


Point(129) = {-3.0, ymax, 0, 1};
Point(130) = {-3.0, -ymax, 0, 1};
Point(131) = {2.0, ymax, 0, 1};
Point(132) = {2.0, -ymax, 0, 1};
Point(133) = {xmax, ymax, 0, 1};
Point(134) = {xmax, -ymax, 0, 1};
Point(135) = {xmax, 0, 0, 1};
Circle(2) = {130, 64, 129};
Line(3) = {57, 129};
Line(4) = {71, 130};
Line(5) = {129, 131};
Line(6) = {130, 132};
Line(7) = {131, 133};
Line(8) = {132, 134};
Line(9) = {135, 134};
Line(10) = {135, 133};
Line(11) = {128, 131};
Line(12) = {128, 132};
Line(13) = {128, 135};
Split Curve {1} Point {57, 71};
Split Curve {15} Point {128};

//Transfinite Curve {2, 14} = n_inlet Using Progression 1;
//Transfinite Curve {3, 11, 10, 4, 12, 9, 9} = n_vertical Using Progression r_vertical;
//Transfinite Curve {17, 16} = n_airfoil Using Bump 0.1;
//Transfinite Curve {5, 6} = n_airfoil Using Bump 2;
//Transfinite Curve {13} = n_wake Using Progression r_wake1;
//Transfinite Curve {7, 8} = n_wake Using Progression r_wake2;

Curve Loop(1) = {2, -3, 14, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {3, 5, -11, 17};
Plane Surface(2) = {2};
Curve Loop(3) = {11, 7, -10, -13};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 6, -12, -16};
Plane Surface(4) = {4};
Curve Loop(5) = {12, 8, -9, -13};
Plane Surface(5) = {5};
//+
//Transfinite Surface {1};
////+
//Transfinite Surface {2};
////+
//Transfinite Surface {3};
////+
//Transfinite Surface {5};
////+
//Transfinite Surface {4};
////+
//Recombine Surface {1, 2, 3, 5, 4};
//+
Extrude {0, 0, 0.05} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Layers {1}; Recombine;
}
//+
Physical Volume("Fluid", 128) = {1, 2, 3, 5, 4};
//+
Physical Surface("Inlet", 129) = {26, 52, 74, 96, 118};
//+
Physical Surface("Outlet", 130) = {78, 122};
//+
Physical Surface("Airfoil", 131) = {60, 104, 34};
//+
Physical Surface("Side", 132) = {83, 127, 105, 39, 61, 1, 2, 3, 5, 4};
//+

Field[131] = BoundaryLayer;
Field[131].AnisoMax = 1000;
Field[131].Thickness = 5;
Field[131].CurvesList = {1};
Field[131].Thickness = 3;
Field[131].NbLayers = 10;
Field[131].PointsList = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128};
Field[131].Quads = 1;
Field[131].Ratio = 1.1;
Field[131].Size = 0.05;
Field[131].SizeFar = 2;

BoundaryLayer Field = 1;

Show "*";
