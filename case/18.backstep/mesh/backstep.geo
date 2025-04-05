Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 1.0 ; // 
//Geometry.Color.Surfaces = Black;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)


// -------------------
// *** (1) section ***
// -------------------
nx_1 = 15;
nx_2 = 215;

ny_1 = 17;
ny_2 = 33;

H = 1.0 ;

Point(1) = {0   ,   H, 0, lc};
Point(2) = {2*H ,   H, 0, lc};
Point(3) = {2*H ,   0, 0, lc};
Point(4) = {30*H, 0.0, 0, lc};
Point(5) = {30*H,   H, 0, lc};
Point(6) = {30*H, 3*H, 0, lc};
Point(7) = {2*H , 3*H, 0, lc};
Point(8) = {0   , 3*H, 0, lc};

l1= newl;
Line(l1) =  {1, 2}; Transfinite Line {l1} = nx_1;

l2= newl;
Line(l2) =  {2, 3}; Transfinite Line {l2} = ny_1;

l3= newl;
Line(l3) =  {3, 4}; Transfinite Line {l3} = nx_2;

l4= newl;
Line(l4) =  {4, 5}; Transfinite Line {l4} = ny_1;

l5= newl;
Line(l5) =  {5, 6}; Transfinite Line {l5} = ny_2;

l6= newl;
Line(l6) =  {6, 7}; Transfinite Line {l6} = nx_2;

l7= newl;
Line(l7) =  {7, 8}; Transfinite Line {l7} = nx_1;

l8= newl;
Line(l8) =  {7, 2}; Transfinite Line {l8} = ny_2;

l9= newl;
Line(l9) =  {8, 1}; Transfinite Line {l9} = ny_2;

l10= newl;
Line(l10) = {2, 5}; Transfinite Line {l10} = nx_2;

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


//
//Curve Loop(1) = {10, 1, 11, 9};
//Plane Surface(1) = {1};
//Transfinite Surface {1};
//Recombine Surface(1);
//
//Curve Loop(2) = {11, -8, -12, -2};
//Plane Surface(2) = {2};
//Transfinite Surface {2};
//Recombine Surface(2);
//
//
//Curve Loop(3) = {12, -7, -13, -3};
//Plane Surface(3) = {3};
//Transfinite Surface {3};
//Recombine Surface(3);
//
//
//Curve Loop(4) = {13, -6, -5, -4};
//Plane Surface(4) = {4};
//Transfinite Surface {4};
//Recombine Surface(4);
//
//Extrude {0, 0, 1.259} {
//  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Layers {5}; Recombine;
//}
//
//Physical Surface("inlet", 1) = {22};
//Physical Surface("outlet", 2) = {96};
//Physical Surface("wall", 3) = {1, 26, 35, 34, 48, 2, 56, 57, 70, 3, 78, 79, 4, 92, 100, 101};
//Physical Volume("fluid", 4) = {1, 2, 3, 4};
////+
//+
Physical Surface("inlet", 1) = {19};
Physical Surface("outlet", 2) = {49, 71};
Physical Surface("top", 3) = {53, 31};
Physical Surface("bot", 4) = {63, 23, 67};
Physical Surface("side2", 5) = {54, 32, 76};
Physical Surface("side1", 6) = {2, 3, 1};
Physical Volume("fluid", 7) = {1, 2, 3};