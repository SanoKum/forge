Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 1.0;
//Geometry.Color.Surfaces = Black;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)

r1 = 0.0100;
r2 = 0.025 ;

theta = 15.0/180.0*Pi;

nr = 80;
nt = 80;

// Place points
Point(1) = {0              ,0.0            ,0,lc};
Point(2) = {-r1*Sin(theta) ,r1*Cos(theta)  ,0,lc};
Point(3) = {-r2*Sin(theta) ,r2*Cos(theta)  ,0,lc};
Point(4) = {-r2            ,0.0            ,0,lc};
Point(5) = {-r2*Sin(theta) ,-r2*Cos(theta) ,0,lc};
Point(6) = {-r1*Sin(theta) ,-r1*Cos(theta) ,0,lc};
Point(7) = {-r1            ,0.0            ,0,lc};

// Create lines from points
Line(1) = {2,3}; Transfinite Line {1} = nr;

Circle(2) = {2, 1, 7}; Transfinite Line {2} = nt;
Circle(3) = {7, 1, 6}; Transfinite Line {3} = nt;

Circle(4) = {3, 1, 4}; Transfinite Line {4} = nt;
Circle(5) = {4, 1, 5}; Transfinite Line {5} = nt;

Line(6) = {5,6}; Transfinite Line {6} = nr;

Line(7) = {4,7}; Transfinite Line {7} = nr;

Curve Loop(1) = {2, -7, -4, -1};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface(1);

Curve Loop(2) = {5, 6, -3, -7};
Plane Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface(2);


Extrude {0, 0, 0.001} {
  Surface{1,2}; Layers{1}; Recombine;
}

Physical Surface("inlet", 1) = {24, 38};
Physical Surface("outlet", 2) = {28, 42};
Physical Surface("wall", 3) = {16, 46};
Physical Surface("perio1", 4) = {29, 51};
Physical Surface("perio2", 5) = {1, 2};
Physical Volume("fluid", 6) = {1, 2};

//
////Physical Surface("inlet", 1) = {31};
////Physical Surface("top", 2) = {27, 49, 71};
////Physical Surface("bottom", 3) = {19, 41, 63};
////Physical Surface("outlet", 4) = {67};
////Physical Surface("slip", 5) = {32, 1, 54, 2, 76, 3};
//+

//+

//+

