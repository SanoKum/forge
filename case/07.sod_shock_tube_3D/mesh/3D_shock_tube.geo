Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 1.0;
//Geometry.Color.Surfaces = Black;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)

nx=100; ny=11;

// Place points
Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(3) = {1,0.05,0,lc};
Point(4) = {0,0.05,0,lc};

// Create lines from points
Line(1) = {1,2}; Transfinite Line {1} = nx;
Line(2) = {2,3}; Transfinite Line {2} = ny;
Line(3) = {3,4}; Transfinite Line {3} = nx;
Line(4) = {4,1}; Transfinite Line {4} = ny;

// Define line loops used to construct surfaces
Line Loop(1) = {1,2,3,4};

// Make surfaces from line loops
Plane Surface(1) = {1}; //using LL #1
Transfinite Surface {1};
Recombine Surface(1);

//+
Extrude {0, 0, 0.05} {
  Surface{1}; Layers{11}; Recombine;
}
//+
Physical Surface("bottom", 1) = {13};
Physical Surface("upper", 2) = {21};
Physical Surface("slip", 3) = {26, 1};
Physical Surface("inlet", 4) = {25};
Physical Surface("outlet", 5) = {17};

Physical Volume("fluid", 6) = {1};
