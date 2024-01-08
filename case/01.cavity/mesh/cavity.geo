Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 1.0;
//Geometry.Color.Surfaces = Black;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)

nx=36; ny=36;

// Place points
Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(3) = {1,1,0,lc};
Point(4) = {0,1,0,lc};

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
Extrude {0, 0, 0.1} {
  Surface{1}; Layers{1}; Recombine;
}
//+
Physical Surface("wall", 1) = {25, 13, 17};
//+
Physical Surface("slip", 2) = {26, 1};
//+
Physical Surface("upper", 3) = {21};
//+
Physical Volume("fluid", 4) = {1};
