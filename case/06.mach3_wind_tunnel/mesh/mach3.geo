Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 1.0;
//Geometry.Color.Surfaces = Black;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)

nx=600/2; ny=200/2;

// Place points
Point(1) = {0  ,0.0,0,lc};
Point(2) = {0.6,0.0,0,lc};
Point(3) = {0.6,0.2,0,lc};
Point(4) = {3.0,0.2,0,lc};
Point(5) = {3.0,1.0,0,lc};
Point(6) = {0.6,1.0,0,lc};
Point(7) = {0.0,1.0,0,lc};
Point(8) = {0.0,0.2,0,lc};

// Create lines from points
Line(1) = {1,2}; Transfinite Line {1} = nx/5;
Line(2) = {2,3}; Transfinite Line {2} = ny/5;
Line(3) = {3,4}; Transfinite Line {3} = nx*4/5;
Line(4) = {4,5}; Transfinite Line {4} = ny*4/5;
Line(5) = {5,6}; Transfinite Line {5} = nx*4/5;
Line(6) = {6,7}; Transfinite Line {6} = nx/5;
Line(7) = {7,8}; Transfinite Line {7} = ny*4/5;
Line(8) = {8,1}; Transfinite Line {8} = ny/5;

Line(9) = {3,6}; Transfinite Line {9} = ny*4/5;
Line(10) = {8,3}; Transfinite Line {10} = nx/5;

// Define line loops used to construct surfaces
Line Loop(1) = {1,2,-10,8};
Line Loop(2) = {3,4,5,-9};
Line Loop(3) = {6,7,9,10};

// Make surfaces from line loops
Plane Surface(1) = {1}; 
Transfinite Surface {1};
Recombine Surface(1);

Plane Surface(2) = {2}; 
Transfinite Surface {2};
Recombine Surface(2);

Plane Surface(3) = {3}; 
Transfinite Surface {3};
Recombine Surface(3);

//+
Extrude {0, 0, 0.04} {
  Surface{1,2,3}; Layers{1}; Recombine;
}

//+
Physical Surface("inlet", 1) = {67, 31};
Physical Surface("top", 2) = {63, 49};
Physical Surface("bottom", 3) = {19};
Physical Surface("wall", 4) = {23, 41};
Physical Surface("outlet", 5) = {45};
Physical Surface("slip1", 6) = {1, 3, 2};
Physical Surface("slip2", 7) = {54, 32, 76};
Physical Volume("fluid", 8) = {1, 2, 3};
