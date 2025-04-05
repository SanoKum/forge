Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 1.0;
//Geometry.Color.Surfaces = Black;

Mesh.Voronoi=1;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)

nx=128; ny=64;

// Place points
Point(1) = {0  ,0.0,0,lc};
Point(2) = {1.0,0.0,0,lc};
Point(3) = {1.5,0.1,0,lc};
Point(4) = {2.0,0.0,0,lc};
Point(5) = {3.0,0.0,0,lc};
Point(6) = {3.0,1.0,0,lc};
Point(7) = {2.0,1.0,0,lc};
Point(8) = {1.0,1.0,0,lc};
Point(9) = {0.0,1.0,0,lc};


// Create lines from points
Line(1) = {1,2}; Transfinite Line {1} = (nx-1)/3;

// circle arc
Point(10) = {1.5,-1.2,0,lc};
Circle(2) = {2, 10, 4}; Transfinite Line {2} = (nx-1)/3;

Line(3) = {4,5}; Transfinite Line {3} = (nx-1)/3;
Line(4) = {5,6}; Transfinite Line {4} = ny;
Line(5) = {6,7}; Transfinite Line {5} = (nx-1)/3;
Line(6) = {7,8}; Transfinite Line {6} = (nx-1)/3;
Line(7) = {8,9}; Transfinite Line {7} = (nx-1)/3;
Line(8) = {9,1}; Transfinite Line {8} = ny      ;

Line(9) = {2,8} ; Transfinite Line {9} = ny      ;
Line(10) = {4,7}; Transfinite Line {10} = ny      ;

//
//// Define line loops used to construct surfaces
Line Loop(1) = {1,9,7,8};
Line Loop(2) = {2,10,6,-9};
Line Loop(3) = {3,4,5,-10};
//
//// Make surfaces from line loops
Plane Surface(1) = {1}; 
//Transfinite Surface {1};
//Recombine Surface(1);

Plane Surface(2) = {2}; 
//Transfinite Surface {2};
//Recombine Surface(2);

Plane Surface(3) = {3}; 
//Transfinite Surface {3};
//Recombine Surface(3);

////+
Extrude {0, 0, 0.02} {
  Surface{1,2,3}; Layers{1}; Recombine;
}

Physical Surface("inlet", 1) = {31};
Physical Surface("top", 2) = {27, 49, 71};
Physical Surface("bottom", 3) = {19, 41, 63};
Physical Surface("outlet", 4) = {67};
Physical Surface("slip", 5) = {32, 1, 54, 2, 76, 3};
Physical Volume("fluid", 6) = {1, 2, 3};
