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
nx_1 = 10;
nx_2 = 100;
nx_3 = 100;

ny_1 = 20;
ny_2 = 50;

rat_y_1 = 1.00;
rat_y_2_1 = 1.08;
rat_y_2_2 = 1.02;
rat_y_2_3 = 1.00;

H = 1.0 ;

Point(1) = {0       , -0.375  , 0, lc};
Point(2) = {0.07    , -0.375  , 0, lc};
Point(3) = {0.96351 , -0.2944 , 0, lc};
Point(4) = {1.98663 , -0.1369 , 0, lc};
Point(5) = {0       , -0.3569 , 0, lc};
Point(6) = {0.06761 , -0.3569 , 0, lc};
Point(7) = {0.93068 , -0.1256 , 0, lc};
Point(8) = {1.91894 , 0.13918 , 0, lc};
Point(9) = {0       , 0.553298, 0, lc};
Point(10) = {0.067615, 0.553298, 0, lc};
Point(11)= {0.85    , 0.553298, 0, lc};
Point(12)= {1.817398, 0.553298, 0, lc};

l1= newl;
Line(l1) =  {1, 2}; Transfinite Line {l1} = nx_1;

l2= newl;
Line(l2) =  {2, 3}; Transfinite Line {l2} = nx_2;

l3= newl;
Line(l3) =  {3, 4}; Transfinite Line {l3} = nx_3;

l4= newl;
Line(l4) =  {5, 6}; Transfinite Line {l4} = nx_1;

l5= newl;
Line(l5) =  {6, 7}; Transfinite Line {l5} = nx_2;

l6= newl;
Line(l6) =  {7, 8}; Transfinite Line {l6} = nx_3;

l7= newl;
Line(l7) =  {9, 10}; Transfinite Line {l7} = nx_1;

l8= newl;
Line(l8) =  {10, 11}; Transfinite Line {l8} = nx_2;

l9= newl;
Line(l9) =  {11, 12}; Transfinite Line {l9} = nx_3;


l10= newl;
Line(l10) =  {1, 5}; Transfinite Line {l10} = ny_1 Using Progression rat_y_1;
//Line(l10) =  {1, 5}; Transfinite Line {l10} = ny_1 Using Progression 1.1;

l11= newl;
Line(l11) =  {2, 6}; Transfinite Line {l11} = ny_1;

l12= newl;
Line(l12) =  {3, 7}; Transfinite Line {l12} = ny_1;

l13= newl;
Line(l13) =  {4, 8}; Transfinite Line {l13} = ny_1;



l14= newl;
Line(l14) =  {5, 9}; Transfinite Line {l14} = ny_2 Using Progression rat_y_2_1;

l15= newl;
Line(l15) =  {6, 10}; Transfinite Line {l15} = ny_2 Using Progression rat_y_2_1;

l16= newl;
Line(l16) =  {7, 11}; Transfinite Line {l16} = ny_2 Using Progression rat_y_2_2;

l17= newl;
Line(l17) =  {8, 12}; Transfinite Line {l17} = ny_2 Using Progression rat_y_2_3;


//Curve Loop(1) = {9, 1, -8, 7};
//Plane Surface(1) = {1};
//Transfinite Surface {1};
//Recombine Surface(1);
//
//Curve Loop(2) = {8, 10, 5, 6};
//Plane Surface(2) = {2};
//Transfinite Surface {2};
//Recombine Surface(2);
//
//
//Curve Loop(3) = {2, 3, 4, -10};
//Plane Surface(3) = {3};
//Transfinite Surface {3};
//Recombine Surface(3);
//
//
//Extrude {0, 0, 4*H} {
//  Surface{1}; Surface{2}; Surface{3}; Layers {80}; Recombine;
//}
//
//
////
////Curve Loop(1) = {10, 1, 11, 9};
////Plane Surface(1) = {1};
////Transfinite Surface {1};
////Recombine Surface(1);
////
////Curve Loop(2) = {11, -8, -12, -2};
////Plane Surface(2) = {2};
////Transfinite Surface {2};
////Recombine Surface(2);
////
////
////Curve Loop(3) = {12, -7, -13, -3};
////Plane Surface(3) = {3};
////Transfinite Surface {3};
////Recombine Surface(3);
////
////
////Curve Loop(4) = {13, -6, -5, -4};
////Plane Surface(4) = {4};
////Transfinite Surface {4};
////Recombine Surface(4);
////
////Extrude {0, 0, 1.259} {
////  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Layers {5}; Recombine;
////}
////
////Physical Surface("inlet", 1) = {22};
////Physical Surface("outlet", 2) = {96};
////Physical Surface("wall", 3) = {1, 26, 35, 34, 48, 2, 56, 57, 70, 3, 78, 79, 4, 92, 100, 101};
////Physical Volume("fluid", 4) = {1, 2, 3, 4};
//////+
////+
//Physical Surface("inlet", 1) = {19};
//Physical Surface("outlet", 2) = {49, 71};
//Physical Surface("top", 3) = {53, 31};
//Physical Surface("bot", 4) = {63, 23, 67};
//Physical Surface("side2", 5) = {54, 32, 76};
//Physical Surface("side1", 6) = {2, 3, 1};
//Physical Volume("fluid", 7) = {1, 2, 3};//+
Curve Loop(1) = {10, 4, -11, -1};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface(1);


Curve Loop(2) = {2, 12, -5, -11};
Plane Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface(2);


Curve Loop(3) = {3, 13, -6, -12};
Plane Surface(3) = {3};
Transfinite Surface {3};
Recombine Surface(3);


Curve Loop(4) = {15, -7, -14, 4};
Plane Surface(4) = {4};
Transfinite Surface {4};
Recombine Surface(4);


Curve Loop(5) = {15, 8, -16, -5};
Plane Surface(5) = {5};
Transfinite Surface {5};
Recombine Surface(5);



Curve Loop(6) = {16, 9, -17, -6};
Plane Surface(6) = {6};
Transfinite Surface {6};
Recombine Surface(6);

Extrude {0, 0, 0.005} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Layers {1}; Recombine;
}


//+
Physical Surface("inlet", 1) = {100, 26};
Physical Surface("outlet", 2) = {144, 74};
Physical Surface("wall", 3) = {70, 48};
Physical Surface("bot", 4) = {38};
Physical Surface("top", 5) = {96, 118, 140};
Physical Surface("side1", 6) = {149, 83, 127, 61, 105, 39};
Physical Surface("side2", 7) = {1, 4, 5, 2, 6, 3};


Physical Volume("fluid", 8) = {1, 2, 3, 4, 5, 6};//+
