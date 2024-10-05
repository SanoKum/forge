Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 0.01; // cm -> m
//Geometry.Color.Surfaces = Black;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)


// -------------------
// *** (1) section ***
// -------------------
x0 = -6.363;
x1 = -4.277;
nx_1 = 70;
dx = (x1-x0);

For i1 In {0 : 1}
	x = x0 + dx*i1;
	y = 1.2450 + 0.00141*x;

	pList1[i1] = newp;

	Point(pList1[i1]) = {x , y , 0.0 , lc};
EndFor

l1= newl;
Spline(l1) =  pList1[]; Transfinite Line {l1} = nx_1;


// -------------------
// *** (2) section ***
// -------------------
x2 = -2.399;
nx_2 = 70;
dx = (x2-x1);

pList2[0] = pList1[1];

//For i2 In {1 : 1}
	x = x1 + dx;
	y = -0.3059 -0.3612*x;

	pList2[1] = newp;

	Point(pList2[1]) = {x , y , 0.0 , lc};
//EndFor

l2 = newl;
Spline(l2) =  pList2[]; Transfinite Line {l2} = nx_2;


// -------------------
// *** (3) section ***
// -------------------
x3 = 0.595;
nx_3 = 70;
dx = (x3-x2)/nx_3;

pList3[0] = pList2[1];

For i3 In {1 : nx_3}
	x = x2 + dx*i3;
	y = 0.2255 + 0.00141*x + 0.025457*x*x -0.013907*x*x*x;

	pList3[i3] = newp;

	Point(pList3[i3]) = {x , y , 0.0 , lc};
EndFor

l3 = newl;
Spline(l3) =  pList3[]; Transfinite Line {l3} = nx_3;


// -------------------
// *** (4) section ***
// -------------------
x4 = 9.45;
nx_4 = 250;
dx = (x4-x3);

pList4[0] = pList3[nx_3];

//For i4 In {1 : nx_4}
	x = x3 + dx;
	y = 0.2224+0.01689*x;

	pList4[1] = newp;

	Point(pList4[1]) = {x , y , 0.0 , lc};
//EndFor

l4 = newl;
Spline(l4) =  pList4[]; Transfinite Line {l4} = nx_4;

// ----------------------
// *** outlet section ***
// ----------------------
nx_5 = 80;
pList5[0] = pList4[1];
pList5[1] = newp;
x = x4;
y = -0.2224 + -0.01689*x;
Point(pList5[1]) = {x , y, 0.0 , lc};

l5 = newl;
Line(l5) = {pList5[0], pList5[1]}; Transfinite Line {l5} = nx_5 Using Bump 0.02;


// -------------------
// *** (-4) section ***
// -------------------
dx = (x4-x3);
pList_m4[0] = pList5[1];

//For i4 In {1 : nx_4}
	x = x4 - dx;
	y = -0.2224 -0.01689*x;

	pList_m4[1] = newp;

	Point(pList_m4[1]) = {x , y , 0.0 , lc};
//EndFor

l_m4 = newl;
Spline(l_m4) =  pList_m4[]; Transfinite Line {l_m4} = nx_4;

// --------------------
// *** (-3) section ***
// --------------------
dx = (x3-x2)/nx_3;

pList_m3[0] = pList_m4[1];

For i3 In {1 : nx_3}
	x = x3 - dx*i3;
	y = -0.2255 + -0.00141*x + -0.025457*x*x +0.013907*x*x*x;

	pList_m3[i3] = newp;

	Point(pList_m3[i3]) = {x , y , 0.0 , lc};
EndFor

l_m3 = newl;
Spline(l_m3) =  pList_m3[]; Transfinite Line {l_m3} = nx_3;




// -------------------
// *** (-2) section ***
// -------------------
dx = (x2-x1);

pList_m2[0] = pList_m3[nx_3];

//For i2 In {1 : 1}
	x = x2 - dx;
	y = +0.3059 +0.3612*x;

	pList_m2[1] = newp;

	Point(pList_m2[1]) = {x , y , 0.0 , lc};
//EndFor

l_m2 = newl;
Spline(l_m2) =  pList_m2[]; Transfinite Line {l_m2} = nx_2;

// -------------------
// *** (-1) section ***
// -------------------
dx = (x1-x0);

pList_m1[0] = pList_m2[1];
//For i1 In {0 : 1}
	x = x1 - dx;
	y = -1.2450 - 0.00141*x;

	pList_m1[1] = newp;

	Point(pList_m1[1]) = {x , y , 0.0 , lc};
//EndFor

l_m1= newl;
Spline(l_m1) =  pList_m1[]; Transfinite Line {l_m1} = nx_1;

l_in = newl;
Line(l_in) = {pList_m1[1], pList1[0]}; Transfinite Line {l_in} = nx_5 Using Bump 0.02;

l_mid1 = newl;
Line(l_mid1) = {pList1[1], pList_m1[0]}; Transfinite Line {l_mid1} = nx_5 Using Bump 0.02;

l_mid2 = newl;
Line(l_mid2) = {pList2[1], pList_m2[0]}; Transfinite Line {l_mid2} = nx_5 Using Bump 0.02;

l_mid3 = newl;
Line(l_mid3) = {pList3[nx_3], pList_m3[0]}; Transfinite Line {l_mid3} = nx_5 Using Bump 0.02;

Curve Loop(1) = {10, 1, 11, 9};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface(1);

Curve Loop(2) = {11, -8, -12, -2};
Plane Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface(2);


Curve Loop(3) = {12, -7, -13, -3};
Plane Surface(3) = {3};
Transfinite Surface {3};
Recombine Surface(3);


Curve Loop(4) = {13, -6, -5, -4};
Plane Surface(4) = {4};
Transfinite Surface {4};
Recombine Surface(4);

Extrude {0, 0, 1.259} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Layers {5}; Recombine;
}

Physical Surface("inlet", 1) = {22};
Physical Surface("outlet", 2) = {96};
Physical Surface("wall", 3) = {1, 26, 35, 34, 48, 2, 56, 57, 70, 3, 78, 79, 4, 92, 100, 101};
Physical Volume("fluid", 4) = {1, 2, 3, 4};
