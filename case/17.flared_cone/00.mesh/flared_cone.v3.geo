Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 0;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 0.001; // mm -> m
Geometry.Color.Surfaces = White;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)

// ---------------------
// *** nose section ***
// ---------------------
r_nose = 0.16;
r_flare = 3000.0;
cone_angle = 7.0;   //deg
flare_angle = 33.0; //deg
L = 470.0;

nx_0 = 80;
nx_1 = 800;
nx_2 = 400;
nx_3 = 800;
nx_4 = 800;
nx_5 = 800;
//nx_4 = 200;
//nx_5 = 30;
//nx_6 = 20;

ny = 160;

pNoseCent = newp;
Point(pNoseCent) = {0.0, 0.0, 0.0, lc};

y0 = 0.08;
x0 = -Sqrt(r_nose^2 - y0^2);
pList1[0] = newp;
Point(pList1[0]) = {x0 , y0 , 0.0 , lc};

x1 = 0.0;
y1 = r_nose;
pList1[1] = newp;
Point(pList1[1]) = {x1 , y1 , 0.0 , lc};

pFlareCent = newp;
Point(pFlareCent) = {0.0, r_flare+r_nose, 0.0, lc};

x2 = (L-r_nose)*1.0/8.0;
y2 = r_flare + r_nose - Sqrt(r_flare^2 - x2^2);
//y2 = 0.0;
pList1[2] = newp;
Point(pList1[2]) = {x2 , y2 , 0.0 , lc};

x3 = (L-r_nose)*1.0/4.0;
y3 = r_flare + r_nose - Sqrt(r_flare^2 - x3^2);
//y2 = 0.0;
pList1[3] = newp;
Point(pList1[3]) = {x3 , y3 , 0.0 , lc};

x4 = (L-r_nose)*2.0/4.0;
y4 = r_flare + r_nose - Sqrt(r_flare^2 - x4^2);
//y2 = 0.0;
pList1[4] = newp;
Point(pList1[4]) = {x4 , y4 , 0.0 , lc};

x5 = (L-r_nose)*3.0/4.0;
y5 = r_flare + r_nose - Sqrt(r_flare^2 - x5^2);
//y2 = 0.0;
pList1[5] = newp;
Point(pList1[5]) = {x5 , y5 , 0.0 , lc};

x6 = (L-r_nose)*4.0/4.0;
y6 = r_flare + r_nose - Sqrt(r_flare^2 - x6^2);
//y2 = 0.0;
pList1[6] = newp;
Point(pList1[6]) = {x6 , y6 , 0.0 , lc};



cone_lines[0] = newl;
Circle(cone_lines[0]) = {pList1[0], pNoseCent, pList1[1]}; Transfinite Line {cone_lines[0]} = nx_0;
cone_lines[1] = newl;
Circle(cone_lines[1]) = {pList1[1], pFlareCent, pList1[2]}; Transfinite Line {cone_lines[1]} = nx_1 Using Progression 1.0055;
cone_lines[2] = newl;
Circle(cone_lines[2]) = {pList1[2], pFlareCent, pList1[3]}; Transfinite Line {cone_lines[2]} = nx_2 Using Progression 1.000;
cone_lines[3] = newl;
Circle(cone_lines[3]) = {pList1[3], pFlareCent, pList1[4]}; Transfinite Line {cone_lines[3]} = nx_3 Using Progression 1.000;
cone_lines[4] = newl;
Circle(cone_lines[4]) = {pList1[4], pFlareCent, pList1[5]}; Transfinite Line {cone_lines[4]} = nx_4 Using Progression 1.000;
cone_lines[5] = newl;
Circle(cone_lines[5]) = {pList1[5], pFlareCent, pList1[6]}; Transfinite Line {cone_lines[5]} = nx_5 Using Progression 1.000;


// ------------------
// *** inlet line ***
// ------------------
//// inlet lines
yin0 = y0;
xin0 = -r_nose*2 ;
pList_in1[0] = newp;
Point(pList_in1[0]) = {xin0 , yin0 , 0.0 , lc};

pList_in1[1] = newp;
Point(pList_in1[1]) = {x1-0.05 , y1*3.5 , 0.0 , lc};

pList_in1[2] = newp;
Point(pList_in1[2]) = {x2-4 , y2*30 , 0.0 , lc};

pList_in1[3] = newp;
Point(pList_in1[3]) = {x3-4 , y3*10 , 0.0 , lc};

pList_in1[4] = newp;
Point(pList_in1[4]) = {x4-6, y4*3.6 , 0.0 , lc};

pList_in1[5] = newp;
Point(pList_in1[5]) = {x5-7, y5*2.3 , 0.0 , lc};

pList_in1[6] = newp;
Point(pList_in1[6]) = {x6-7, y6*1.8 , 0.0 , lc};


// *** inlet line 1 ***
pList_inner[0] = newp;
Point(pList_inner[0]) = {-0.32, 0.25, 0.0, 1.0};

INLET_LINES[0] = newl;
Bezier(newl) = {pList_in1[0], pList_inner[0], pList_in1[1]}; Transfinite Line {INLET_LINES[0]} = nx_0;

// *** inlet line 2 ***
pList_inner[1] = newp;
Point(pList_inner[1]) = {2, 3, 0.0, 1.0};

pList_inner[2] = newp;
Point(pList_inner[2]) = {20, 20, 0.0, 1.0};

INLET_LINES[1] = newl;
Bezier(newl) = {pList_in1[1], pList_inner[1], pList_inner[2], pList_in1[2]}; Transfinite Line {INLET_LINES[1]} = nx_1 Using Progression 1.005;


// *** inlet line 3 ***
INLET_LINES[2] = newl;
Bezier(newl) = {pList_in1[2], pList_in1[3]}; Transfinite Line {INLET_LINES[2]} = nx_2;

// *** inlet line 4 ***
INLET_LINES[3] = newl;
Bezier(newl) = {pList_in1[3], pList_in1[4]}; Transfinite Line {INLET_LINES[3]} = nx_3;

// *** inlet line 5 ***
INLET_LINES[4] = newl;
Bezier(newl) = {pList_in1[4], pList_in1[5]}; Transfinite Line {INLET_LINES[4]} = nx_4;

// *** inlet line 5 ***
INLET_LINES[5] = newl;
Bezier(newl) = {pList_in1[5], pList_in1[6]}; Transfinite Line {INLET_LINES[5]} = nx_5;


//////////span lines
span_lines[0] = newl;
Line(span_lines[0]) = {pList1[0], pList_in1[0]}; Transfinite Line {span_lines[0]} = ny Using Progression 1.00;
span_lines[1] = newl;
Line(span_lines[1]) = {pList1[1], pList_in1[1]}; Transfinite Line {span_lines[1]} = ny Using Progression 1.01;
span_lines[2] = newl;
Line(span_lines[2]) = {pList1[2], pList_in1[2]}; Transfinite Line {span_lines[2]} = ny Using Progression 1.02;
span_lines[3] = newl;
Line(span_lines[3]) = {pList1[3], pList_in1[3]}; Transfinite Line {span_lines[3]} = ny Using Progression 1.02;
span_lines[4] = newl;
Line(span_lines[4]) = {pList1[4], pList_in1[4]}; Transfinite Line {span_lines[4]} = ny Using Progression 1.02;
span_lines[5] = newl;
Line(span_lines[5]) = {pList1[5], pList_in1[5]}; Transfinite Line {span_lines[5]} = ny Using Progression 1.02;
span_lines[6] = newl;
Line(span_lines[6]) = {pList1[6], pList_in1[6]}; Transfinite Line {span_lines[6]} = ny Using Progression 1.02;


For idom In {0 : 5}
	cl = cone_lines[idom];
	sl1 = span_lines[idom+1];
	il = INLET_LINES[idom];
	sl0 = span_lines[idom];

	surfs[idom] = news;
	Curve Loop(surfs[idom]) = {cl, sl1, -il, -sl0};
	Plane Surface(surfs[idom]) = {surfs[idom]};
	Transfinite Surface {surfs[idom]};
	Recombine Surface(surfs[idom]);
EndFor

Extrude {{1, 0, 0}, {0, 0, 0}, Pi/40} {
  Surface{surfs[0]}; 
  Surface{surfs[1]}; 
  Surface{surfs[2]}; 
  Surface{surfs[3]}; 
  Surface{surfs[4]}; 
  Surface{surfs[5]}; 
  Layers {1}; Recombine;
}
//Physical Surface("inlet", 1) = {26, 48};
//Physical Surface("outlet", 2) = {44};
//Physical Surface("wall", 3) = {40, 18};
//Physical Surface("perio1", 4) = {31, 53};
//Physical Surface("perio2", 5) = {9, 8};
//Physical Surface("slip", 6) = {30};
//Physical Volume("fluid", 7) = {1, 2};
//
//Point(92) = {-2.5, -0, 292.7, 1.0};
//Point(93) = {15.7, 13.3, 292.7, 1.0};
//Point(94) = {66.8, 16.2, 292.7, 1.0};
//Point(95) = {133, 19.3, 292.7, 1.0};
//Point(96) = {259.3, 31.6, 292.7, 1.0};
//Point(97) = {387.2, 47, 292.7, 1.0};
//Point(98) = {462.3, 62.4, 292.7, 1.0};
//Spline(44) = {92, 93, 94, 95, 96, 97, 98};
//+
//+
//+
Physical Surface("inlet", 1) = {42, 64, 86, 108, 130, 152};
//+
Physical Surface("outlet", 2) = {148};
//+
Physical Surface("wall", 3) = {144, 122, 100, 78, 56, 34};
//+
Physical Surface("axis", 4) = {46};
Physical Surface("side1", 5) = {25, 24, 23, 22, 21, 20};
Physical Surface("side2",6) = {69, 91, 113, 135, 157, 47};
//+
Physical Volume("fluid", 7) = {1, 2, 3, 4, 5, 6};
//+
Recursive Delete {
  Surface{157}; 
}
//+
Recursive Delete {
  Surface{157}; 
}
//+
Recursive Delete {
  Surface{157}; 
}
