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
r_nose = 2.5;
cone_angle = 7.0;   //deg
flare_angle = 33.0; //deg

nx_0 = 100;
nx_1 = 5;
nx_2 = 150;
nx_3 = 500;
nx_4 = 200;
nx_5 = 30;
nx_6 = 20;

ny = 120;

pNoseCent = newp;
Point(pNoseCent) = {r_nose, 0.0, 0.0, lc};

y0 = 0.2;
x0 = 2.5-Sqrt(r_nose^2 - y0^2);
pList1[0] = newp;
Point(pList1[0]) = {x0 , y0 , 0.0 , lc};

x1 = 1.15;
dx = r_nose-x1;
y1 = Sqrt(r_nose*r_nose - dx*dx);
pList1[1] = newp;
Point(pList1[1]) = {x1 , y1 , 0.0 , lc};

x2 = 1.29989;
dx = r_nose-x2;
y2 = Sqrt(r_nose*r_nose - dx*dx);
pList1[2] = newp;
Point(pList1[2]) = {x2 , y2 , 0.0 , lc};

x3 = 1.5;
y3 = (x3-x2)*Tan(cone_angle/360.*2*Pi) + y2; 
pList1[3] = newp;
Point(pList1[3]) = {x3 , y3 , 0.0 , lc};

x4 = 60;
y4 = (x4-x3)*Tan(cone_angle/360.*2*Pi) + y3; 
pList1[4] = newp;
Point(pList1[4]) = {x4 , y4 , 0.0 , lc};

x5 = 1101.25;
y5 = (x5-x4)*Tan(cone_angle/360.*2*Pi) + y4; // target: 274.55/2

pList1[5] = newp;
Point(pList1[5]) = {x5 , y5 , 0.0 , lc};

x6 = 1601.30;
y6 = y5;
pList1[6] = newp;
Point(pList1[6]) = {x6 , y6 , 0.0 , lc};

x7 = 1663.76;
y7 = (x7-x6)*Tan(flare_angle/370.*2*Pi) + y6;
pList1[7] = newp;
Point(pList1[7]) = {x7 , y7 , 0.0 , lc};

x8 = 1700.00;
y8 = y7;
pList1[8] = newp;
Point(pList1[8]) = {x8 , y8 , 0.0 , lc};


cone_lines[0] = newl;
Circle(cone_lines[0]) = {pList1[0], pNoseCent, pList1[1]}; Transfinite Line {cone_lines[0]} = nx_0;
cone_lines[1] = newl;
Bezier(cone_lines[1]) = {pList1[1], pList1[2], pList1[3]}; Transfinite Line {cone_lines[1]} = nx_1;
cone_lines[2] = newl;
Line(cone_lines[2]) = {pList1[3], pList1[4]}; Transfinite Line {cone_lines[2]} = nx_2 Using Progression 1.012;
cone_lines[3] = newl;
Line(cone_lines[3]) = {pList1[4], pList1[5]}; Transfinite Line {cone_lines[3]} = nx_3 Using Bump 0.25;
cone_lines[4] = newl;
Line(cone_lines[4]) = {pList1[5], pList1[6]}; Transfinite Line {cone_lines[4]} = nx_4 Using Bump 0.35;
cone_lines[5] = newl;
Line(cone_lines[5]) = {pList1[6], pList1[7]}; Transfinite Line {cone_lines[5]} = nx_5 Using Bump 0.3;
cone_lines[6] = newl;
Line(cone_lines[6]) = {pList1[7], pList1[8]}; Transfinite Line {cone_lines[6]} = nx_6 Using Progression 1.08;



// inlet
r_le_in = 1.0;
x_le_in = 2.5;
pInletCent = newp;
Point(pInletCent) = {x_le_in, 0.0, 0.0, lc};
a = 5.0;
b = 8.0;

majorAxisPoint = newp;
Point(majorAxisPoint) = {x_le_in, b, 0.0, lc};

inlet_cone_angle = 14.0;

yin0 = y0;
xin0 = x_le_in - a*Sqrt(r_le_in^2 -(yin0*a/b)^2 ) ;
pList_in1[0] = newp;
Point(pList_in1[0]) = {xin0 , yin0 , 0.0 , lc};

xin1 = -0.4;
yin1 = b*Sqrt(r_le_in^2 -((xin1-x_le_in)/a)^2 ) ;
//yin1 = Sqrt(r_le_in^2-(xin1-x_le_in)^2);
pList_in1[1] = newp;
Point(pList_in1[1]) = {xin1 , yin1 , 0.0 , lc};

xin2 = -0.1;
yin2 = b*Sqrt(r_le_in^2 -((xin2-x_le_in)/a)^2 ) ;
pList_in1[2] = newp;
Point(pList_in1[2]) = {xin2 , yin2 , 0.0 , lc};

xin3 = 0.3;
yin3 = b*Sqrt(r_le_in^2 -((xin3-x_le_in)/a)^2 ) ;
pList_in1[3] = newp;
Point(pList_in1[3]) = {xin3 , yin3 , 0.0 , lc};

xin3_1 = 4 ;
yin3_1 = 10.5 ;
dummy_p1 = newp;
Point(dummy_p1) = {xin3_1 , yin3_1 , 0.0 , lc};

xin3_2 = 25.0;
yin3_2 = 22.0 ;
dummy_p2 = newp;
Point(dummy_p2) = {xin3_2 , yin3_2 , 0.0 , lc};

xin4 = x4-5.0;
yin4 = 30.0;
pList_in1[4] = newp;
Point(pList_in1[4]) = {xin4 , yin4 , 0.0 , lc};

xin5 = x5-15.0;
yin5 = (xin5-xin3)*Tan(inlet_cone_angle/360.*2*Pi) + yin3;
pList_in1[5] = newp;
Point(pList_in1[5]) = {xin5 , yin5 , 0.0 , lc};


xin6 = x6-30.0;
yin6 = (xin6-xin4)*Tan(inlet_cone_angle/360.*2*Pi) + yin4;
pList_in1[6] = newp;
Point(pList_in1[6]) = {xin6 , yin6 , 0.0 , lc};

xin7 = x7-30.0;
yin7 = (xin7-xin6)*Tan(inlet_cone_angle/360.*2*Pi) + yin6;
pList_in1[7] = newp;
Point(pList_in1[7]) = {xin7 , yin7 , 0.0 , lc};

xin8 = x8-30.0;
yin8 = (xin8-xin7)*Tan(inlet_cone_angle/360.*2*Pi) + yin7;
pList_in1[8] = newp;
Point(pList_in1[8]) = {xin8 , yin8 , 0.0 , lc};



// inlet lines
inlet_lines[0] = newl;
Ellipse(inlet_lines[0]) = {pList_in1[0], pInletCent, majorAxisPoint, pList_in1[1]}; Transfinite Line {inlet_lines[0]} = nx_0;
inlet_lines[1] = newl;
Ellipse(inlet_lines[1]) = {pList_in1[1], pInletCent, majorAxisPoint, pList_in1[3]}; Transfinite Line {inlet_lines[1]} = nx_1;
inlet_lines[2] = newl;
//Line(inlet_lines[2]) = {pList_in1[3], pList_in1[4]}; Transfinite Line {inlet_lines[2]} = nx_2 Using Progression 1.0025;
Bezier(inlet_lines[2]) = {pList_in1[3], dummy_p1, dummy_p2, pList_in1[4]}; Transfinite Line {inlet_lines[2]} = nx_2 Using Progression 1.009;
inlet_lines[3] = newl;
Line(inlet_lines[3]) = {pList_in1[4], pList_in1[5]}; Transfinite Line {inlet_lines[3]} = nx_3 Using Bump 0.25;
inlet_lines[4] = newl;
Line(inlet_lines[4]) = {pList_in1[5], pList_in1[6]}; Transfinite Line {inlet_lines[4]} = nx_4 Using Bump 0.35;
inlet_lines[5] = newl;
Line(inlet_lines[5]) = {pList_in1[6], pList_in1[7]}; Transfinite Line {inlet_lines[5]} = nx_5 Using Bump 0.3;
inlet_lines[6] = newl;
Line(inlet_lines[6]) = {pList_in1[7], pList_in1[8]}; Transfinite Line {inlet_lines[6]} = nx_6 Using Progression 1.08;



//////span lines
span_lines[0] = newl;
Line(span_lines[0]) = {pList1[0], pList_in1[0]}; Transfinite Line {span_lines[0]} = ny Using Progression 1.02;
span_lines[1] = newl;
Line(span_lines[1]) = {pList1[1], pList_in1[1]}; Transfinite Line {span_lines[1]} = ny Using Progression 1.02;
//span_lines[2] = newl;
//Line(span_lines[2]) = {pList1[2], pList_in1[2]}; Transfinite Line {span_lines[2]} = ny Using Progression 1.1;
span_lines[2] = newl;
Line(span_lines[2]) = {pList1[3], pList_in1[3]}; Transfinite Line {span_lines[2]} = ny Using Progression 1.02;
span_lines[3] = newl;
Line(span_lines[3]) = {pList1[4], pList_in1[4]}; Transfinite Line {span_lines[3]} = ny Using Progression 1.07;
span_lines[4] = newl;
Line(span_lines[4]) = {pList1[5], pList_in1[5]}; Transfinite Line {span_lines[4]} = ny Using Progression 1.07;
span_lines[5] = newl;
Line(span_lines[5]) = {pList1[6], pList_in1[6]}; Transfinite Line {span_lines[5]} = ny Using Progression 1.07;
span_lines[6] = newl;
Line(span_lines[6]) = {pList1[7], pList_in1[7]}; Transfinite Line {span_lines[6]} = ny Using Progression 1.07;
span_lines[7] = newl;
Line(span_lines[7]) = {pList1[8], pList_in1[8]}; Transfinite Line {span_lines[7]} = ny Using Progression 1.07;


For idom In {0 : 6}
	cl = cone_lines[idom];
	sl1 = span_lines[idom+1];
	il = inlet_lines[idom];
	sl0 = span_lines[idom];

	surfs[idom] = news;
	Curve Loop(surfs[idom]) = {cl, sl1, -il, -sl0};
	Plane Surface(surfs[idom]) = {surfs[idom]};
	Transfinite Surface {surfs[idom]};
	Recombine Surface(surfs[idom]);
EndFor

Extrude {{1, 0, 0}, {0, 0, 0}, Pi/32} {
  Surface{surfs[0]}; 
  Surface{surfs[1]}; 
  Surface{surfs[2]}; 
  Surface{surfs[3]}; 
  Surface{surfs[4]}; 
  Surface{surfs[5]}; 
  Surface{surfs[6]}; 
  Layers {1}; Recombine;
}
//temp Physical Surface("inlet", 1) = {42, 86, 64, 108, 130, 152};
//temp Physical Surface("outlet", 2) = {148};
//temp Physical Surface("wall", 3) = {144, 122, 100, 78, 34, 56};
//temp Physical Surface("perio1", 4) = {47, 69, 91, 113, 135, 157};
//temp Physical Surface("perio2", 5) = {25, 24, 23, 22, 21, 20};
//temp Physical Surface("slip", 6) = {46};
//temp Physical Volume("fluid", 7) = {1, 2, 3, 4, 5, 6};
//temp //+
//temp //+
//+
//+
//+
Physical Surface("inlet", 1) = {46, 68, 90, 112, 134, 156, 178};
Physical Surface("outlet", 2) = {174};
Physical Surface("wall", 3) = {170, 148, 126, 104, 82, 38, 60};
Physical Surface("perio1", 4) = {183, 161, 139, 117, 95, 73, 51};
Physical Surface("perio2", 5) = {23, 24, 25, 26, 27, 28, 29};
Physical Surface("slip", 6) = {50};
Physical Volume("fluid", 7) = {1, 2, 3, 4, 5, 6, 7};
