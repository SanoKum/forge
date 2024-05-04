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

nx_0 = 25;
nx_1 = 5;
nx_2 = 1000;
nx_3 = 120;
nx_4 = 30;
nx_5 = 20;

ny = 80;

pNoseCent = newp;
Point(pNoseCent) = {r_nose, 0.0, 0.0, lc};

y0 = 0.1;
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

x4 = 1101.25;
y4 = (x4-x3)*Tan(cone_angle/360.*2*Pi) + y3; // target: 274.55/2
Printf("y4=%f",y4);

pList1[4] = newp;
Point(pList1[4]) = {x4 , y4 , 0.0 , lc};

x5 = 1601.30;
y5 = y4;
pList1[5] = newp;
Point(pList1[5]) = {x5 , y5 , 0.0 , lc};

x6 = 1663.76;
y6 = (x6-x5)*Tan(flare_angle/360.*2*Pi) + y5;
pList1[6] = newp;
Point(pList1[6]) = {x6 , y6 , 0.0 , lc};

x7 = 1700.00;
y7 = y6;
pList1[7] = newp;
Point(pList1[7]) = {x7 , y7 , 0.0 , lc};

cone_lines[0] = newl;
Circle(cone_lines[0]) = {pList1[0], pNoseCent, pList1[1]}; Transfinite Line {cone_lines[0]} = nx_0;
cone_lines[1] = newl;
Bezier(cone_lines[1]) = {pList1[1], pList1[2], pList1[3]}; Transfinite Line {cone_lines[1]} = nx_1;
cone_lines[2] = newl;
Line(cone_lines[2]) = {pList1[3], pList1[4]}; Transfinite Line {cone_lines[2]} = nx_2 Using Progression 1.0025;
cone_lines[3] = newl;
Line(cone_lines[3]) = {pList1[4], pList1[5]}; Transfinite Line {cone_lines[3]} = nx_3 Using Progression 0.995;
cone_lines[4] = newl;
Line(cone_lines[4]) = {pList1[5], pList1[6]}; Transfinite Line {cone_lines[4]} = nx_4 Using Bump 0.3;
cone_lines[5] = newl;
Line(cone_lines[5]) = {pList1[6], pList1[7]}; Transfinite Line {cone_lines[5]} = nx_5 Using Progression 1.2;


// inlet
r_le_in = 1.0;
x_le_in = 1.5;
pInletCent = newp;
Point(pInletCent) = {x_le_in, 0.0, 0.0, lc};
a = 4.0;
b = 7.0;

majorAxisPoint = newp;
Point(majorAxisPoint) = {x_le_in, b, 0.0, lc};

inlet_cone_angle = 14.0;

yin0 = 0.1;
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

xin4 = x4-25.0;
yin4 = (xin4-xin3)*Tan(inlet_cone_angle/360.*2*Pi) + yin3;
pList_in1[4] = newp;
Point(pList_in1[4]) = {xin4 , yin4 , 0.0 , lc};

xin5 = x5-25.0;
yin5 = (xin5-xin4)*Tan(inlet_cone_angle/360.*2*Pi) + yin4;
pList_in1[5] = newp;
Point(pList_in1[5]) = {xin5 , yin5 , 0.0 , lc};

xin6 = x6-25.0;
yin6 = (xin6-xin5)*Tan(inlet_cone_angle/360.*2*Pi) + yin5;
pList_in1[6] = newp;
Point(pList_in1[6]) = {xin6 , yin6 , 0.0 , lc};

xin7 = x7-25.0;
yin7 = (xin7-xin6)*Tan(inlet_cone_angle/360.*2*Pi) + yin6;
pList_in1[7] = newp;
Point(pList_in1[7]) = {xin7 , yin7 , 0.0 , lc};



// inlet lines
inlet_lines[0] = newl;
Ellipse(inlet_lines[0]) = {pList_in1[0], pInletCent, majorAxisPoint, pList_in1[1]}; Transfinite Line {inlet_lines[0]} = nx_0;
inlet_lines[1] = newl;
Ellipse(inlet_lines[1]) = {pList_in1[1], pInletCent, majorAxisPoint, pList_in1[3]}; Transfinite Line {inlet_lines[1]} = nx_1;
inlet_lines[2] = newl;
Line(inlet_lines[2]) = {pList_in1[3], pList_in1[4]}; Transfinite Line {inlet_lines[2]} = nx_2 Using Progression 1.0025;
inlet_lines[3] = newl;
Line(inlet_lines[3]) = {pList_in1[4], pList_in1[5]}; Transfinite Line {inlet_lines[3]} = nx_3 Using Progression 0.995;
inlet_lines[4] = newl;
Line(inlet_lines[4]) = {pList_in1[5], pList_in1[6]}; Transfinite Line {inlet_lines[4]} = nx_4 Using Bump 0.3;
inlet_lines[5] = newl;
Line(inlet_lines[5]) = {pList_in1[6], pList_in1[7]}; Transfinite Line {inlet_lines[5]} = nx_5 Using Progression 1.2;

////span lines
span_lines[0] = newl;
Line(span_lines[0]) = {pList1[0], pList_in1[0]}; Transfinite Line {span_lines[0]} = ny Using Progression 1.02;
span_lines[1] = newl;
Line(span_lines[1]) = {pList1[1], pList_in1[1]}; Transfinite Line {span_lines[1]} = ny Using Progression 1.02;
//span_lines[2] = newl;
//Line(span_lines[2]) = {pList1[2], pList_in1[2]}; Transfinite Line {span_lines[2]} = ny Using Progression 1.1;
span_lines[2] = newl;
Line(span_lines[2]) = {pList1[3], pList_in1[3]}; Transfinite Line {span_lines[2]} = ny Using Progression 1.05;
span_lines[3] = newl;
Line(span_lines[3]) = {pList1[4], pList_in1[4]}; Transfinite Line {span_lines[3]} = ny Using Progression 1.05;
span_lines[4] = newl;
Line(span_lines[4]) = {pList1[5], pList_in1[5]}; Transfinite Line {span_lines[4]} = ny Using Progression 1.05;
span_lines[5] = newl;
Line(span_lines[5]) = {pList1[6], pList_in1[6]}; Transfinite Line {span_lines[5]} = ny Using Progression 1.05;
span_lines[6] = newl;
Line(span_lines[6]) = {pList1[7], pList_in1[7]}; Transfinite Line {span_lines[6]} = ny Using Progression 1.05;

For idom In {0 : 5}
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

Extrude {{1, 0, 0}, {0, 0, 0}, Pi/16} {
  Surface{surfs[0]}; 
  Surface{surfs[1]}; 
  Surface{surfs[2]}; 
  Surface{surfs[3]}; 
  Surface{surfs[4]}; 
  Surface{surfs[5]}; 
  Layers {3}; Recombine;
}
Physical Surface("inlet", 1) = {42, 86, 64, 108, 130, 152};
Physical Surface("outlet", 2) = {148};
Physical Surface("wall", 3) = {144, 122, 100, 78, 34, 56};
Physical Surface("perio1", 4) = {47, 69, 91, 113, 135, 157};
Physical Surface("perio2", 5) = {25, 24, 23, 22, 21, 20};
Physical Surface("slip", 6) = {46};
Physical Volume("fluid", 7) = {1, 2, 3, 4, 5, 6};
//+
