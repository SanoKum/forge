Geometry.SurfaceNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.LineNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
Mesh.ScalingFactor = 0.001;
//Geometry.Color.Surfaces = Black;

lc=0.1; //characteristic mesh size (optional make smaller to refine mesh)

rcr  = 8.5;   // mm
rin  = 136.0/2; // mm
rout = 42.0/2; // mm

rcyl  = 0.1; // mm


Throat  = 17.0;

nx_strght = 85*2;
L_strght = 100.0;

nx_conv = 100*2;
L_conv  = 238.0;

nx_conv = 100*2;
L_conv  = 238.0;

nx_div = 100*2;
L_div = 238.0;


nx_inl = 20*2;
nx_out = 20*2;



// ------------------------
// *** Straight section ***
// ------------------------

Point(1) = {-L_strght , rin ,0.0, lc};
Point(2) = {0.0,rin,0.0, lc};
Line(1) = {1,2}; Transfinite Line {1} = nx_strght;

// --------------------------
// *** Convergent section ***
// --------------------------
pListC[0] = 2;

For iconv In {1 : nx_conv}
	x = 0.0 + L_conv*iconv/(nx_conv);
	y = rcr/Sqrt(1.0-(1.0-(rcr/rin)^2)*(1-(x/L_conv)^2)^2/(1.0+((x/L_conv)^2)/3.0)^3);

	pListC[iconv] = newp;

	Point(pListC[iconv]) = {x , y , 0.0 , lc};
EndFor

lconv = newl;
Spline(lconv) =  pListC[]; Transfinite Line {lconv} = nx_conv;

// -------------------------
// *** Divergent section ***
// -------------------------
pListD[0] = pListC[nx_conv];
pListD[1] = newp;
x = L_conv + L_div;
y = rout;
Point(pListD[1]) = {x , y , 0.0 , lc};

ldiv = newl;
Line(ldiv) = {pListD[0], pListD[1]}; Transfinite Line {ldiv} = nx_div;

// ----------------------
// *** outlet section ***
// ----------------------
pListO[0] = pListD[1];
pListO[1] = newp;
x = L_conv + L_div;
y = -rout;
//Point(pListO[1]) = {x , y , 0.0 , lc};
Point(pListO[1]) = {x , rcyl , 0.0 , lc};

lout = newl;
Line(lout) = {pListO[0], pListO[1]}; Transfinite Line {lout} = nx_out Using Bump 0.1;

// --------------------------
// *** Divergent section2 ***
// --------------------------
pListD2[0] = pListO[1];
pListD2[1] = newp;
x = L_conv ;
y = -17.0/2;
//Point(pListD2[1]) = {x , y , 0.0 , lc};
Point(pListD2[1]) = {x , rcyl , 0.0 , lc};

ldiv2 = newl;
Line(ldiv2) = {pListD2[0], pListD2[1]}; Transfinite Line {ldiv2} = nx_div;

// --------------------------
// *** Convergent section2 ***
// --------------------------
pListC2[0] = pListD2[1];

For iconv In {1 : nx_conv}
	x = 0.0 + L_conv*(nx_conv-iconv)/(nx_conv);
	y = -rcr/Sqrt(1.0-(1.0-(rcr/rin)^2)*(1-(x/L_conv)^2)^2/(1.0+((x/L_conv)^2)/3.0)^3);

	pListC2[iconv] = newp;

	//Point(pListC2[iconv]) = {x , y , 0.0 , lc};
	Point(pListC2[iconv]) = {x , rcyl , 0.0 , lc};
EndFor

lconv2 = newl;
Spline(lconv2) =  pListC2[]; Transfinite Line {lconv2} = nx_conv;

// -------------------------
// *** Straight section2 ***
// -------------------------

pListS2[0] = pListC2[nx_conv];
pListS2[1] = newp;

//Point(pListS2[1]) = {-L_strght,-rin,0.0, lc};
Point(pListS2[1]) = {-L_strght,rcyl,0.0, lc};
lst2 = newl;
Line(lst2) = {pListS2[0],pListS2[1]}; Transfinite Line {lst2} = nx_strght;


// -------------
// *** Inlet ***
// -------------
linl = newl;
Line(linl) = {pListS2[1],1}; Transfinite Curve {linl} = nx_inl Using Bump 0.1;

// mid
lmid1 = newl;
Line(lmid1) = {pListC[0],pListS2[0]}; Transfinite Curve {lmid1} = nx_inl Using Bump 0.1;

// mid2
lmid2 = newl;
Line(lmid2) = {pListD[0],pListC2[0]}; Transfinite Curve {lmid2} = nx_inl Using Bump 0.1;

Line Loop(1) = {1,lmid1,lst2,linl};
Plane Surface(1) = {1}; //using LL #1
Transfinite Surface {1};
Recombine Surface(1);

Line Loop(2) = {lconv,lmid2,lconv2,-lmid1};
Plane Surface(2) = {2}; //using LL #1
Transfinite Surface {2};
Recombine Surface(2);

Line Loop(3) = {ldiv,lout,ldiv2,-lmid2};
Plane Surface(3) = {3}; //using LL #1
Transfinite Surface {3};
Recombine Surface(3);

//
//Extrude {0, 0, 5} {
//  Surface{1,2,3}; Layers{1}; Recombine;
//}
//
//Extrude {0, 0, 5} {
//  Surface{2}; Layers{2}; Recombine;
//}




//Line Loop(1) = {1,lconv,ldiv,lout,ldiv2,lconv2,lst2,linl};



Extrude {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Surface{1,2,3}; Layers{5}; Recombine;
}

//Physical Surface("inlet", 1) = {26};
//Physical Surface("outlet", 2) = {56};
//Physical Surface("wall", 3) = {19, 36, 53};
//Physical Surface("slip", 4) = {27, 1, 44, 2, 61, 3};
//Physical Volume("fluid", 5) = {1, 2, 3};
//+
Physical Surface("inlet", 1) = {31};
Physical Surface("outlet", 2) = {67};
Physical Surface("wall", 3) = {19, 41, 63};
Physical Surface("bot", 4) = {27, 49, 71};
Physical Surface("per1", 5) = {32, 54, 76};
Physical Surface("per2", 6) = {1, 2, 3};
Physical Volume("fluid", 7) = {1, 2, 3};
