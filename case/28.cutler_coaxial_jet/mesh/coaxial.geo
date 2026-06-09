// Cutler 超音速同軸ジェット (軸対称 2D, y=r)。Cutler et al. 2006 (J. Propulsion & Power) Fig.1 + 本文。
//   出口面 (x=0) 径方向寸法 (実機値):
//     中心ジェット出口 内径 10.00mm  -> Rj  = 5.00mm   (He/O2, Mach1.8)
//     center body lip 厚 0.25mm       -> tlip= 0.25mm   (lip 壁; 外縁 r=5.25mm)
//     coflow ノズル出口 内径 60.47mm  -> Ro  = 30.235mm (空気, Mach1.8)
//     周囲(静止)空気                   -> r > Ro, 遠方場 Rfar=60mm まで
//   計測は x=0.13mm から ~261mm まで及ぶため下流長 L=280mm。ScalingFactor=0.001 (mm->m)。
Mesh.ScalingFactor = 0.001;
SetFactory("OpenCASCADE");

Rj    = 5.0;       // 中心ジェット半径 [mm]
tlip  = 0.25;      // center body lip 厚 [mm] (Cutler2006 本文)
Ro    = 30.235;    // coflow 外半径 [mm] (= 60.47mm ID / 2)
Rfar  = 60.0;      // 遠方場半径 [mm]
L     = 280.0;     // 下流長 [mm] (計測 ~261mm をカバー)
xin   = -5.0;      // 入口を少し上流に (ノズル exit を x=0 に)

// y(=r) バンド境界点 (x=xin 入口面)
Point(1)={xin,0,0};      Point(2)={xin,Rj,0};
Point(3)={xin,Rj+tlip,0};Point(4)={xin,Ro,0};
Point(5)={xin,Rfar,0};
// 下流端 x=L
Point(6)={L,0,0};        Point(7)={L,Rj,0};
Point(8)={L,Rj+tlip,0};  Point(9)={L,Ro,0};
Point(10)={L,Rfar,0};

// 軸 (y=0)
Line(1)={1,6};
// 出口面 (x=L) 4 区間
Line(2)={6,7}; Line(3)={7,8}; Line(4)={8,9}; Line(5)={9,10};
// 上端 (y=Rfar)
Line(6)={10,5};
// 入口面 (x=xin) 4 区間 (上から下)
Line(7)={5,4}; Line(8)={4,3}; Line(9)={3,2}; Line(10)={2,1};
// 内部水平線 (バンド分割; lip を境界化するため)
Line(11)={2,7};   // r=Rj
Line(12)={3,8};   // r=Rj+tlip
Line(13)={4,9};   // r=Ro

// 4 サブ領域 (center jet / lip / coflow / ambient)
Curve Loop(1)={1,2,-11,10};            Plane Surface(1)={1};  // center jet band
Curve Loop(2)={11,3,-12,9};            Plane Surface(2)={2};  // lip band
Curve Loop(3)={12,4,-13,8};            Plane Surface(3)={3};  // coflow band
Curve Loop(4)={13,5,6,7};              Plane Surface(4)={4};  // ambient band

// 構造格子。coflow バンドが ~25mm と広いので径方向セルを増やす。
// 剪断層 (r=Rj 近傍) はバンド境界で自然に細かい side に揃う。lip は薄いので少数。
nx=320; ncj=64; nlip=3; ncof=120; namb=56;
Transfinite Curve{1,11,12,13,6}=nx Using Progression 1.0;
Transfinite Curve{10,-2}=ncj;  Transfinite Curve{9,-3}=nlip;
Transfinite Curve{8,-4}=ncof;  Transfinite Curve{7,-5}=namb;
Transfinite Surface{1,2,3,4}; Recombine Surface{1,2,3,4};

// Physical (boundary)
Physical Curve("axis",1)   = {1};
Physical Curve("jet_in",2) = {10};               // center jet inlet
Physical Curve("lip",3)    = {9};                // lip wall (inlet face of lip band)
Physical Curve("coflow_in",4)={8};               // coflow inlet
Physical Curve("amb_in",5) = {7};                // ambient inlet/farfield (left-top)
Physical Curve("top",6)    = {6};                // far-field top
Physical Curve("outflow",7)= {2,3,4,5};          // outflow (x=L)
Physical Surface("fluid",8)= {1,2,3,4};
