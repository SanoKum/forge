# TGV Re=1600 DNS 参照データ

`TGV_Re1600.dat`: Incompact3d 512³ DNS (非圧縮, (2π)³, ν=1/1600) の統計時系列 (t=0..20, 1001点)。
列: 1=t, 2=Ek (Ek(0)=0.125), 3=ε_t=−dEk/dt, 4=ε(勾配ベース), 5=enstrophy, 6-17=各成分二乗平均。

出典 (引用義務): Dairay, Lamballais, Laizet, Vassilicos,
"Numerical dissipation vs. subgrid-scale modelling for large eddy simulation",
JCP 337 (2017) 252-274. https://doi.org/10.1016/j.jcp.2017.02.035
取得元: https://github.com/xcompact3d/Incompact3d/blob/master/examples/TGV-Taylor-Green-vortex/TGV_Re1600.dat

forge との対応: forge は M0.4 圧縮性 (V0=0.4, c=1, ro0=1, L=1)。t* = t·V0/L が DNS の t に対応。
K/K0 = Ek/0.125、ε* = −d(K/K0)/dt* = ε_t/0.125 (ピーク ≈ 0.0908 @ t≈9.0, K/K0(10)=0.596)。
