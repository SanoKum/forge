#pragma once
// 一般EOS flux Jacobian の固有系 (Method A: 数値 L=R^-1)。
// plan time_integration-general-eos-jacobian.md / docs/time_integration。
// CPG では h=c²/(γ-1)・χ=0 で従来の閉形式に簡約されるが、TP (NASA) では別物。
// 法線方向 n の Euler 流束 Jacobian A_n = R Λ L を、実 H_t・κ・χ から構築する。
//   音響:  r∓ = [1, u∓c n, H_t ∓ c U_n]
//   接触:  r_c = [1, u, H_t - c²/κ]          (CPG で H_t-c²/κ = K)
//   せん断: r_sk = [0, t_k, u·t_k]
// L=R^-1 を部分ピボット付き double 5×5 で数値的に解き、A± = R Λ± L を返す。

#include "cuda_forge/thermo_d.cuh"  // THERMO_HD

// 5×5 の逆行列 (Gauss-Jordan, 部分ピボット, double)。
// 戻り値 = 最小ピボット絶対値 (≈0 で特異)。Rinv に R^-1 を書く。
THERMO_HD double eos_inv5(const double Rin[5][5], double Rinv[5][5])
{
    double A[5][10];
    for (int i=0;i<5;i++){
        for (int j=0;j<5;j++) A[i][j]   = Rin[i][j];
        for (int j=0;j<5;j++) A[i][5+j] = (i==j) ? 1.0 : 0.0;
    }
    double minpiv = 1.0e300;
    for (int col=0; col<5; ++col){
        // 部分ピボット: 列 col で |A[r][col]| 最大の行を選ぶ
        int piv = col; double pmax = fabs(A[col][col]);
        for (int r=col+1; r<5; ++r){ double v=fabs(A[r][col]); if (v>pmax){pmax=v; piv=r;} }
        if (piv != col){ for (int j=0;j<10;j++){ double t=A[col][j]; A[col][j]=A[piv][j]; A[piv][j]=t; } }
        const double d = A[col][col];
        if (fabs(d) < minpiv) minpiv = fabs(d);
        if (fabs(d) < 1.0e-300) return 0.0;  // 特異
        // 他行を消去
        for (int r=0; r<5; ++r){
            if (r==col) continue;
            const double f = A[r][col]/d;
            if (f != 0.0) for (int j=col;j<10;j++) A[r][j] -= f*A[col][j];
        }
    }
    for (int i=0;i<5;i++){ const double d=A[i][i]; for (int j=0;j<5;j++) Rinv[i][j]=A[i][5+j]/d; }
    return minpiv;
}

// n に直交する正規直交接線 t1,t2 を作る (python 検証と同一手順)。
THERMO_HD void eos_tangents(double nx,double ny,double nz,
                            double t1[3], double t2[3])
{
    double ax,ay,az;
    if (fabs(nx) < 0.9){ ax=1.0; ay=0.0; az=0.0; } else { ax=0.0; ay=1.0; az=0.0; }
    // t1 = normalize(n × a)
    double c1x = ny*az - nz*ay;
    double c1y = nz*ax - nx*az;
    double c1z = nx*ay - ny*ax;
    double inv = 1.0/fmax(sqrt(c1x*c1x+c1y*c1y+c1z*c1z), 1.0e-300);
    t1[0]=c1x*inv; t1[1]=c1y*inv; t1[2]=c1z*inv;
    // t2 = n × t1
    t2[0] = ny*t1[2] - nz*t1[1];
    t2[1] = nz*t1[0] - nx*t1[2];
    t2[2] = nx*t1[1] - ny*t1[0];
}

// 一般EOS右固有ベクトル行列 R (列=固有ベクトル) と固有値 lambda を作る。
// 入力: 速度 (ux,uy,uz)、単位法線 (nx,ny,nz)、音速 c、全エンタルピー Ht、κ。
THERMO_HD void eos_eigvecs_general(double ux,double uy,double uz,
                                   double nx,double ny,double nz,
                                   double c, double Ht, double kappa,
                                   double R[5][5], double lambda[5])
{
    const double Un = ux*nx + uy*ny + uz*nz;
    double t1[3], t2[3]; eos_tangents(nx,ny,nz,t1,t2);
    const double ut1 = ux*t1[0]+uy*t1[1]+uz*t1[2];
    const double ut2 = ux*t2[0]+uy*t2[1]+uz*t2[2];
    const double Ec  = Ht - c*c/kappa;   // 接触波エネルギー成分 (= K - χ/κ; CPG で K)
    // col0 = r-  (Un-c)
    R[0][0]=1.0;        R[1][0]=ux-c*nx; R[2][0]=uy-c*ny; R[3][0]=uz-c*nz; R[4][0]=Ht-c*Un;
    // col1 = r_s1 (Un)
    R[0][1]=0.0;        R[1][1]=t1[0];   R[2][1]=t1[1];   R[3][1]=t1[2];   R[4][1]=ut1;
    // col2 = r_s2 (Un)
    R[0][2]=0.0;        R[1][2]=t2[0];   R[2][2]=t2[1];   R[3][2]=t2[2];   R[4][2]=ut2;
    // col3 = r_c  (Un)
    R[0][3]=1.0;        R[1][3]=ux;      R[2][3]=uy;      R[3][3]=uz;      R[4][3]=Ec;
    // col4 = r+  (Un+c)
    R[0][4]=1.0;        R[1][4]=ux+c*nx; R[2][4]=uy+c*ny; R[3][4]=uz+c*nz; R[4][4]=Ht+c*Un;
    lambda[0]=Un-c; lambda[1]=Un; lambda[2]=Un; lambda[3]=Un; lambda[4]=Un+c;
}

// A± = R Λ± L (Λ+ = max(λ,0), Λ- = min(λ,0))。Aplus/Aminus に書く。
// 戻り値 = LU 最小ピボット (診断用)。Aminus が nullptr なら A+ のみ。
THERMO_HD double eos_split_jacobian_general(double ux,double uy,double uz,
                                            double nx,double ny,double nz,
                                            double c, double Ht, double kappa,
                                            double Aplus[5][5], double Aminus[5][5])
{
    double R[5][5], L[5][5], lam[5];
    eos_eigvecs_general(ux,uy,uz,nx,ny,nz,c,Ht,kappa,R,lam);
    const double minpiv = eos_inv5(R, L);
    for (int i=0;i<5;i++){
        for (int j=0;j<5;j++){
            double ap=0.0, am=0.0;
            for (int k=0;k<5;k++){
                const double rl = R[i][k]*L[k][j];
                const double lp = lam[k] > 0.0 ? lam[k] : 0.0;
                const double lm = lam[k] < 0.0 ? lam[k] : 0.0;
                ap += rl*lp; am += rl*lm;
            }
            Aplus[i][j] = ap;
            if (Aminus) Aminus[i][j] = am;
        }
    }
    return minpiv;
}
