// 一般EOS flux Jacobian 固有系の単体検証 (Level1/2)。
// plan time_integration-general-eos-jacobian.md §6。__CUDACC__ 未定義で純 host ビルド可。
//   Level1: ||LR-I||, ||RΛL - A_FD||/||A_FD||  (CPG/TP, 複数温度/Mach/法線)
//   Level2: A+ + A- = A,  法線反転  A_{-n}^+ = -A_n^- ,  A_{-n}^- = -A_n^+
// ビルド: g++ -O2 -I solver_density_cuda solver_density_cuda/tools/test_eos_jacobian.cpp -o /tmp/teij
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "../cuda_forge/thermo_d.cuh"
#include "../cuda_forge/eos_jacobian_d.cuh"

static SpeciesThermo mkN2(){
    const double lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    SpeciesThermo s; s.MW=0.0280134; s.sigma_LJ=3.621; s.eps_kB=97.53;
    s.Tlo=200.0; s.Tmid=1000.0; s.Thi=6000.0;
    for(int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; }
    return s;
}
static SpeciesThermo N2 = mkN2();
static const double Y1[1]={1.0};

struct Prim { double ro,ux,uy,uz,P,T,c,h,Ht,kappa,e; };

// mode: 0=CPG (gamma,cp), 1=TP (NASA N2)
static Prim prim_from_Q(const double Q[5], int mode){
    Prim p; p.ro=Q[0]; p.ux=Q[1]/Q[0]; p.uy=Q[2]/Q[0]; p.uz=Q[3]/Q[0];
    double ek=0.5*(p.ux*p.ux+p.uy*p.uy+p.uz*p.uz);
    p.e = Q[4]/Q[0] - ek;
    if (mode==0){
        const double gam=1.4, cp=1039.0; double R=cp*(gam-1.0)/gam, cv=cp-R;
        p.T=p.e/cv; p.P=p.ro*R*p.T; p.h=cp*p.T; p.c=sqrt(gam*R*p.T); p.kappa=gam-1.0;
    } else {
        p.T = thermo_T_from_e(&N2,1,Y1,p.e, 300.0, 50.0, 6000.0);
        ThermoDerivatives D; thermo_derivatives_mix(&N2,1,Y1,p.T,&D);
        p.P=p.ro*D.R*p.T; p.h=D.h; p.c=sqrt(D.a2); p.kappa=D.kappa;
    }
    p.Ht = p.h + ek;
    return p;
}
static void Fn(const double Q[5], const double n[3], int mode, double F[5]){
    Prim p=prim_from_Q(Q,mode); double Un=p.ux*n[0]+p.uy*n[1]+p.uz*n[2];
    F[0]=p.ro*Un;
    F[1]=p.ro*p.ux*Un+p.P*n[0];
    F[2]=p.ro*p.uy*Un+p.P*n[1];
    F[3]=p.ro*p.uz*Un+p.P*n[2];
    F[4]=(Q[4]+p.P)*Un;
}
static void A_fd(const double Q[5], const double n[3], int mode, double A[5][5]){
    for(int j=0;j<5;j++){
        double h=1e-6*std::max(std::fabs(Q[j]),1e-3);
        double Qp[5],Qm[5],Fp[5],Fm[5];
        for(int k=0;k<5;k++){Qp[k]=Q[k];Qm[k]=Q[k];} Qp[j]+=h; Qm[j]-=h;
        Fn(Qp,n,mode,Fp); Fn(Qm,n,mode,Fm);
        for(int i=0;i<5;i++) A[i][j]=(Fp[i]-Fm[i])/(2*h);
    }
}
static double fro(double A[5][5]){ double s=0; for(int i=0;i<5;i++)for(int j=0;j<5;j++)s+=A[i][j]*A[i][j]; return sqrt(s); }
static double frodiff(double A[5][5],double B[5][5]){ double s=0; for(int i=0;i<5;i++)for(int j=0;j<5;j++){double d=A[i][j]-B[i][j];s+=d*d;} return sqrt(s); }

// 全 Jacobian A = R Λ L (split せず) と L=R^-1 を作り、||LR-I|| を返す。
static double build_A_and_LRerr(const Prim&p,const double n[3],double A[5][5]){
    double R[5][5],L[5][5],lam[5];
    eos_eigvecs_general(p.ux,p.uy,p.uz,n[0],n[1],n[2],p.c,p.Ht,p.kappa,R,lam);
    eos_inv5(R,L);
    for(int i=0;i<5;i++)for(int j=0;j<5;j++){double s=0;for(int k=0;k<5;k++)s+=R[i][k]*lam[k]*L[k][j];A[i][j]=s;}
    // ||LR-I||
    double e=0; for(int i=0;i<5;i++)for(int j=0;j<5;j++){double s=0;for(int k=0;k<5;k++)s+=L[i][k]*R[k][j];double d=s-(i==j?1.0:0.0);e+=d*d;}
    return sqrt(e);
}

int main(){
    struct Case{const char*name;int mode;double T;double M;double dir;};
    // dir: 0=軸平行 (1,0,0)、1=一般方向
    Case cs[]={
        {"CPG T250 M1.5 axis", 0,250,1.5,0},
        {"CPG T250 M1.5 genN", 0,250,1.5,1},
        {"TP  T250 M0.0",      1,250,0.0,1},
        {"TP  T250 M0.5",      1,250,0.5,1},
        {"TP  T250 M0.98",     1,250,0.98,1},
        {"TP  T250 M1.5",      1,250,1.5,1},
        {"TP  T250 M3.0",      1,250,3.0,1},
        {"TP  T1200 M1.5",     1,1200,1.5,1},  // hi 係数側 (Tmid=1000 の seam は FD 基準が不連続なので避ける)
        {"TP  T2500 M1.5",     1,2500,1.5,1},  // 高温 NASA 区間 (hi 係数)
    };
    printf("%-22s %12s %14s %10s\n","case","||LR-I||","||RLL-Afd||rel","minpiv");
    bool ok=true;
    for(auto&c:cs){
        // build Q
        double ro=0.4;
        double T=c.T;
        double R,cc,e;
        if(c.mode==0){ const double gam=1.4,cp=1039.0; R=cp*(gam-1)/gam; cc=sqrt(gam*R*T); e=(cp-R)*T; }
        else { ThermoDerivatives D; thermo_derivatives_mix(&N2,1,Y1,T,&D); R=D.R; cc=sqrt(D.a2); e=D.e; }
        double n[3]; if(c.dir==0){n[0]=1;n[1]=0;n[2]=0;} else {n[0]=0.6;n[1]=0.7;n[2]=0.39; double l=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); n[0]/=l;n[1]/=l;n[2]/=l;}
        double sp=c.M*cc; // speed along x
        double Q[5]={ro, ro*sp, ro*0.05*cc, 0.0, ro*(e+0.5*((sp)*(sp)+(0.05*cc)*(0.05*cc)))};
        Prim p=prim_from_Q(Q,c.mode);
        double Aeig[5][5]; double lrerr=build_A_and_LRerr(p,n,Aeig);
        double Afd[5][5]; A_fd(Q,n,c.mode,Afd);
        double rel=frodiff(Aeig,Afd)/fro(Afd);
        // minpiv
        double Ap[5][5],Am[5][5]; double mp=eos_split_jacobian_general(p.ux,p.uy,p.uz,n[0],n[1],n[2],p.c,p.Ht,p.kappa,Ap,Am);
        printf("%-22s %12.2e %14.2e %10.2e %s\n",c.name,lrerr,rel,mp,(rel<1e-5&&lrerr<1e-9)?"":" <-- FAIL");
        if(!(rel<1e-5&&lrerr<1e-9)) ok=false;
        // Level2: A+ + A- = A,  flip
        double Asum[5][5]; for(int i=0;i<5;i++)for(int j=0;j<5;j++)Asum[i][j]=Ap[i][j]+Am[i][j];
        double splitErr=frodiff(Asum,Aeig)/fro(Aeig);
        double nn[3]={-n[0],-n[1],-n[2]}; double Bp[5][5],Bm[5][5];
        eos_split_jacobian_general(p.ux,p.uy,p.uz,nn[0],nn[1],nn[2],p.c,p.Ht,p.kappa,Bp,Bm);
        // A_{-n}^+ = -A_n^-
        double flipErr=0; { double M1[5][5]; for(int i=0;i<5;i++)for(int j=0;j<5;j++)M1[i][j]=Bp[i][j]+Am[i][j]; flipErr=fro(M1)/fmax(fro(Aeig),1e-30); }
        if(splitErr>1e-9||flipErr>1e-9){ printf("    Level2: split=%.1e flip=%.1e <-- FAIL\n",splitErr,flipErr); ok=false; }
    }
    printf("\nLevel1/2 %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
