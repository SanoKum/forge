#include "solveNavierStokes.hpp"

using namespace std;
using namespace Eigen;

void calcRoeAverage (flow_float ga, flow_float roL, flow_float roR, flow_float uL, flow_float uR, 
                     flow_float vL, flow_float vR, 
                     flow_float wL, flow_float wR, 
                     flow_float HL, flow_float HR, flow_float cL, flow_float cR,
                     flow_float& ro_a, flow_float& Ux_a, flow_float& Uy_a, flow_float& Uz_a, flow_float& H_a, flow_float& c_a){

    ro_a = sqrt(roL*roR);
    Ux_a= (sqrt(roL)*uL + sqrt(roR)*uR)/(sqrt(roL)+sqrt(roL));
    Uy_a= (sqrt(roL)*vL + sqrt(roR)*vR)/(sqrt(roL)+sqrt(roL));
    Uz_a= (sqrt(roL)*wL + sqrt(roR)*wR)/(sqrt(roL)+sqrt(roL));
    H_a= (sqrt(roL)*HL + sqrt(roR)*HR)/(sqrt(roL)+sqrt(roL));
    c_a= sqrt((ga-1.0)*(H_a-0.5*(Ux_a*Ux_a+Uy_a*Uy_a+Uz_a*Uz_a)));

};


void calcFluxJacobianMatrix (flow_float ga, flow_float nx, flow_float ny, flow_float nz, 
                             flow_float u , flow_float v , flow_float w,
                             //flow_float p , flow_float ro, flow_float H, flow_float c,
                             flow_float H, flow_float c,
                             Eigen::MatrixXd& A,Eigen::MatrixXd& lambda, 
                             Eigen::MatrixXd& R,Eigen::MatrixXd& L) {

    flow_float ek = 0.5*(u*u+v*v+w*w);
    flow_float U  = u*nx+v*ny+w*nz;
    flow_float chi= (ga-1.0)/c;

    A(0,0) = 0.0              ; A(0,1) = nx                ; A(0,2) = ny                ; A(0,3) = nz                ; A(0,4) = 0.0;
    A(1,0) =(ga-1.0)*ek*nx-u*U; A(1,1) = U-(ga-2.0)*u*nx   ; A(1,2) = u*ny-(ga-1.0)*v*nx; A(1,3) = u*nz-(ga-1.0)*w*nx; A(1,4) = (ga-1.0)*nx;
    A(2,0) =(ga-1.0)*ek*ny-v*U; A(2,1) = v*nx-(ga-1.0)*u*ny; A(2,2) = U-(ga-2.0)*v*ny   ; A(2,3) = v*nz-(ga-1.0)*w*ny; A(2,4) = (ga-1.0)*ny;
    A(3,0) =(ga-1.0)*ek*nz-w*U; A(3,1) = w*nx-(ga-1.0)*u*nz; A(3,2) = w*ny-(ga-1.0)*v*nz; A(3,3) = U-(ga-2.0)*w*nz   ; A(3,4) = (ga-1.0)*nz;
    A(4,0) =((ga-1.0)*ek-H)*U ; A(4,1) = H*nx -(ga-1.0)*u*U; A(4,2) = H*ny -(ga-1.0)*v*U; A(4,3) = H*nz -(ga-1.0)*w*U; A(4,4) = ga*U;

    //lambda(0,0) = abs(U-c); lambda(0,1) = 0.0   ; lambda(0,2) = 0.0      ; lambda(0,3) = 0.0   ; lambda(0,4) = 0.0   ;
    //lambda(1,0) = 0.0     ; lambda(1,1) = abs(U); lambda(1,2) = 0.0      ; lambda(1,3) = 0.0   ; lambda(1,4) = 0.0   ;
    //lambda(2,0) = 0.0     ; lambda(2,1) = 0.0   ; lambda(2,2) = abs(U+c) ; lambda(2,3) = 0.0   ; lambda(2,4) = 0.0   ;
    //lambda(3,0) = 0.0     ; lambda(3,1) = 0.0   ; lambda(3,2) = 0.0      ; lambda(3,3) = abs(U); lambda(3,4) = 0.0   ;
    //lambda(4,0) = 0.0     ; lambda(4,1) = 0.0   ; lambda(4,2) = 0.0      ; lambda(4,3) = 0.0   ; lambda(4,4) = abs(U);

    //R(0,0) = 1.0     ; R(0,1) = 1.0   ; R(0,2) = 1.0      ; R(0,3) = 0.0       ; R(0,4) = 0.0      ;
    //R(1,0) = u-c*nx  ; R(1,1) = u     ; R(1,2) = u+c*nx   ; R(1,3) = ny        ; R(1,4) = -nz      ;
    //R(2,0) = v-c*nx  ; R(2,1) = v     ; R(2,2) = v+c*ny   ; R(2,3) = -nx       ; R(2,4) = 0.0      ;
    //R(3,0) = w-c*nx  ; R(3,1) = w     ; R(3,2) = w+c*nz   ; R(3,3) = 0.0       ; R(3,4) = nx       ;
    //R(4,0) = H-c*U*ek; R(4,1) = ek    ; R(4,2) = H+c*U    ; R(4,3) = u*ny-v*nx ; R(4,4) = w*nx-u*nz;

    //Rinv(0,0) = ((ga-1.0)*ek+c*U)/(2*c*c) ; Rinv(0,1) = ((1.0-ga)*u-c*nx)/(2*c*c) ; Rinv(0,2) = ((1.0-ga)*v-c*ny)/(2*c*c) ; Rinv(0,3) = ((1.0-ga)*w-c*nz)/(2*c*c) ; Rinv(0,4) = (ga-1.0)/(2*c*c);
    //Rinv(1,0) = (c*c-(ga-1.0)*ek)/(2*c*c) ; Rinv(1,1) = (ga-1.0)*u/(2*c*c)        ; Rinv(1,2) = (ga-1.0)*v/(2*c*c)        ; Rinv(1,3) = (ga-1.0)*w/(2*c*c)        ; Rinv(1,4) = (1.0-ga)*u/(c*c);
    //Rinv(2,0) = ((ga-1.0)*ek-c*U)/(2*c*c) ; Rinv(2,1) = ((1.0-ga)*u+c*nx)/(2*c*c) ; Rinv(2,2) = ((1.0-ga)*v+c*ny)/(2*c*c) ; Rinv(2,3) = ((1.0-ga)*w+c*nz)/(2*c*c) ; Rinv(2,4) = (ga-1.0)/(2*c*c);
    //Rinv(3,0) = (v-U*ny)/nx               ; Rinv(3,1) = ny                        ; Rinv(3,2) = (ny*ny-1.0)/nx            ; Rinv(3,3) = ny*nz/nx                  ; Rinv(3,4) = 0.0;
    //Rinv(4,0) = (U*nz-w)/nx               ; Rinv(4,1) =-nz                        ; Rinv(4,2) = -ny*nz/nx                 ; Rinv(4,3) = (1-nz*nz)/nx              ; Rinv(4,4) = 0.0;

    lambda(0,0) = abs(U+c); lambda(0,1) = 0.0   ; lambda(0,2) = 0.0      ; lambda(0,3) = 0.0   ; lambda(0,4) = 0.0     ;
    lambda(1,0) = 0.0     ; lambda(1,1) = abs(U); lambda(1,2) = 0.0      ; lambda(1,3) = 0.0   ; lambda(1,4) = 0.0     ;
    lambda(2,0) = 0.0     ; lambda(2,1) = 0.0   ; lambda(2,2) = abs(U)   ; lambda(2,3) = 0.0   ; lambda(2,4) = 0.0     ;
    lambda(3,0) = 0.0     ; lambda(3,1) = 0.0   ; lambda(3,2) = 0.0      ; lambda(3,3) = abs(U); lambda(3,4) = 0.0     ;
    lambda(4,0) = 0.0     ; lambda(4,1) = 0.0   ; lambda(4,2) = 0.0      ; lambda(4,3) = 0.0   ; lambda(4,4) = abs(U-c);

    R(0,0) = 1.0/(sqrt(2.0)*c)                ; R(0,1) = ny/c     ; R(0,2) = nz/c      ; R(0,3) = nx/c     ; R(0,4) = 1.0/(sqrt(2.0)*c)   ;
    R(1,0) = (u/c+nx)/(sqrt(2.0))             ; R(1,1) = u*ny/c+nz; R(1,2) = u*nz/c-ny ; R(1,3) = u*nx/c   ; R(1,4) = (u/c-nx)/(sqrt(2.0)) ;
    R(2,0) = (v/c+ny)/(sqrt(2.0))             ; R(2,1) = v*ny/c   ; R(2,2) = v*nz/c+nx ; R(2,3) = v*nx/c-nz; R(2,4) = (v/c-ny)/(sqrt(2.0)) ;
    R(3,0) = (w/c+nz)/(sqrt(2.0))             ; R(3,1) = w*ny/c-nx; R(3,2) = w*nz/c    ; R(3,3) = w*nx/c+ny; R(3,4) = (w/c-nz)/(sqrt(2.0)) ;
    R(4,0) = (ek/c+1.0/chi+U)/(sqrt(2.0)); R(4,1) = (ek*ny/c+nz*u-nx*w)   ; R(4,2) = ek*nz/c +nx*v-ny*u ; R(4,3) = ek*nx/c +ny*w -nz*v ; R(4,4) = (ek/c+1.0/chi-U)/sqrt(2.0);

    L(0,0) = (chi*ek-U)/sqrt(2.0)             ; L(0,1) = (-chi*u+nx)/sqrt(2.0)  ; L(0,2) = (-chi*v+ny)/sqrt(2.0) ; L(0,3) =  (-chi*w+nz)/sqrt(2.0) ; L(0,4) = chi/sqrt(2.0)  ;
    L(1,0) = ny*(-chi*ek+c)-nz*u+nx*w         ; L(1,1) = ny*chi*u+nz            ; L(1,2) = ny*chi*v              ; L(1,3) = ny*chi*w-nx            ; L(1,4) = -ny*chi ;
    L(2,0) = nz*(-chi*ek+c)-nx*v+ny*u         ; L(2,1) = nz*chi*u-ny            ; L(2,2) = nz*chi*v+nx           ; L(2,3) = nz*chi*w               ; L(2,4) = -nz*chi ;
    L(3,0) = nx*(-chi*ek+c)-ny*w+nz*v         ; L(3,1) = nx*chi*u               ; L(3,2) = nx*chi*v-nz           ; L(3,3) = nx*chi*w+ny            ; L(3,4) = -nx*chi ;
    L(4,0) = (chi*ek+U)/sqrt(2.0)             ; L(4,1) = (-chi*u-nx)/sqrt(2.0)  ; L(4,2) = (-chi*v-ny)/sqrt(2.0) ; L(4,3) =  (-chi*w-nz)/sqrt(2.0) ; L(4,4) = chi/sqrt(2.0)  ;


};

void timeIntegration(int iloop , solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
{
    // time integration
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        geom_float vol = msh.cells[ic].volume;

        vector<flow_float>& ro = v.c["ro"];
        vector<flow_float>& roUx = v.c["roUx"];
        vector<flow_float>& roUy = v.c["roUy"];
        vector<flow_float>& roUz = v.c["roUz"];
        vector<flow_float>& roe  = v.c["roe"];

        vector<flow_float>& roN = v.c["roN"];
        vector<flow_float>& roUxN = v.c["roUxN"];
        vector<flow_float>& roUyN = v.c["roUyN"];
        vector<flow_float>& roUzN = v.c["roUzN"];
        vector<flow_float>& roeN = v.c["roeN"];

        vector<flow_float>& res_ro   = v.c["res_ro"];
        vector<flow_float>& res_roUx = v.c["res_roUx"];
        vector<flow_float>& res_roUy = v.c["res_roUy"];
        vector<flow_float>& res_roUz = v.c["res_roUz"];
        vector<flow_float>& res_roe  = v.c["res_roe"];


        res_ro[ic]   = 0.0;
        res_roUx[ic] = 0.0; 
        res_roUy[ic] = 0.0; 
        res_roUz[ic] = 0.0; 
        res_roe[ic]  = 0.0; 

        // Explicit
        std::string s_mloop;
        for (int mloop=0; mloop<iloop+1; mloop++){
            s_mloop = std::to_string(mloop);

            vector<flow_float>& res_ro_mRK   = v.c["res_ro"  +s_mloop];
            vector<flow_float>& res_roUx_mRK = v.c["res_roUx"+s_mloop];
            vector<flow_float>& res_roUy_mRK = v.c["res_roUy"+s_mloop];
            vector<flow_float>& res_roUz_mRK = v.c["res_roUz"+s_mloop];
            vector<flow_float>& res_roe_mRK  = v.c["res_roe" +s_mloop];

            res_ro[ic]   += cfg.beta[iloop][mloop]*res_ro_mRK[ic];
            res_roUx[ic] += cfg.beta[iloop][mloop]*res_roUx_mRK[ic];
            res_roUy[ic] += cfg.beta[iloop][mloop]*res_roUy_mRK[ic];
            res_roUz[ic] += cfg.beta[iloop][mloop]*res_roUz_mRK[ic];
            res_roe[ic]  += cfg.beta[iloop][mloop]*res_roe_mRK[ic];
        }
        ro[ic]   = roN[ic]   + res_ro[ic]*cfg.dt/vol;
        roUx[ic] = roUxN[ic] + res_roUx[ic]*cfg.dt/vol;
        roUy[ic] = roUyN[ic] + res_roUy[ic]*cfg.dt/vol;
        roUz[ic] = roUzN[ic] + res_roUz[ic]*cfg.dt/vol;
        roe[ic]  = roeN[ic]  + res_roe[ic]*cfg.dt/vol;
    }
}

void convectiveFlux(int iRK , solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
{
    geom_float ss;
    vector<geom_float> sv(3);
    vector<geom_float> nv(3); // normal

    vector<geom_float> pcent(3);
    vector<geom_float> c0cent(3);
    vector<geom_float> c1cent(3);
    geom_int ic0;
    geom_int ic1;

    flow_float US;
    flow_float roUS;

    vector<flow_float>& Ux = v.c["Ux"];
    vector<flow_float>& Uy = v.c["Uy"];
    vector<flow_float>& Uz = v.c["Uz"];

    vector<flow_float>& dUxdx = v.c["dUxdx"];
    vector<flow_float>& dUydx = v.c["dUydx"];
    vector<flow_float>& dUzdx = v.c["dUzdx"];

    vector<flow_float>& dUxdy = v.c["dUxdy"];
    vector<flow_float>& dUydy = v.c["dUydy"];
    vector<flow_float>& dUzdy = v.c["dUzdy"];

    vector<flow_float>& dUxdz = v.c["dUxdz"];
    vector<flow_float>& dUydz = v.c["dUydz"];
    vector<flow_float>& dUzdz = v.c["dUzdz"];

    vector<flow_float>& dPdx = v.c["dPdx"];
    vector<flow_float>& dPdy = v.c["dPdy"];
    vector<flow_float>& dPdz = v.c["dPdz"];

    vector<flow_float>& UxN = v.c["UxN"];
    vector<flow_float>& UyN = v.c["UyN"];
    vector<flow_float>& UzN = v.c["UzN"];

//    vector<flow_float>& divU_vol = v.c["divU*vol"];
//    vector<flow_float>& divU_star = v.c["divU_star"];
//
    vector<flow_float>& P = v.c["P"];
    vector<flow_float>& T = v.c["T"];

    vector<flow_float>& ro = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];

    vector<flow_float>& sonic = v.c["sonic"];
    vector<flow_float>& Ht    = v.c["Ht"];

    vector<flow_float>& roN = v.c["roN"];
    vector<flow_float>& roUxN = v.c["roUxN"];
    vector<flow_float>& roUyN = v.c["roUyN"];
    vector<flow_float>& roUzN = v.c["roUzN"];
    vector<flow_float>& roeN = v.c["roeN"];

    vector<flow_float>& res_ro   = v.c["res_ro"];
    vector<flow_float>& res_roUx = v.c["res_roUx"];
    vector<flow_float>& res_roUy = v.c["res_roUy"];
    vector<flow_float>& res_roUz = v.c["res_roUz"];
    vector<flow_float>& res_roe  = v.c["res_roe"];

    // residuals for runge kutta 
    std::string s_iRK;

    s_iRK = std::to_string(iRK);

    vector<flow_float>& res_ro_iRK   = v.c["res_ro"  +s_iRK];
    vector<flow_float>& res_roUx_iRK = v.c["res_roUx"+s_iRK];
    vector<flow_float>& res_roUy_iRK = v.c["res_roUy"+s_iRK];
    vector<flow_float>& res_roUz_iRK = v.c["res_roUz"+s_iRK];
    vector<flow_float>& res_roe_iRK  = v.c["res_roe" +s_iRK];

    vector<flow_float>& Uxp = v.p["Ux"];
    vector<flow_float>& Uyp = v.p["Uy"];
    vector<flow_float>& Uzp = v.p["Uz"];

    vector<flow_float>& rop = v.p["ro"];
    vector<flow_float>& USp = v.p["US"];
    vector<flow_float>& Pp  = v.p["P"];
    vector<flow_float>& Tp  = v.p["T"];

    vector<flow_float>& fxp  = v.p["fx"];

    flow_float ro_a;
    flow_float Ux_a;
    flow_float Uy_a;
    flow_float Uz_a;
    flow_float H_a;
    flow_float c_a;

    Eigen::MatrixXd A(5,5);
    Eigen::MatrixXd lambda(5,5);
    Eigen::MatrixXd R(5,5);
    Eigen::MatrixXd L(5,5);
    Eigen::VectorXd delQ(5);
    Eigen::VectorXd difQ(5);

    vector<flow_float> E1;
    vector<flow_float> E2;

    E1.resize(5); E2.resize(5);

    flow_float Un;

    flow_float Ht0 , Ht1;

    for (geom_int ic=0 ; ic<msh.nCells_all; ic++)
    {
        res_ro_iRK[ic] = 0.0;
        res_roUx_iRK[ic] = 0.0;
        res_roUy_iRK[ic] = 0.0;
        res_roUz_iRK[ic] = 0.0;
        res_roe_iRK[ic] = 0.0;
    }

    // normal plane
    for (geom_int ip=0 ; ip<msh.nPlanes ; ip++)
    {

        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        pcent   = msh.planes[ip].centCoords;

        c0cent  = msh.cells[ic0].centCoords;
        c1cent  = msh.cells[ic1].centCoords;

        if (cfg.convMethod == 1){ // Roe scheme

            calcRoeAverage(cfg.gamma, ro[ic0]    , ro[ic1]  , Ux[ic0]   , Ux[ic1], 
                                      Uy[ic0]    , Uy[ic1]  , Uz[ic0]   , Uz[ic1],
                                      Ht[ic0]    , Ht[ic1]  , sonic[ic0], sonic[ic1],
                                      ro_a, Ux_a, Uy_a, Uz_a, H_a, c_a) ;

            calcFluxJacobianMatrix(cfg.gamma, sv[0]/ss, sv[1]/ss, sv[2]/ss, Ux_a, Uy_a, Uz_a, H_a, c_a, 
                                   A, lambda, R, L);


            delQ(0) = ro[ic1]  -ro[ic0];
            delQ(1) = roUx[ic1]-roUx[ic0];
            delQ(2) = roUy[ic1]-roUy[ic0];
            delQ(3) = roUz[ic1]-roUz[ic0];
            delQ(4) = roe[ic1] -roe[ic0];
            difQ = L*delQ;
            difQ = lambda*difQ;
            difQ = R*difQ;

            //temp Un =( (fxp[ip]*Ux[ic0] +(1.0-fxp[ip])*Ux[ic1])*sv[0]
            //temp      +(fxp[ip]*Uy[ic0] +(1.0-fxp[ip])*Uy[ic1])*sv[1]
            //temp      +(fxp[ip]*Uz[ic0] +(1.0-fxp[ip])*Uz[ic1])*sv[2])/ss;
            Un =0.5*( (Ux[ic0] +Ux[ic1])*sv[0]
                     +(Uy[ic0] +Uy[ic1])*sv[1]
                     +(Uz[ic0] +Uz[ic1])*sv[2])/ss;



            res_ro_iRK[ic0]   -= 0.5*(ro[ic0]  +ro[ic1])*Un*ss - difQ(0)*ss;
            res_roUx_iRK[ic0] -= 0.5*(roUx[ic0]+roUx[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[0] - difQ(1)*ss;
            res_roUy_iRK[ic0] -= 0.5*(roUy[ic0]+roUy[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[1] - difQ(2)*ss;
            res_roUz_iRK[ic0] -= 0.5*(roUz[ic0]+roUz[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[2] - difQ(3)*ss;
            res_roe_iRK[ic0]  -= 0.5*(ro[ic0]*Ht[ic0]+ro[ic1]*Ht[ic1])*Un*ss - difQ(4)*ss;


            res_ro_iRK[ic1]   += 0.5*(ro[ic0]  +ro[ic1])*Un*ss - difQ(0)*ss;
            res_roUx_iRK[ic1] += 0.5*(roUx[ic0]+roUx[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[0] - difQ(1)*ss;
            res_roUy_iRK[ic1] += 0.5*(roUy[ic0]+roUy[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[1] - difQ(2)*ss;
            res_roUz_iRK[ic1] += 0.5*(roUz[ic0]+roUz[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[2] - difQ(3)*ss;
            res_roe_iRK[ic1]  += 0.5*(ro[ic0]*Ht[ic0]+ro[ic1]*Ht[ic1])*Un*ss - difQ(4)*ss;

        } else if (cfg.convMethod == 2){ // SD-SLAU
            flow_float u_p = Ux[ic0];
            flow_float v_p = Uy[ic0];
            flow_float w_p = Uz[ic0];
            flow_float h_p = Ht[ic0];
            flow_float u_m = Ux[ic1];
            flow_float v_m = Uy[ic1];
            flow_float w_m = Uz[ic1];
            flow_float h_m = Ht[ic1];

            flow_float U_p = ((Ux[ic0])*sv[0] +(Uy[ic0])*sv[1] +(Uz[ic0])*sv[2])/ss;
            flow_float U_m = ((Ux[ic1])*sv[0] +(Uy[ic1])*sv[1] +(Uz[ic1])*sv[2])/ss;

            flow_float roU_p = ro[ic0]*U_p;
            flow_float roU_m = ro[ic1]*U_m;

            flow_float P_p = P[ic0];
            flow_float ro_p= ro[ic0];
            flow_float P_m = P[ic1];
            flow_float ro_m = ro[ic1];
            flow_float c_bar = sqrt(cfg.gamma*(P_p+P_m)/(ro_p+ro_m));
            flow_float U_bar = 0.5*(U_p + U_m);
            flow_float ro_bar = 0.5*(ro_p + ro_m);
            flow_float M_bar = U_bar/c_bar;
            flow_float M_p = U_p/c_bar;
            flow_float M_m = U_m/c_bar;

            flow_float ro_del= ro_p - ro_m;
            flow_float U_del= U_p - U_m;
            flow_float P_del= P_p - P_m;

            flow_float beta_p, beta_m;

            if (abs(M_m)<=1.0){
                beta_m = 0.25*pow((M_m-1.0)*(2.0+M_m),2.0);
            } else {
                if (-M_m>0) {
                    beta_m = 1.0;
                } else {
                    beta_m = 0.0;
                };
            }

            if (abs(M_p)<=1.0){
                beta_p = 0.25*pow((M_p+1.0)*(2.0-M_p),2.0);
            } else {
                if (M_p>0) {
                    beta_p = 1.0;
                } else {
                    beta_p = 0.0;
                };
            }

            flow_float zero = 0.0;
            flow_float one  = 1.0;
            flow_float half = 0.5;

            flow_float g = -std::max(min(M_p,zero),-one)*min(max(M_m,zero),one);
            flow_float U_bar_abs = (ro_m*abs(U_m) + ro_p*abs(U_p))/(ro_m+ro_p);
            flow_float U_bar_p_abs = (one-g)*abs(U_bar)+g*abs(U_p) ;
            flow_float U_bar_m_abs = (one-g)*abs(U_bar)+g*abs(U_m) ;
            flow_float M_hat = std::min(one, sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))/c_bar);
            flow_float chi = (1.0-M_hat)*(1.0-M_hat);

            flow_float pnew = half*(P_p+P_m) +half*(beta_p-beta_m)*(P_p-P_m)
                              +(one-chi)*(beta_p+beta_m-one)*half*(P_m+P_p); 

            flow_float m = 0.5*(ro_p*(U_p+U_bar_p_abs)+ro_m*(U_m-U_bar_m_abs)-chi/c_bar*P_del);

            res_ro_iRK[ic0]   -= 0.5*(m+abs(m))     +0.5*(m-abs(m));
            res_roUx_iRK[ic0] -= 0.5*(m+abs(m))*u_p +0.5*(m-abs(m))*u_m +pnew*sv[0]/ss;
            res_roUy_iRK[ic0] -= 0.5*(m+abs(m))*v_p +0.5*(m-abs(m))*v_m +pnew*sv[1]/ss;
            res_roUz_iRK[ic0] -= 0.5*(m+abs(m))*w_p +0.5*(m-abs(m))*w_m +pnew*sv[2]/ss;
            res_roe_iRK[ic0]  -= 0.5*(m+abs(m))*h_p +0.5*(m-abs(m))*h_m;

            res_ro_iRK[ic1]   += 0.5*(m+abs(m))     +0.5*(m-abs(m));
            res_roUx_iRK[ic1] += 0.5*(m+abs(m))*u_p +0.5*(m-abs(m))*u_m +pnew*sv[0]/ss;
            res_roUy_iRK[ic1] += 0.5*(m+abs(m))*v_p +0.5*(m-abs(m))*v_m +pnew*sv[1]/ss;
            res_roUz_iRK[ic1] += 0.5*(m+abs(m))*w_p +0.5*(m-abs(m))*w_m +pnew*sv[2]/ss;
            res_roe_iRK[ic1]  += 0.5*(m+abs(m))*h_p +0.5*(m-abs(m))*h_m;

        } else if (cfg.convMethod == 3){ // KEEP scheme
            flow_float Ctilde  = 0.5*(ro[ic0]+ro[ic1])*0.5*( (Ux[ic0]+Ux[ic1])*sv[0]
                                                            +(Uy[ic0]+Uy[ic1])*sv[1]
                                                            +(Uz[ic0]+Uz[ic1])*sv[2] );
            flow_float Mtildex = Ctilde*(Ux[ic0]+Ux[ic1])*0.5;
            flow_float Mtildey = Ctilde*(Uy[ic0]+Uy[ic1])*0.5;
            flow_float Mtildez = Ctilde*(Uz[ic0]+Uz[ic1])*0.5;

            flow_float Ktilde = Ctilde*0.5*(Ux[ic0]*Ux[ic1] +Uy[ic0]*Uy[ic1] +Uz[ic0]*Uz[ic1]);
            flow_float Itilde = Ctilde*0.5*(P[ic0]/ro[ic0] +P[ic1]/ro[ic1])/(cfg.gamma-1.0);

            flow_float Gtildex = 0.5*(P[ic0]+P[ic1])*sv[0];
            flow_float Gtildey = 0.5*(P[ic0]+P[ic1])*sv[1];
            flow_float Gtildez = 0.5*(P[ic0]+P[ic1])*sv[2];

            flow_float Ptilde = 0.5*((Ux[ic0]*P[ic1] + Ux[ic1]*P[ic0])*sv[0]
                                    +(Uy[ic0]*P[ic1] + Uy[ic1]*P[ic0])*sv[1]
                                    +(Uz[ic0]*P[ic1] + Uz[ic1]*P[ic0])*sv[2]);

            res_ro_iRK[ic0]   -= Ctilde;
            res_roUx_iRK[ic0] -= Mtildex + Gtildex;
            res_roUy_iRK[ic0] -= Mtildey + Gtildey;
            res_roUz_iRK[ic0] -= Mtildez + Gtildez;
            res_roe_iRK[ic0]  -= Ktilde + Itilde + Ptilde;

            res_ro_iRK[ic1]   += Ctilde;
            res_roUx_iRK[ic1] += Mtildex + Gtildex;
            res_roUy_iRK[ic1] += Mtildey + Gtildey;
            res_roUz_iRK[ic1] += Mtildez + Gtildez;
            res_roe_iRK[ic1]  += Ktilde + Itilde + Ptilde;
        }
    }
}


void solveNavierStokes(solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
{
    geom_float ss;
    vector<geom_float> sv(3);
    vector<geom_float> nv(3); // normal

    vector<geom_float> sv_dia(3); // diagonal
    vector<geom_float> sv_nodia(3); // non-diagonal

    vector<geom_float> pcent(3);
    vector<geom_float> c0cent(3);
    vector<geom_float> c1cent(3);
    geom_int ic0;
    geom_int ic1;
    geom_float dn;
    flow_float temp;


    geom_int ip_loc0;
    geom_int ip_loc1;
    
    vector<geom_float> dccv(3);
    geom_float dcc;

    vector<geom_float> dc0pv(3);
    geom_float dc0p;

    vector<geom_float> dc1pv(3);
    geom_float dc1p;

    geom_float dccv_dot_sv;

    vector<geom_float> deltav(3);
    geom_float delta;

    //geom_float temp_ndia;
    //geom_float temp_ndia_x;
    //geom_float temp_ndia_y;
    //geom_float temp_ndia_z;

    flow_float US;
    flow_float roUS;
    flow_float rof;
    flow_float volf;
    flow_float mass;

    flow_float sonicf;
    flow_float mach;

    geom_float cosT;

    vector<flow_float>& Ux = v.c["Ux"];
    vector<flow_float>& Uy = v.c["Uy"];
    vector<flow_float>& Uz = v.c["Uz"];

    vector<flow_float>& dUxdx = v.c["dUxdx"];
    vector<flow_float>& dUydx = v.c["dUydx"];
    vector<flow_float>& dUzdx = v.c["dUzdx"];

    vector<flow_float>& dUxdy = v.c["dUxdy"];
    vector<flow_float>& dUydy = v.c["dUydy"];
    vector<flow_float>& dUzdy = v.c["dUzdy"];

    vector<flow_float>& dUxdz = v.c["dUxdz"];
    vector<flow_float>& dUydz = v.c["dUydz"];
    vector<flow_float>& dUzdz = v.c["dUzdz"];

    vector<flow_float>& dPdx = v.c["dPdx"];
    vector<flow_float>& dPdy = v.c["dPdy"];
    vector<flow_float>& dPdz = v.c["dPdz"];

    vector<flow_float>& UxN = v.c["UxN"];
    vector<flow_float>& UyN = v.c["UyN"];
    vector<flow_float>& UzN = v.c["UzN"];

//    vector<flow_float>& divU_vol = v.c["divU*vol"];
//    vector<flow_float>& divU_star = v.c["divU_star"];
//
    vector<flow_float>& P = v.c["P"];
    vector<flow_float>& T = v.c["T"];

    vector<flow_float>& ro = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];

    vector<flow_float>& sonic = v.c["sonic"];
    vector<flow_float>& Ht    = v.c["Ht"];

    vector<flow_float>& roN = v.c["roN"];
    vector<flow_float>& roUxN = v.c["roUxN"];
    vector<flow_float>& roUyN = v.c["roUyN"];
    vector<flow_float>& roUzN = v.c["roUzN"];
    vector<flow_float>& roeN = v.c["roeN"];

    vector<flow_float>& res_ro   = v.c["res_ro"];
    vector<flow_float>& res_roUx = v.c["res_roUx"];
    vector<flow_float>& res_roUy = v.c["res_roUy"];
    vector<flow_float>& res_roUz = v.c["res_roUz"];
    vector<flow_float>& res_roe  = v.c["res_roe"];

    vector<flow_float>& Uxp = v.p["Ux"];
    vector<flow_float>& Uyp = v.p["Uy"];
    vector<flow_float>& Uzp = v.p["Uz"];

    vector<flow_float>& rop = v.p["ro"];
    vector<flow_float>& USp = v.p["US"];
    vector<flow_float>& Pp  = v.p["P"];
    vector<flow_float>& Tp  = v.p["T"];

    vector<flow_float>& fxp  = v.p["fx"];

    flow_float ro_a;
    flow_float Ux_a;
    flow_float Uy_a;
    flow_float Uz_a;
    flow_float H_a;
    flow_float c_a;

    Eigen::MatrixXd A(5,5);
    Eigen::MatrixXd lambda(5,5);
    Eigen::MatrixXd R(5,5);
    Eigen::MatrixXd L(5,5);
    Eigen::VectorXd delQ(5);
    Eigen::VectorXd difQ(5);

    vector<flow_float> E1;
    vector<flow_float> E2;

    E1.resize(5); E2.resize(5);

    flow_float Un;

    flow_float Ht0 , Ht1;

    for (geom_int ic=0 ; ic<msh.nCells_all; ic++)
    {
        res_ro[ic] = 0.0;
        res_roUx[ic] = 0.0;
        res_roUy[ic] = 0.0;
        res_roUz[ic] = 0.0;
        res_roe[ic] = 0.0;
        
    }

    // normal plane
//ghost>
    //for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    for (geom_int ip=0 ; ip<msh.nPlanes ; ip++)
//ghost<
    {

        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        pcent   = msh.planes[ip].centCoords;

        c0cent  = msh.cells[ic0].centCoords;
        c1cent  = msh.cells[ic1].centCoords;

        if (cfg.convMethod == 1){ // Roe scheme
//            if (ic1==6042){
//                cout << "ic0 U " << ic0 << " "<< Ux[ic0] << " " << Uy[ic0] << " " << Uz[ic0] << endl;
//                cout << "ic1 U " << ic1 << " "<< Ux[ic1] << " " << Uy[ic1] << " " << Uz[ic1] << endl;
//                cout << "    Ht,c,ro "  << " "<< Ht[ic1] << " " << sonic[ic1] << " " << ro[ic1] << endl;
//                cout << "    fxp "  << " "<< fxp[ip] << endl;
//            }

            calcRoeAverage(cfg.gamma, ro[ic0]    , ro[ic1]  , Ux[ic0]   , Ux[ic1], 
                                      Uy[ic0]    , Uy[ic1]  , Uz[ic0]   , Uz[ic1],
                                      Ht[ic0]    , Ht[ic1]  , sonic[ic0], sonic[ic1],
                                      ro_a, Ux_a, Uy_a, Uz_a, H_a, c_a) ;

            calcFluxJacobianMatrix(cfg.gamma, sv[0]/ss, sv[1]/ss, sv[2]/ss, Ux_a, Uy_a, Uz_a, H_a, c_a, 
                                   A, lambda, R, L);


            delQ(0) = ro[ic1]  -ro[ic0];
            delQ(1) = roUx[ic1]-roUx[ic0];
            delQ(2) = roUy[ic1]-roUy[ic0];
            delQ(3) = roUz[ic1]-roUz[ic0];
            delQ(4) = roe[ic1] -roe[ic0];
            difQ = L*delQ;
            difQ = lambda*difQ;
            difQ = R*difQ;

            //temp Un =( (fxp[ip]*Ux[ic0] +(1.0-fxp[ip])*Ux[ic1])*sv[0]
            //temp      +(fxp[ip]*Uy[ic0] +(1.0-fxp[ip])*Uy[ic1])*sv[1]
            //temp      +(fxp[ip]*Uz[ic0] +(1.0-fxp[ip])*Uz[ic1])*sv[2])/ss;
            Un =0.5*( (Ux[ic0] +Ux[ic1])*sv[0]
                     +(Uy[ic0] +Uy[ic1])*sv[1]
                     +(Uz[ic0] +Uz[ic1])*sv[2])/ss;



            res_ro[ic0]   -= 0.5*(ro[ic0]  +ro[ic1])*Un*ss - difQ(0)*ss;
            res_roUx[ic0] -= 0.5*(roUx[ic0]+roUx[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[0] - difQ(1)*ss;
            res_roUy[ic0] -= 0.5*(roUy[ic0]+roUy[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[1] - difQ(2)*ss;
            res_roUz[ic0] -= 0.5*(roUz[ic0]+roUz[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[2] - difQ(3)*ss;
            res_roe[ic0]  -= 0.5*(ro[ic0]*Ht[ic0]+ro[ic1]*Ht[ic1])*Un*ss - difQ(4)*ss;


            res_ro[ic1]   += 0.5*(ro[ic0]  +ro[ic1])*Un*ss - difQ(0)*ss;
            res_roUx[ic1] += 0.5*(roUx[ic0]+roUx[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[0] - difQ(1)*ss;
            res_roUy[ic1] += 0.5*(roUy[ic0]+roUy[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[1] - difQ(2)*ss;
            res_roUz[ic1] += 0.5*(roUz[ic0]+roUz[ic1])*Un*ss +0.5*(P[ic0]+P[ic1])*sv[2] - difQ(3)*ss;
            res_roe[ic1]  += 0.5*(ro[ic0]*Ht[ic0]+ro[ic1]*Ht[ic1])*Un*ss - difQ(4)*ss;

        } else if (cfg.convMethod == 2){ // SD-SLAU
            flow_float u_p = Ux[ic0];
            flow_float v_p = Uy[ic0];
            flow_float w_p = Uz[ic0];
            flow_float h_p = Ht[ic0];
            flow_float u_m = Ux[ic1];
            flow_float v_m = Uy[ic1];
            flow_float w_m = Uz[ic1];
            flow_float h_m = Ht[ic1];

            flow_float U_p = ((Ux[ic0])*sv[0] +(Uy[ic0])*sv[1] +(Uz[ic0])*sv[2])/ss;
            flow_float U_m = ((Ux[ic1])*sv[0] +(Uy[ic1])*sv[1] +(Uz[ic1])*sv[2])/ss;

            flow_float roU_p = ro[ic0]*U_p;
            flow_float roU_m = ro[ic1]*U_m;

            flow_float P_p = P[ic0];
            flow_float ro_p= ro[ic0];
            flow_float P_m = P[ic1];
            flow_float ro_m = ro[ic1];
            flow_float c_bar = sqrt(cfg.gamma*(P_p+P_m)/(ro_p+ro_m));
            flow_float U_bar = 0.5*(U_p + U_m);
            flow_float ro_bar = 0.5*(ro_p + ro_m);
            flow_float M_bar = U_bar/c_bar;
            flow_float M_p = U_p/c_bar;
            flow_float M_m = U_m/c_bar;

            flow_float ro_del= ro_p - ro_m;
            flow_float U_del= U_p - U_m;
            flow_float P_del= P_p - P_m;

            flow_float beta_p, beta_m;

            if (abs(M_m)<=1.0){
                beta_m = 0.25*pow((M_m-1.0)*(2.0+M_m),2.0);
            } else {
                if (-M_m>0) {
                    beta_m = 1.0;
                } else {
                    beta_m = 0.0;
                };
            }

            if (abs(M_p)<=1.0){
                beta_p = 0.25*pow((M_p+1.0)*(2.0-M_p),2.0);
            } else {
                if (M_p>0) {
                    beta_p = 1.0;
                } else {
                    beta_p = 0.0;
                };
            }

            flow_float zero = 0.0;
            flow_float one  = 1.0;
            flow_float half = 0.5;

            flow_float g = -std::max(min(M_p,zero),-one)*min(max(M_m,zero),one);
            flow_float U_bar_abs = (ro_m*abs(U_m) + ro_p*abs(U_p))/(ro_m+ro_p);
            flow_float U_bar_p_abs = (one-g)*abs(U_bar)+g*abs(U_p) ;
            flow_float U_bar_m_abs = (one-g)*abs(U_bar)+g*abs(U_m) ;
            flow_float M_hat = std::min(one, sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))/c_bar);
            flow_float chi = (1.0-M_hat)*(1.0-M_hat);

            flow_float pnew = half*(P_p+P_m) +half*(beta_p-beta_m)*(P_p-P_m)
                              +(one-chi)*(beta_p+beta_m-one)*half*(P_m+P_p); 

            flow_float m = 0.5*(ro_p*(U_p+U_bar_p_abs)+ro_m*(U_m-U_bar_m_abs)-chi/c_bar*P_del);

            res_ro[ic0]   -= 0.5*(m+abs(m))     +0.5*(m-abs(m));
            res_roUx[ic0] -= 0.5*(m+abs(m))*u_p +0.5*(m-abs(m))*u_m +pnew*sv[0]/ss;
            res_roUy[ic0] -= 0.5*(m+abs(m))*v_p +0.5*(m-abs(m))*v_m +pnew*sv[1]/ss;
            res_roUz[ic0] -= 0.5*(m+abs(m))*w_p +0.5*(m-abs(m))*w_m +pnew*sv[2]/ss;
            res_roe[ic0]  -= 0.5*(m+abs(m))*h_p +0.5*(m-abs(m))*h_m;

            res_ro[ic1]   += 0.5*(m+abs(m))     +0.5*(m-abs(m));
            res_roUx[ic1] += 0.5*(m+abs(m))*u_p +0.5*(m-abs(m))*u_m +pnew*sv[0]/ss;
            res_roUy[ic1] += 0.5*(m+abs(m))*v_p +0.5*(m-abs(m))*v_m +pnew*sv[1]/ss;
            res_roUz[ic1] += 0.5*(m+abs(m))*w_p +0.5*(m-abs(m))*w_m +pnew*sv[2]/ss;
            res_roe[ic1]  += 0.5*(m+abs(m))*h_p +0.5*(m-abs(m))*h_m;

        } else if (cfg.convMethod == 3){ // KEEP scheme
            flow_float Ctilde  = 0.5*(ro[ic0]+ro[ic1])*0.5*( (Ux[ic0]+Ux[ic1])*sv[0]
                                                            +(Uy[ic0]+Uy[ic1])*sv[1]
                                                            +(Uz[ic0]+Uz[ic1])*sv[2] );
            flow_float Mtildex = Ctilde*(Ux[ic0]+Ux[ic1])*0.5;
            flow_float Mtildey = Ctilde*(Uy[ic0]+Uy[ic1])*0.5;
            flow_float Mtildez = Ctilde*(Uz[ic0]+Uz[ic1])*0.5;

            flow_float Ktilde = Ctilde*0.5*(Ux[ic0]*Ux[ic1] +Uy[ic0]*Uy[ic1] +Uz[ic0]*Uz[ic1]);
            flow_float Itilde = Ctilde*0.5*(P[ic0]/ro[ic0] +P[ic1]/ro[ic1])/(cfg.gamma-1.0);

            flow_float Gtildex = 0.5*(P[ic0]+P[ic1])*sv[0];
            flow_float Gtildey = 0.5*(P[ic0]+P[ic1])*sv[1];
            flow_float Gtildez = 0.5*(P[ic0]+P[ic1])*sv[2];

            flow_float Ptilde = 0.5*((Ux[ic0]*P[ic1] + Ux[ic1]*P[ic0])*sv[0]
                                    +(Uy[ic0]*P[ic1] + Uy[ic1]*P[ic0])*sv[1]
                                    +(Uz[ic0]*P[ic1] + Uz[ic1]*P[ic0])*sv[2]);

            res_ro[ic0]   -= Ctilde;
            res_roUx[ic0] -= Mtildex + Gtildex;
            res_roUy[ic0] -= Mtildey + Gtildey;
            res_roUz[ic0] -= Mtildez + Gtildez;
            res_roe[ic0]  -= Ktilde + Itilde + Ptilde;

            res_ro[ic1]   += Ctilde;
            res_roUx[ic1] += Mtildex + Gtildex;
            res_roUy[ic1] += Mtildey + Gtildey;
            res_roUz[ic1] += Mtildez + Gtildez;
            res_roe[ic1]  += Ktilde + Itilde + Ptilde;
        }
    }

//ghost>    // boundary conditions
//ghost>    for (auto& bc : msh.bconds)
//ghost>    {
//ghost>        vector<flow_float>& rob   = bc.bvar["ro"];
//ghost>        vector<flow_float>& Uxb   = bc.bvar["Ux"];
//ghost>        vector<flow_float>& Uyb   = bc.bvar["Uy"];
//ghost>        vector<flow_float>& Uzb   = bc.bvar["Uz"];
//ghost>        vector<flow_float>& roUxb = bc.bvar["roUx"];
//ghost>        vector<flow_float>& roUyb = bc.bvar["roUy"];
//ghost>        vector<flow_float>& roUzb = bc.bvar["roUz"];
//ghost>        vector<flow_float>& roeb  = bc.bvar["roe"];
//ghost>        vector<flow_float>& Psb   = bc.bvar["Ps"];
//ghost>
//ghost>        geom_int ip0;
//ghost>
//ghost>        //for (geom_int& ip : bc.iPlanes)
//ghost>        for (geom_int ib=0 ; ib<bc.iBPlanes.size() ; ib++)
//ghost>        {
//ghost>            ip0   = bc.iPlanes[ib];
//ghost>            ic0   = bc.iCells[ib];
//ghost>            sv    = msh.planes[ip0].surfVect;
//ghost>            ss    = msh.planes[ip0].surfArea;
//ghost>
//ghost>            US = (Uxb[ib]*sv[0] + Uyb[ib]*sv[1] + Uzb[ib]*sv[2]);
//ghost>
//ghost>            if (cfg.convMethod == 1){ // roe scheme
//ghost>
//ghost>                // keep
//ghost>                //flow_float Ctilde  = 0.5*(ro[ic0]+rob[ib])*0.5*( (Ux[ic0]+Uxb[ib])*sv[0]
//ghost>                //                                                +(Uy[ic0]+Uyb[ib])*sv[1]
//ghost>                //                                                +(Uz[ic0]+Uzb[ib])*sv[2] );
//ghost>                //flow_float Mtildex = Ctilde*(Ux[ic0]+Uxb[ib])*0.5;
//ghost>                //flow_float Mtildey = Ctilde*(Uy[ic0]+Uyb[ib])*0.5;
//ghost>                //flow_float Mtildez = Ctilde*(Uz[ic0]+Uzb[ib])*0.5;
//ghost>
//ghost>                //flow_float Ktilde = Ctilde*0.5*(Ux[ic0]*Uxb[ib] +Uy[ic0]*Uyb[ib] +Uz[ic0]*Uzb[ib]);
//ghost>                //flow_float Itilde = Ctilde*0.5*(P[ic0]/ro[ic0] +Psb[ib]/rob[ib])/(cfg.gamma-1.0);
//ghost>
//ghost>                //flow_float Gtildex = 0.5*(P[ic0]+Psb[ib])*sv[0];
//ghost>                //flow_float Gtildey = 0.5*(P[ic0]+Psb[ib])*sv[1];
//ghost>                //flow_float Gtildez = 0.5*(P[ic0]+Psb[ib])*sv[2];
//ghost>
//ghost>                //flow_float Ptilde = 0.5*((Ux[ic0]*Psb[ib] + Uxb[ib]*P[ic0])*sv[0]
//ghost>                //                        +(Uy[ic0]*Psb[ib] + Uyb[ib]*P[ic0])*sv[1]
//ghost>                //                        +(Uz[ic0]*Psb[ib] + Uzb[ib]*P[ic0])*sv[2]);
//ghost>
//ghost>                //res_ro[ic0]   -= Ctilde;
//ghost>                //res_roUx[ic0] -= Mtildex + Gtildex;
//ghost>                //res_roUy[ic0] -= Mtildey + Gtildey;
//ghost>                //res_roUz[ic0] -= Mtildez + Gtildez;
//ghost>                //res_roe[ic0]  -= Ktilde + Itilde + Ptilde;
//ghost>
//ghost>                flow_float Ctilde  = rob[ib]*( (Uxb[ib])*sv[0]
//ghost>                                              +(Uyb[ib])*sv[1]
//ghost>                                              +(Uzb[ib])*sv[2] );
//ghost>                flow_float Mtildex = Ctilde*(Uxb[ib]);
//ghost>                flow_float Mtildey = Ctilde*(Uyb[ib]);
//ghost>                flow_float Mtildez = Ctilde*(Uzb[ib]);
//ghost>
//ghost>                flow_float Ktilde = Ctilde*0.5*(Uxb[ib]*Uxb[ib] +Uyb[ib]*Uyb[ib] +Uzb[ib]*Uzb[ib]);
//ghost>                flow_float Itilde = Ctilde*(Psb[ib]/rob[ib])/(cfg.gamma-1.0);
//ghost>
//ghost>                flow_float Gtildex = (Psb[ib])*sv[0];
//ghost>                flow_float Gtildey = (Psb[ib])*sv[1];
//ghost>                flow_float Gtildez = (Psb[ib])*sv[2];
//ghost>
//ghost>                flow_float Ptilde = ((Uxb[ib]*Psb[ib])*sv[0]
//ghost>                                    +(Uyb[ib]*Psb[ib])*sv[1]
//ghost>                                    +(Uzb[ib]*Psb[ib])*sv[2]);
//ghost>
//ghost>                res_ro[ic0]   -= Ctilde;
//ghost>                res_roUx[ic0] -= Mtildex + Gtildex;
//ghost>                res_roUy[ic0] -= Mtildey + Gtildey;
//ghost>                res_roUz[ic0] -= Mtildez + Gtildez;
//ghost>                res_roe[ic0]  -= Ktilde + Itilde + Ptilde;
//ghost>
//ghost>
//ghost>//                flow_float Ht_b;
//ghost>//                flow_float roe_b;
//ghost>//                flow_float sonic_b;
//ghost>//
//ghost>//                Ht_b = (roeb[ib] + Psb[ib])/rob[ib];
//ghost>//                
//ghost>//                sonic_b = sqrt(cfg.gamma*Psb[ib]/rob[ib]);
//ghost>//
//ghost>//                calcRoeAverage(cfg.gamma, ro[ic0]    , rob[ib], Ux[ic0]   , Uxb[ib], 
//ghost>//                                          Uy[ic0]    , Uyb[ib], Uz[ic0]   , Uzb[ib],
//ghost>//                                          Ht[ic0] , Ht_b   , sonic[ic0], sonic_b,
//ghost>//                                          ro_a, Ux_a , Uy_a   , Uz_a, H_a , c_a) ;
//ghost>//
//ghost>//                calcFluxJacobianMatrix(cfg.gamma, sv[0]/ss, sv[1]/ss, sv[2]/ss, Ux_a, Uy_a, Uz_a, H_a, c_a, 
//ghost>//                                       A, lambda, R, L);
//ghost>//
//ghost>//
//ghost>//                delQ(0) = rob[ib]         -ro[ic0]  ;
//ghost>//                delQ(1) = rob[ib]*Uxb[ib] -roUx[ic0];
//ghost>//                delQ(2) = rob[ib]*Uyb[ib] -roUy[ic0];
//ghost>//                delQ(3) = rob[ib]*Uzb[ib] -roUz[ic0];
//ghost>//                delQ(4) = roeb[ib]        -roe[ic0];
//ghost>//                          
//ghost>//
//ghost>//                difQ = L*delQ;
//ghost>//                difQ = lambda*difQ;
//ghost>//                difQ = R*difQ;
//ghost>//
//ghost>//                Un =0.5*( (Ux[ic0] +Uxb[ib])*sv[0]
//ghost>//                         +(Uy[ic0] +Uyb[ib])*sv[1]
//ghost>//                         +(Uz[ic0] +Uzb[ib])*sv[2])/ss;
//ghost>//
//ghost>//                res_ro[ic0]   -= 0.5*(ro[ic0]  +rob[ib])*Un*ss - difQ(0)*ss;
//ghost>//                res_roUx[ic0] -= 0.5*(roUx[ic0]+rob[ib]*Uxb[ib])*Un*ss +0.5*(P[ic0]+Psb[ib])*sv[0] - difQ(1)*ss;
//ghost>//                res_roUy[ic0] -= 0.5*(roUy[ic0]+rob[ib]*Uyb[ib])*Un*ss +0.5*(P[ic0]+Psb[ib])*sv[1] - difQ(2)*ss;
//ghost>//                res_roUz[ic0] -= 0.5*(roUz[ic0]+rob[ib]*Uzb[ib])*Un*ss +0.5*(P[ic0]+Psb[ib])*sv[2] - difQ(3)*ss;
//ghost>//                res_roe[ic0]  -= 0.5*(ro[ic0]*Ht[ic0]+rob[ib]*Ht_b)*Un*ss - difQ(4)*ss;
//ghost>//
//ghost>            } else if (cfg.convMethod == 2){ // SLAU
//ghost>                flow_float u_p = Ux[ic0];
//ghost>                flow_float v_p = Uy[ic0];
//ghost>                flow_float w_p = Uz[ic0];
//ghost>                flow_float h_p = Ht[ic0];
//ghost>                flow_float u_m = Uxb[ib];
//ghost>                flow_float v_m = Uyb[ib];
//ghost>                flow_float w_m = Uzb[ib];
//ghost>                flow_float h_m = (roeb[ib] + Psb[ib])/rob[ib];
//ghost>
//ghost>                flow_float U_p = ((Ux[ic0])*sv[0] +(Uy[ic0])*sv[1] +(Uz[ic0])*sv[2])/ss;
//ghost>                flow_float U_m = ((Uxb[ib])*sv[0] +(Uyb[ib])*sv[1] +(Uzb[ib])*sv[2])/ss;
//ghost>
//ghost>                flow_float roU_p = ro[ic0]*U_p;
//ghost>                flow_float roU_m = rob[ib]*U_m;
//ghost>
//ghost>                flow_float P_p = P[ic0];
//ghost>                flow_float P_m = Psb[ib];
//ghost>                flow_float ro_p = ro[ic0];
//ghost>                flow_float ro_m = rob[ib];
//ghost>                flow_float c_bar = sqrt(cfg.gamma*(P_p+P_m)/(ro_p+ro_m));
//ghost>                flow_float U_bar = 0.5*(U_p + U_m);
//ghost>                flow_float ro_bar = 0.5*(ro_p + ro_m);
//ghost>                flow_float M_bar = U_bar/c_bar;
//ghost>                flow_float M_p = U_p/c_bar;
//ghost>                flow_float M_m = U_m/c_bar;
//ghost>
//ghost>                flow_float P_del= P_p - P_m;
//ghost>
//ghost>                flow_float beta_p, beta_m;
//ghost>
//ghost>                if (abs(M_m)<=1.0){
//ghost>                    beta_m = 0.25*pow((M_m-1.0)*(2.0+M_m),2.0);
//ghost>                } else {
//ghost>                    if (-M_m>0) {
//ghost>                        beta_m = 1.0;
//ghost>                    } else {
//ghost>                        beta_m = 0.0;
//ghost>                    };
//ghost>                }
//ghost>
//ghost>                if (abs(M_p)<=1.0){
//ghost>                    beta_p = 0.25*pow((M_p+1.0)*(2.0-M_p),2.0);
//ghost>                } else {
//ghost>                    if (M_p>0) {
//ghost>                        beta_p = 1.0;
//ghost>                    } else {
//ghost>                        beta_p = 0.0;
//ghost>                    };
//ghost>                }
//ghost>
//ghost>                flow_float zero = 0.0;
//ghost>                flow_float one  = 1.0;
//ghost>                flow_float half = 0.5;
//ghost>
//ghost>                flow_float g = -std::max(min(M_p,zero),-one)*min(max(M_m,zero),one);
//ghost>                flow_float U_bar_abs = (ro_m*abs(U_m) + ro_p*abs(U_p))/(ro_m+ro_p);
//ghost>                flow_float U_bar_p_abs = (one-g)*abs(U_bar)+g*abs(U_p) ;
//ghost>                flow_float U_bar_m_abs = (one-g)*abs(U_bar)+g*abs(U_m) ;
//ghost>                flow_float M_hat = std::min(one, sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))/c_bar);
//ghost>                flow_float chi = (1.0-M_hat)*(1.0-M_hat);
//ghost>
//ghost>                flow_float pnew = half*(P_m+P_p) +half*(beta_p-beta_m)*(P_p-P_m)
//ghost>                                  +(one-chi)*(beta_p+beta_m-one)*half*(P_m+P_p); 
//ghost>
//ghost>                flow_float m = 0.5*(ro_p*(U_p+U_bar_p_abs)+ro_m*(U_m-U_bar_m_abs)-chi/c_bar*P_del);
//ghost>
//ghost>                res_ro[ic0]   -= 0.5*(m+abs(m))     +0.5*(m-abs(m));
//ghost>                res_roUx[ic0] -= 0.5*(m+abs(m))*u_p +0.5*(m-abs(m))*u_m +pnew*sv[0]/ss;
//ghost>                res_roUy[ic0] -= 0.5*(m+abs(m))*v_p +0.5*(m-abs(m))*v_m +pnew*sv[1]/ss;
//ghost>                res_roUz[ic0] -= 0.5*(m+abs(m))*w_p +0.5*(m-abs(m))*w_m +pnew*sv[2]/ss;
//ghost>                res_roe[ic0]  -= 0.5*(m+abs(m))*h_p +0.5*(m-abs(m))*h_m;
//ghost>
//ghost>            } else if (cfg.convMethod == 3){ // KEEP scheme
//ghost>                flow_float Ctilde  = 0.5*(ro[ic0]+rob[ib])*0.5*( (Ux[ic0]+Uxb[ib])*sv[0]
//ghost>                                                                +(Uy[ic0]+Uyb[ib])*sv[1]
//ghost>                                                                +(Uz[ic0]+Uzb[ib])*sv[2] );
//ghost>                flow_float Mtildex = Ctilde*(Ux[ic0]+Uxb[ib])*0.5;
//ghost>                flow_float Mtildey = Ctilde*(Uy[ic0]+Uyb[ib])*0.5;
//ghost>                flow_float Mtildez = Ctilde*(Uz[ic0]+Uzb[ib])*0.5;
//ghost>
//ghost>                flow_float Ktilde = Ctilde*0.5*(Ux[ic0]*Uxb[ib] +Uy[ic0]*Uyb[ib] +Uz[ic0]*Uzb[ib]);
//ghost>                flow_float Itilde = Ctilde*0.5*(P[ic0]/ro[ic0] +Psb[ib]/rob[ib])/(cfg.gamma-1.0);
//ghost>
//ghost>                flow_float Gtildex = 0.5*(P[ic0]+Psb[ib])*sv[0];
//ghost>                flow_float Gtildey = 0.5*(P[ic0]+Psb[ib])*sv[1];
//ghost>                flow_float Gtildez = 0.5*(P[ic0]+Psb[ib])*sv[2];
//ghost>
//ghost>                flow_float Ptilde = 0.5*((Ux[ic0]*Psb[ib] + Uxb[ib]*P[ic0])*sv[0]
//ghost>                                        +(Uy[ic0]*Psb[ib] + Uyb[ib]*P[ic0])*sv[1]
//ghost>                                        +(Uz[ic0]*Psb[ib] + Uzb[ib]*P[ic0])*sv[2]);
//ghost>
//ghost>                res_ro[ic0]   -= Ctilde;
//ghost>                res_roUx[ic0] -= Mtildex + Gtildex;
//ghost>                res_roUy[ic0] -= Mtildey + Gtildey;
//ghost>                res_roUz[ic0] -= Mtildez + Gtildez;
//ghost>                res_roe[ic0]  -= Ktilde + Itilde + Ptilde;
//ghost>            }
//ghost>        } 
//ghost>       }
//ghost>   }

    // time integration
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        geom_float vol = msh.cells[ic].volume;

        //cout << "res_ro=" << ic << " " << res_ro[ic] << endl;

        ro[ic]   = roN[ic]   + res_ro[ic]*cfg.dt/vol;
        roUx[ic] = roUxN[ic] + res_roUx[ic]*cfg.dt/vol;
        roUy[ic] = roUyN[ic] + res_roUy[ic]*cfg.dt/vol;
        roUz[ic] = roUzN[ic] + res_roUz[ic]*cfg.dt/vol;
        roe[ic]  = roeN[ic]  + res_roe[ic]*cfg.dt/vol;

        ////Ux[ic] = UxN[ic] -dPdx[ic]/ro[ic]*cfg.dt + (-convx[ic] + diffx[ic] )*cfg.dt/vol/ro[ic];
        //Uy[ic] = UyN[ic] -dPdy[ic]/ro[ic]*cfg.dt + (-convy[ic] + diffy[ic] )*cfg.dt/vol/ro[ic];
        //Uz[ic] = UzN[ic] -dPdz[ic]/ro[ic]*cfg.dt + (-convz[ic] + diffz[ic] )*cfg.dt/vol/ro[ic];
        //printf("%d UxN=%f, dt=%f, vol=%f, ro=%f, convx=%f, diffx=%f\n", ic, UxN[ic], cfg.dt, vol, ro[ic], convx[ic], diffx[ic]);
    }

}
