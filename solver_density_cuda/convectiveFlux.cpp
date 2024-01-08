#include "convectiveFlux.hpp"

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

void convectiveFlux(int iloop , solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
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

    //vector<flow_float>& dUxdx = v.c["dUxdx"];
    //vector<flow_float>& dUydx = v.c["dUydx"];
    //vector<flow_float>& dUzdx = v.c["dUzdx"];

    //vector<flow_float>& dUxdy = v.c["dUxdy"];
    //vector<flow_float>& dUydy = v.c["dUydy"];
    //vector<flow_float>& dUzdy = v.c["dUzdy"];

    //vector<flow_float>& dUxdz = v.c["dUxdz"];
    //vector<flow_float>& dUydz = v.c["dUydz"];
    //vector<flow_float>& dUzdz = v.c["dUzdz"];

    //vector<flow_float>& dPdx = v.c["dPdx"];
    //vector<flow_float>& dPdy = v.c["dPdy"];
    //vector<flow_float>& dPdz = v.c["dPdz"];

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

            Un =( Ux_a*sv[0] +Uy_a*sv[1] +Uz_a*sv[2] )/ss;

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


        } else if (cfg.convMethod == 2){ // SLAU
            flow_float ro_p= ro[ic0];
            flow_float ro_m= ro[ic1];

            flow_float ro_L= ro_p;
            flow_float ro_R= ro_m;

            flow_float u_p = Ux[ic0];
            flow_float v_p = Uy[ic0];
            flow_float w_p = Uz[ic0];
            flow_float h_p = Ht[ic0];
            flow_float u_m = Ux[ic1];
            flow_float v_m = Uy[ic1];
            flow_float w_m = Uz[ic1];
            flow_float h_m = Ht[ic1];

            flow_float P_p = P[ic0];
            flow_float P_m = P[ic1];
            flow_float Vn_p = ((Ux[ic0])*sv[0] +(Uy[ic0])*sv[1] +(Uz[ic0])*sv[2])/ss;
            flow_float Vn_m = ((Ux[ic1])*sv[0] +(Uy[ic1])*sv[1] +(Uz[ic1])*sv[2])/ss;

            flow_float roVn_p = (roUx[ic0]*sv[0] + roUy[ic0]*sv[1] + roUz[ic0]*sv[2])/ss;
            flow_float roVn_m = (roUx[ic1]*sv[0] + roUy[ic1]*sv[1] + roUz[ic1]*sv[2])/ss;

            flow_float VnL = ((Ux[ic0])*sv[0] +(Uy[ic0])*sv[1] +(Uz[ic0])*sv[2])/ss;
            flow_float VnR = ((Ux[ic1])*sv[0] +(Uy[ic1])*sv[1] +(Uz[ic1])*sv[2])/ss;

            flow_float VnL_abs = abs(VnL);
            flow_float VnR_abs = abs(VnR);

            flow_float c_hat = 0.5*(sonic[ic0] + sonic[ic1]);

            flow_float M_p = Vn_p/c_hat;
            flow_float M_m = Vn_m/c_hat;

            flow_float ro_del= ro_m - ro_p;
            flow_float Vn_del= Vn_m - Vn_p;
            flow_float P_del = P_m  - P_p;


            flow_float beta_p, beta_m;

            if (abs(M_p)>=1.0){
                beta_p = 0.5*(M_p + abs(M_p))/M_p;
            } else {
                beta_p = 0.25*pow((M_p+1.0),2.0)*(2.0-M_p);
            }

            if (abs(M_m)>=1.0){
                beta_m = 0.5*(M_m - abs(M_m))/M_m;
            } else {
                beta_m = 0.25*pow((M_m-1.0),2.0)*(2.0+M_m);
            }

            flow_float zero = 0.0;
            flow_float one  = 1.0;
            flow_float half = 0.5;

            flow_float g = -max(min(M_p,zero),-one)*min(max(M_m,zero),one);
            flow_float Vn_hat_abs   = (ro_p*abs(Vn_p) + ro_m*abs(Vn_m))/(ro_p + ro_m);
            flow_float Vn_hat_p_abs = (1.0-g)*Vn_hat_abs + g*VnL_abs;
            flow_float Vn_hat_m_abs = (1.0-g)*Vn_hat_abs + g*VnR_abs;

            flow_float M_hat = std::min(one, sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))/c_hat);
            flow_float chi = (1.0-M_hat)*(1.0-M_hat);

            flow_float p_tilde = half*(P_p+P_m) +half*(beta_p-beta_m)*(P_p-P_m)
                                 +(one-chi)*(beta_p+beta_m-one)*half*(P_m+P_p); 
                                 //+sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))*(beta_p+beta_m-one)*half*(P_m+P_p)*c_hat; 


            //flow_float mdot = ss*(0.5*(ro_L*VnL+ro_R*VnR -Vn_hat_abs*ro_del) -0.5*chi/(c_hat)*P_del);
            flow_float mdot = ss*0.5*((ro_L*(VnL+Vn_hat_p_abs)+ro_R*(VnR-Vn_hat_m_abs)) -chi/(c_hat)*P_del);

            res_ro[ic0]   -= mdot;
            res_roUx[ic0] -= 0.5*(mdot+abs(mdot))*u_p +0.5*(mdot-abs(mdot))*u_m +p_tilde*sv[0];
            res_roUy[ic0] -= 0.5*(mdot+abs(mdot))*v_p +0.5*(mdot-abs(mdot))*v_m +p_tilde*sv[1];
            res_roUz[ic0] -= 0.5*(mdot+abs(mdot))*w_p +0.5*(mdot-abs(mdot))*w_m +p_tilde*sv[2];
            res_roe[ic0]  -= 0.5*(mdot+abs(mdot))*h_p +0.5*(mdot-abs(mdot))*h_m ;

            res_ro[ic1]   += mdot;
            res_roUx[ic1] += 0.5*(mdot+abs(mdot))*u_p +0.5*(mdot-abs(mdot))*u_m +p_tilde*sv[0];
            res_roUy[ic1] += 0.5*(mdot+abs(mdot))*v_p +0.5*(mdot-abs(mdot))*v_m +p_tilde*sv[1];
            res_roUz[ic1] += 0.5*(mdot+abs(mdot))*w_p +0.5*(mdot-abs(mdot))*w_m +p_tilde*sv[2];
            res_roe[ic1]  += 0.5*(mdot+abs(mdot))*h_p +0.5*(mdot-abs(mdot))*h_m ;


        } else if (cfg.convMethod == 3){ // KEEP scheme
            flow_float nx = sv[0]/ss;
            flow_float ny = sv[1]/ss;
            flow_float nz = sv[2]/ss;

            flow_float Ctilde  = 0.5*(ro[ic0]+ro[ic1])*0.5*( (Ux[ic0]+Ux[ic1])*nx
                                                            +(Uy[ic0]+Uy[ic1])*ny
                                                            +(Uz[ic0]+Uz[ic1])*nz );
            flow_float Mtildex = Ctilde*(Ux[ic0]+Ux[ic1])*0.5;
            flow_float Mtildey = Ctilde*(Uy[ic0]+Uy[ic1])*0.5;
            flow_float Mtildez = Ctilde*(Uz[ic0]+Uz[ic1])*0.5;

            flow_float Ktilde = Ctilde*0.5*(Ux[ic0]*Ux[ic1] +Uy[ic0]*Uy[ic1] +Uz[ic0]*Uz[ic1]);
            flow_float Itilde = Ctilde*0.5*(P[ic0]/ro[ic0] +P[ic1]/ro[ic1])/(cfg.gamma-1.0);

            flow_float Gtildex = 0.5*(P[ic0]+P[ic1])*nx;
            flow_float Gtildey = 0.5*(P[ic0]+P[ic1])*ny;
            flow_float Gtildez = 0.5*(P[ic0]+P[ic1])*nz;

            flow_float Ptilde = 0.5*((Ux[ic0]*P[ic1] + Ux[ic1]*P[ic0])*nx
                                    +(Uy[ic0]*P[ic1] + Uy[ic1]*P[ic0])*ny
                                    +(Uz[ic0]*P[ic1] + Uz[ic1]*P[ic0])*nz);

            res_ro[ic0]   -= Ctilde*ss;
            res_roUx[ic0] -= (Mtildex + Gtildex)*ss;
            res_roUy[ic0] -= (Mtildey + Gtildey)*ss;
            res_roUz[ic0] -= (Mtildez + Gtildez)*ss;
            res_roe[ic0]  -= (Ktilde + Itilde + Ptilde)*ss;

            res_ro[ic1]   += Ctilde*ss;
            res_roUx[ic1] += (Mtildex + Gtildex)*ss;
            res_roUy[ic1] += (Mtildey + Gtildey)*ss;
            res_roUz[ic1] += (Mtildez + Gtildez)*ss;
            res_roe[ic1]  += (Ktilde + Itilde + Ptilde)*ss;

        } else if (cfg.convMethod == 4){ // AUSM+
            flow_float Ht_L = Ht[ic0];
            flow_float Ht_R = Ht[ic1];

            flow_float P_L = P[ic0];
            flow_float P_R = P[ic1];

            flow_float ro_L= ro[ic0];
            flow_float ro_R= ro[ic1];

            //flow_float roUx_L = roUx[ic0];
            //flow_float roUy_L = roUy[ic0];
            //flow_float roUz_L = roUz[ic0];

            //flow_float roUx_R = roUx[ic1];
            //flow_float roUy_R = roUy[ic1];
            //flow_float roUz_R = roUz[ic1];

            flow_float Ux_L = Ux[ic0];
            flow_float Uy_L = Uy[ic0];
            flow_float Uz_L = Uz[ic0];

            flow_float Ux_R = Ux[ic1];
            flow_float Uy_R = Uy[ic1];
            flow_float Uz_R = Uz[ic1];

            flow_float Vn_L = ((Ux[ic0])*sv[0] +(Uy[ic0])*sv[1] +(Uz[ic0])*sv[2])/ss; //ok
            flow_float Vn_R = ((Ux[ic1])*sv[0] +(Uy[ic1])*sv[1] +(Uz[ic1])*sv[2])/ss; //ok

            // A SEQUEL TO AUSM: AUSM1 (Liou)
            // A1
            flow_float a_star_L = sqrt(Ht_L*2.0*(cfg.gamma-1.0)/(cfg.gamma+1.0)); //ok
            flow_float a_star_R = sqrt(Ht_R*2.0*(cfg.gamma-1.0)/(cfg.gamma+1.0)); //ok
            flow_float a_tilde_L = pow(a_star_L,2.0)/max(a_star_L, abs(Vn_L)); //ok
            flow_float a_tilde_R = pow(a_star_R,2.0)/max(a_star_R, abs(Vn_R)); //ok
            flow_float a_half = min(a_tilde_L, a_tilde_R); //ok

            flow_float ML = Vn_L/a_half; //ok
            flow_float MR = Vn_R/a_half; //ok

            flow_float beta  = 1.0/8.0;
            //flow_float beta  = 0.5;
            flow_float alpha = 3.0/16.0;

            // A2 calc m(j+1/2)
            flow_float MM_pls ;
            flow_float MM_mns ;

            flow_float PP_pls;
            flow_float PP_mns;

            //flow_float MMbeta_pls = +0.5*pow((ML+1),2.0) +beta*pow((ML*ML-1),2.0); //ok
            //flow_float MMbeta_mns = -0.5*pow((MR-1),2.0) -beta*pow((MR*MR-1),2.0); //ok
            flow_float MMbeta_pls = +0.25*pow((ML+1),2.0) +beta*pow((ML*ML-1),2.0); //ok
            flow_float MMbeta_mns = -0.25*pow((MR-1),2.0) -beta*pow((MR*MR-1),2.0); //ok

            flow_float PPalpha_pls = +0.25*(ML+1)*(ML+1)*(2.0-ML) +alpha*ML*(ML*ML-1)*(ML*ML-1);//ok
            flow_float PPalpha_mns = +0.25*(MR-1)*(MR-1)*(2.0+MR) -alpha*MR*(MR*MR-1)*(MR*MR-1);//ok

            if (abs(ML)>=1.0){ 
                MM_pls = 0.5*(ML+abs(ML)); //ok
                PP_pls = 0.5*(ML+abs(ML))/ML; //ok
            } else {
                MM_pls = MMbeta_pls; //ok
                PP_pls = PPalpha_pls;//ok
            }

            if (abs(MR)>=1.0){
                MM_mns = 0.5*(MR-abs(MR)); //ok
                PP_mns = 0.5*(MR-abs(MR))/MR;//ok
            } else {
                MM_mns = MMbeta_mns; //ok
                PP_mns = PPalpha_mns;//ok
            }
           
            flow_float P_half= PP_pls*P_L + PP_mns*P_R; //ok

            flow_float m_half = MM_pls + MM_mns ; //ok
            flow_float m_half_pls = 0.5*(m_half + abs(m_half)); //ok
            flow_float m_half_mns = 0.5*(m_half - abs(m_half)); //ok

            res_ro[ic0]   -= a_half*(m_half_pls*ro_L      + m_half_mns*ro_R)*ss ;
            res_roUx[ic0] -= a_half*(m_half_pls*ro_L*Ux_L + m_half_mns*ro_R*Ux_R)*ss +P_half*sv[0] ;
            res_roUy[ic0] -= a_half*(m_half_pls*ro_L*Uy_L + m_half_mns*ro_R*Uy_R)*ss +P_half*sv[1] ;
            res_roUz[ic0] -= a_half*(m_half_pls*ro_L*Uz_L + m_half_mns*ro_R*Uz_R)*ss +P_half*sv[2] ;
            res_roe[ic0]  -= a_half*(m_half_pls*ro_L*Ht_L + m_half_mns*ro_R*Ht_R)*ss ;

            res_ro[ic1]   += a_half*(m_half_pls*ro_L      + m_half_mns*ro_R)*ss ;
            res_roUx[ic1] += a_half*(m_half_pls*ro_L*Ux_L + m_half_mns*ro_R*Ux_R)*ss +P_half*sv[0] ;
            res_roUy[ic1] += a_half*(m_half_pls*ro_L*Uy_L + m_half_mns*ro_R*Uy_R)*ss +P_half*sv[1] ;
            res_roUz[ic1] += a_half*(m_half_pls*ro_L*Uz_L + m_half_mns*ro_R*Uz_R)*ss +P_half*sv[2] ;
            res_roe[ic1]  += a_half*(m_half_pls*ro_L*Ht_L + m_half_mns*ro_R*Ht_R)*ss ;

        } else if (cfg.convMethod == 5){ // AUSM+UP
            flow_float Ht_L = Ht[ic0];
            flow_float Ht_R = Ht[ic1];

            flow_float P_L = P[ic0];
            flow_float P_R = P[ic1];

            flow_float ro_L= ro[ic0];
            flow_float ro_R= ro[ic1];

            flow_float Ux_L = Ux[ic0];
            flow_float Uy_L = Uy[ic0];
            flow_float Uz_L = Uz[ic0];

            flow_float Ux_R = Ux[ic1];
            flow_float Uy_R = Uy[ic1];
            flow_float Uz_R = Uz[ic1];

            flow_float Vn_L = ((Ux[ic0])*sv[0] +(Uy[ic0])*sv[1] +(Uz[ic0])*sv[2])/ss; //ok
            flow_float Vn_R = ((Ux[ic1])*sv[0] +(Uy[ic1])*sv[1] +(Uz[ic1])*sv[2])/ss; //ok

            // A SEQUEL TO AUSM: AUSM1 (Liou)
            // A1
            flow_float a_star_L = sqrt(Ht_L*2.0*(cfg.gamma-1.0)/(cfg.gamma+1.0)); //ok
            flow_float a_star_R = sqrt(Ht_R*2.0*(cfg.gamma-1.0)/(cfg.gamma+1.0)); //ok
            flow_float a_tilde_L = pow(a_star_L,2.0)/max(a_star_L, abs(Vn_L)); //ok
            flow_float a_tilde_R = pow(a_star_R,2.0)/max(a_star_R, abs(Vn_R)); //ok
            flow_float a_half = min(a_tilde_L, a_tilde_R); //ok

            flow_float ML = Vn_L/a_half; //ok
            flow_float MR = Vn_R/a_half; //ok

            flow_float beta  = 1.0/8.0;
            flow_float alpha = 3.0/16.0;

            flow_float Mp_L_1 = 0.5*(ML+abs(ML)); //ok
            flow_float Mp_R_1 = 0.5*(MR+abs(MR)); //ok

            flow_float Mm_L_1 = 0.5*(ML-abs(ML)); //ok
            flow_float Mm_R_1 = 0.5*(MR-abs(MR)); //ok

            flow_float Mp_R_2 = 0.25*pow(MR+1.0, 2.0); //ok
            flow_float Mp_L_2 = 0.25*pow(ML+1.0, 2.0); //ok

            flow_float Mm_L_2 =-0.25*pow(ML-1.0, 2.0); //ok
            flow_float Mm_R_2 =-0.25*pow(MR-1.0, 2.0); //ok

            flow_float Mp_L_4;
            flow_float Mm_R_4;


            if (abs(ML)>=1) {
                Mp_L_4 = Mp_L_1;
            } else {
                Mp_L_4 = Mp_L_2*(1.0-16*beta*Mm_L_2);
            }

            if (abs(MR)>=1) {
                Mm_R_4 = Mm_R_1;
            } else {
                Mm_R_4 = Mm_R_2*(1.0+16*beta*Mp_R_2);
            }

            flow_float M_half = Mp_L_4 + Mm_R_4;


            flow_float Pp_L_5;
            flow_float Pp_R_5;

            flow_float Pm_L_5;
            flow_float Pm_R_5;


            if (abs(ML)>=1.0){ 
                Pp_L_5 = Mp_L_1/ML;
                Pm_L_5 = Mm_L_1/ML;
            } else {
                Pp_L_5 = Mp_L_2*((+2.0-ML)-16.0*alpha*ML*Mm_L_2);
                Pm_L_5 = Mm_L_1*((-2.0-ML)+16.0*alpha*ML*Mp_L_2);
            }

            if (abs(MR)>=1.0){ 
                Pp_R_5 = Mp_R_1/MR;
                Pm_R_5 = Mm_R_1/MR;
            } else {
                Pp_R_5 = Mp_R_1*((+2.0-MR)-16.0*alpha*MR*Mm_R_2);
                Pm_R_5 = Mm_R_2*((-2.0-MR)+16.0*alpha*MR*Mp_R_2);
            }

            flow_float Kp    = 0.5; // 0<Kp<1
            flow_float Ku    = 0.5; // 0<Ku<1
            flow_float sigma = 0.5; // ?
            flow_float Mbar = sqrt(0.5*(ML*ML + MR*MR));
            flow_float Mp = -Kp*max(1.0-sigma*Mbar*Mbar, 0.0)*(P_R-P_L)/(0.5*(ro_L+ro_R))/(a_half*a_half);
            flow_float Pu = -Ku*Pp_L_5*Pm_R_5*(ro_L+ro_R)*a_half*(Vn_R-Vn_L);
          
            M_half += Mp;
            flow_float P_half= Pp_L_5*P_L + Pm_R_5*P_R + Pu; //ok

            flow_float M_half_pls = 0.5*(M_half + abs(M_half)); //ok
            flow_float M_half_mns = 0.5*(M_half - abs(M_half)); //ok

            res_ro[ic0]   -= a_half*(M_half_pls*ro_L      + M_half_mns*ro_R)*ss ;
            res_roUx[ic0] -= a_half*(M_half_pls*ro_L*Ux_L + M_half_mns*ro_R*Ux_R)*ss +P_half*sv[0] ;
            res_roUy[ic0] -= a_half*(M_half_pls*ro_L*Uy_L + M_half_mns*ro_R*Uy_R)*ss +P_half*sv[1] ;
            res_roUz[ic0] -= a_half*(M_half_pls*ro_L*Uz_L + M_half_mns*ro_R*Uz_R)*ss +P_half*sv[2] ;
            res_roe[ic0]  -= a_half*(M_half_pls*ro_L*Ht_L + M_half_mns*ro_R*Ht_R)*ss ;

            res_ro[ic1]   += a_half*(M_half_pls*ro_L      + M_half_mns*ro_R)*ss ;
            res_roUx[ic1] += a_half*(M_half_pls*ro_L*Ux_L + M_half_mns*ro_R*Ux_R)*ss +P_half*sv[0] ;
            res_roUy[ic1] += a_half*(M_half_pls*ro_L*Uy_L + M_half_mns*ro_R*Uy_R)*ss +P_half*sv[1] ;
            res_roUz[ic1] += a_half*(M_half_pls*ro_L*Uz_L + M_half_mns*ro_R*Uz_R)*ss +P_half*sv[2] ;
            res_roe[ic1]  += a_half*(M_half_pls*ro_L*Ht_L + M_half_mns*ro_R*Ht_R)*ss ;
        }


    }
}
