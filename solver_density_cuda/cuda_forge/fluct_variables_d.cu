#include "cuda_forge/cudaConfig.cuh"
#include "fluct_variables_d.cuh"
#include "cuda_forge/random.h"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"


fluct_variables::fluct_variables() {};
fluct_variables::~fluct_variables() {
    cudaWrapper::cudaFree_wrapper(this->ubar_d);
    cudaWrapper::cudaFree_wrapper(this->sigma_x_d);
    cudaWrapper::cudaFree_wrapper(this->sigma_y_d);
    cudaWrapper::cudaFree_wrapper(this->sigma_z_d);

    cudaWrapper::cudaFree_wrapper(this->k_vec_x_d);
    cudaWrapper::cudaFree_wrapper(this->k_vec_y_d);
    cudaWrapper::cudaFree_wrapper(this->k_vec_z_d);

    cudaWrapper::cudaFree_wrapper(this->Psi_d);
};

void fluct_variables::allocVariables()
{

    this->ubar.resize   (this->knum_max);
    this->sigma_x.resize(this->knum_max);
    this->sigma_y.resize(this->knum_max);
    this->sigma_z.resize(this->knum_max);
    this->k_vec_x.resize(this->knum_max);
    this->k_vec_y.resize(this->knum_max);
    this->k_vec_z.resize(this->knum_max);
    this->Psi.resize    (this->knum_max);

    cudaMalloc((void**) &(this->ubar_d)   , this->knum_max*sizeof(flow_float));
    cudaMalloc((void**) &(this->sigma_x_d), this->knum_max*sizeof(flow_float));
    cudaMalloc((void**) &(this->sigma_y_d), this->knum_max*sizeof(flow_float));
    cudaMalloc((void**) &(this->sigma_z_d), this->knum_max*sizeof(flow_float));
    cudaMalloc((void**) &(this->k_vec_x_d), this->knum_max*sizeof(flow_float));
    cudaMalloc((void**) &(this->k_vec_y_d), this->knum_max*sizeof(flow_float));
    cudaMalloc((void**) &(this->k_vec_z_d), this->knum_max*sizeof(flow_float));
    cudaMalloc((void**) &(this->Psi_d)    , this->knum_max*sizeof(flow_float));
}

void fluct_variables::set_fluctVelocity(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var )
{
    flow_float k[knum_max];
    flow_float Ek[knum_max];
    flow_float delta_k[knum_max];

    flow_float Phi[knum_max];
    flow_float Theta[knum_max];
    flow_float Alpha[knum_max];
    flow_float Omega[knum_max];

    flow_float acore = (1.+0.5*(cfg.gamma-1.0)*Mach*Mach);
    flow_float Rgas = cfg.cp - cfg.cp/cfg.gamma;

    flow_float Ts = T0/acore;
    flow_float Ps = P0/pow(acore,(cfg.gamma)/(cfg.gamma-1.));
    flow_float rhos = Ps/Rgas/Ts;
    flow_float cc = sqrt(cfg.gamma*Ps/rhos);
    flow_float mue_ref = 1.458*pow(Ts,1.5)*1e-6/(Ts+110.4);
    flow_float uinf = Mach*cc;
    flow_float uu_ref = uinf*IntTurb/100.;

    flow_float tke_sum = 0.;
    flow_float tdf_sum = 0.;

    int i;
    flow_float dkl;

    flow_float aa = 1.453;
    flow_float bb = 0.3019;
    flow_float tke_input = 2./3.*pow(uu_ref, 2.);
    flow_float tdf_input = pow(uu_ref,3.)/ramda;
    flow_float k_kol = pow(tdf_input, 1./4.)*pow(mue_ref, -3./4.);
    flow_float ke = 9.*M_PI*aa/(55.*ramda);
    flow_float ke_max = sqrt(12./5.)*ke;

    for (i=0; i<knum_max; i++) {
        dkl = (log(k_kol) - log(ke_max))/(flow_float)(knum_max - 1);
        k[i] = exp(log(ke_max) + (flow_float)i * dkl);
        Ek[i] = aa*pow(uu_ref, 2.)/ke*pow(k[i]/ke,4.)/pow(1.+pow(k[i]/ke,2.), 17./6.)*exp(-2.*pow(k[i]/k_kol,2.));
    }

    // Calc ubar
    delta_k[0] = k[1]-k[0];
    delta_k[knum_max-1] = k[knum_max-1]-k[knum_max-2];
    for (i=0; i<knum_max-2; i++) {
        delta_k[i+1] = 0.5*(k[i+2] - k[i]);
    }

    for (i=0; i<knum_max; i++) {
        ubar[i] = sqrt(Ek[i]*delta_k[i]);
        tke_sum += pow(ubar[i],2.);
        tdf_sum += 2.*mue_ref*pow(k[i],2.)*pow(ubar[i],2.);
    }

    // generate random number
    for (i=0; i<knum_max; i++) {
        Phi[i] = Uniform()*M_PI;
        Theta[i] = acos(1.-2.*Uniform()) ;
        Alpha[i] = Uniform()*2.*M_PI;
        Psi[i] = Uniform()*M_PI;
    }

    flow_float k_vec[3];
    flow_float kmag;
    flow_float sigma[3];
    flow_float k1_prime[3];
    flow_float k2_prime[3];
    flow_float k3_prime[3];
    flow_float k1_d_prime[3];

    for (i=0; i<knum_max; i++) {
        k_vec[0] = k[i] * sin(Theta[i]) * cos(Phi[i]);
        k_vec[1] = k[i] * sin(Theta[i]) * sin(Phi[i]);
        k_vec[2] = k[i] * cos(Theta[i]);

        k_vec_x[i] = k_vec[0];
        k_vec_y[i] = k_vec[1];
        k_vec_z[i] = k_vec[2];

        kmag = sqrt(pow(k_vec[0],2.)+pow(k_vec[1],2.)+pow(k_vec[2], 2.));

        k3_prime[0] = k_vec[0]/k[i];
        k3_prime[1] = k_vec[1]/k[i];
        k3_prime[2] = k_vec[2]/k[i];

        k1_d_prime[0] = cos(Phi[i]);
        k1_d_prime[1] = sin(Phi[i]);
        k1_d_prime[2] = 0.0;

        k2_prime[0] = k3_prime[1]*k1_d_prime[2] - k2_prime[2]*k1_d_prime[1];
        k2_prime[1] = k3_prime[2]*k1_d_prime[0] - k2_prime[0]*k1_d_prime[2];
        k2_prime[2] = k3_prime[0]*k1_d_prime[1] - k2_prime[1]*k1_d_prime[0];

        k1_prime[0] = k2_prime[1]*k3_prime[2] - k2_prime[2]*k3_prime[1];
        k1_prime[1] = k2_prime[2]*k3_prime[0] - k2_prime[0]*k3_prime[2];
        k1_prime[2] = k2_prime[0]*k3_prime[1] - k2_prime[1]*k3_prime[0];

        sigma[0] = cos(Theta[i])*cos(Phi[i])*cos(Alpha[i]) - sin(Phi[i])*sin(Alpha[i]);
        sigma[1] = cos(Theta[i])*sin(Phi[i])*cos(Alpha[i]) + cos(Phi[i])*sin(Alpha[i]);
        sigma[2] = -sin(Theta[i])*cos(Alpha[i]);

        sigma_x[i] = sigma[0];
        sigma_y[i] = sigma[1];
        sigma_z[i] = sigma[2];
    }


    cudaWrapper::cudaMemcpy_H2D_wrapper(this->ubar.data() , this->ubar_d, this->ubar.size());
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->sigma_x.data() , this->sigma_x_d, this->sigma_x.size());
    cudaWrapper::cudaMemcpy_H2D_wrapper(sigma_x.data() , sigma_x_d, sigma_x.size());
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->sigma_y.data() , this->sigma_y_d, this->sigma_y.size());
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->sigma_z.data() , this->sigma_z_d, this->sigma_z.size());

    cudaWrapper::cudaMemcpy_H2D_wrapper(this->k_vec_x.data() , this->k_vec_x_d, this->k_vec_x.size());
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->k_vec_y.data() , this->k_vec_y_d, this->k_vec_y.size());
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->k_vec_z.data() , this->k_vec_z_d, this->k_vec_z.size());

    cudaWrapper::cudaMemcpy_H2D_wrapper(this->Psi.data() , this->Psi_d, this->Psi.size());
}

__global__ 
void inlet_fluctVelocity_d
( 
 // gas 
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ux0b ,
 flow_float* Uy0b ,
 flow_float* Uz0b ,
 //flow_float* Ttb ,
 //flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb ,

 flow_float time,
 int knum_max,

 flow_float* ubar,

 flow_float* sigma_x,
 flow_float* sigma_y,
 flow_float* sigma_z,

 flow_float* k_vec_x,
 flow_float* k_vec_y,
 flow_float* k_vec_z,

 flow_float* Psi,

 flow_float* px,
 flow_float* py,
 flow_float* pz
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        flow_float s_vec[3];

        s_vec[0] = px[ip] - time*Ux0b[ib];
        s_vec[1] = py[ip] - time*Uy0b[ib];
        s_vec[2] = pz[ip] - time*Uz0b[ib];

        //s_vec[0] = px[ip] - time*Uxb[ib];
        //s_vec[1] = py[ip] - time*Uyb[ib];
        //s_vec[2] = pz[ip] - time*Uzb[ib];


        flow_float U_fluc = 0.;
        flow_float V_fluc = 0.;
        flow_float W_fluc = 0.;

        for (int j=0; j<knum_max; j++) {
            flow_float k_prod = k_vec_x[j]*s_vec[0] + k_vec_y[j]*s_vec[1] + k_vec_z[j]*s_vec[2];

            U_fluc = U_fluc + 2.*ubar[j]*sigma_x[j]*cos(k_prod + Psi[j]);
            V_fluc = V_fluc + 2.*ubar[j]*sigma_y[j]*cos(k_prod + Psi[j]);
            W_fluc = W_fluc + 2.*ubar[j]*sigma_z[j]*cos(k_prod + Psi[j]);
        }

        flow_float Uxbb = Ux0b[ib]+U_fluc;
        flow_float Uybb = Uy0b[ib]+V_fluc;
        flow_float Uzbb = Uz0b[ib]+W_fluc;
        //flow_float Uxbb = Uxb[ib]+U_fluc;
        //flow_float Uybb = Uyb[ib]+V_fluc;
        //flow_float Uzbb = Uzb[ib]+W_fluc;

        //ro[ig]    = rob[ib];
        //roUx[ig]  = rob[ib]*Uxbb;
        //roUy[ig]  = rob[ib]*Uybb;
        //roUz[ig]  = rob[ib]*Uzbb;
        //roe[ig]   = Psb[ib]/(ga-1.0) + 0.5*rob[ib]*(Uxbb*Uxbb +Uybb*Uybb +Uzbb*Uzbb);
        //P[ig]     = Psb[ib];
        //Ux[ig]    = Uxbb;
        //Uy[ig]    = Uybb;
        //Uz[ig]    = Uzbb;
        //Ht[ig]    = roe[ig]/rob[ib] + Psb[ib]/rob[ib];
        //sonic[ig] = sqrt(ga*Psb[ib]/rob[ib]);

        //rob[ib]    = rob[ib];
        roUxb[ib]  = rob[ib]*Uxbb;
        roUyb[ib]  = rob[ib]*Uybb;
        roUzb[ib]  = rob[ib]*Uzbb;
        roeb[ib]   = Psb[ib]/(ga-1.0) + 0.5*rob[ib]*(Uxbb*Uxbb +Uybb*Uybb +Uzbb*Uzbb);
        //Psb[ib]    = Psb[ib];
        Tsb[ib]    = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
        Uxb[ib]    = Uxbb;
        Uyb[ib]    = Uybb;
        Uzb[ib]    = Uzbb;
 
    }
};

void inlet_fluctVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p , fluct_variables& fluct)
{
    inlet_fluctVelocity_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  , 
        var.c_d["P"]  , 
        var.c_d["Ht"]  , 
        var.c_d["sonic"]  , 

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Ux0"],
        bc.bvar_d["Uy0"],
        bc.bvar_d["Uz0"],
        //bc.bvar_d["Tt"],
        //bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"],

        cfg.totalTime,
        fluct.knum_max,        

        fluct.ubar_d,

        fluct.sigma_x_d,
        fluct.sigma_y_d,
        fluct.sigma_z_d,

        fluct.k_vec_x_d,
        fluct.k_vec_y_d,
        fluct.k_vec_z_d,

        fluct.Psi_d,

        var.p_d["pcx"],
        var.p_d["pcy"],
        var.p_d["pcz"]

    ) ;

}
