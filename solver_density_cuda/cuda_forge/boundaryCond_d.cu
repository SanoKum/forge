#include "cuda_forge/boundaryCond_d.cuh"

__global__ 
void slip_d 
( 
 // gas properties
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
 flow_float* sonic  ,

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        geom_float Un =  (sx[ip]*Ux[ic] + sy[ip]*Uy[ic] + sz[ip]*Uz[ic])/ss[ip];

        ro[ig]   = ro[ic];
        P[ig]    = P[ic];
        Ux[ig]   = Ux[ic]- 2*Un*sx[ip]/ss[ip];
        Uy[ig]   = Uy[ic]- 2*Un*sy[ip]/ss[ip];
        Uz[ig]   = Uz[ic]- 2*Un*sz[ip]/ss[ip];


        roUx[ig] = ro[ic]*(Ux[ic]- 2*Un*sx[ip]/ss[ip]);
        roUy[ig] = ro[ic]*(Uy[ic]- 2*Un*sy[ip]/ss[ip]);
        roUz[ig] = ro[ic]*(Uz[ic]- 2*Un*sz[ip]/ss[ip]);

        roe[ig]  = P[ic]/(ga-1.0) + 0.5*ro[ic]*(Ux[ic]*Ux[ic]+Uy[ic]*Uy[ic]+Uz[ic]*Uz[ic]);

        Ht[ig]   = (roe[ig] + P[ig])/ro[ig];
        sonic[ig]= sqrt(ga*P[ig]/ro[ig]);


        rob[ib]   = ro[ic];
        Psb[ib]   = P[ic];
        Uxb[ib]   = Ux[ic]- Un*sx[ip]/ss[ip];
        Uyb[ib]   = Uy[ic]- Un*sy[ip]/ss[ip];
        Uzb[ib]   = Uz[ic]- Un*sz[ip]/ss[ip];

        roUxb[ib] = ro[ic]*(Ux[ic]- Un*sx[ip]/ss[ip]);
        roUyb[ib] = ro[ic]*(Uy[ic]- Un*sy[ip]/ss[ip]);
        roUzb[ib] = ro[ic]*(Uz[ic]- Un*sz[ip]/ss[ip]);

        roeb[ib]  = P[ic]/(ga-1.0) + 0.5*ro[ic]*(Ux[ic]*Ux[ic]+Uy[ic]*Uy[ic]+Uz[ic]*Uz[ic]);
        Tsb[ib]   = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
    }
};

void slip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    slip_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    );
}

__global__ 
void wall_d 
( 

 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 geom_float* fx  , 

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
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 

)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];


        ro[ig]   = ro[ic];
        Ux[ig]   = -roUx[ic]/ro[ic];
        Uy[ig]   = -roUy[ic]/ro[ic];
        Uz[ig]   = -roUz[ic]/ro[ic];
        roUx[ig] = -roUx[ic];
        roUy[ig] = -roUy[ic];
        roUz[ig] = -roUz[ic];
        roe[ig]  = roe[ic];

        flow_float ek = 0.5*(Ux[ig]*Ux[ig] +Uy[ig]*Uy[ig] +Uz[ig]*Uz[ig]);
        P[ig] =(ga-1.0)*(roe[ig]-ro[ig]*ek);

        Ht[ig]   = ga*roe[ig]/ro[ig] + (1.0-ga)*ek;
        sonic[ig]= sqrt(ga*P[ig]/ro[ig]);

        flow_float Ux_b = Uxb[ib];
        flow_float Uy_b = Uyb[ib];
        flow_float Uz_b = Uzb[ib];
        ek = 0.5*(Ux_b*Ux_b +Uy_b*Uy_b +Uz_b*Uz_b);

        flow_float R = (ga-1.0)/ga*cp;
        flow_float Tc = P[ic]/(ro[ic]*R);
        Tsb[ib]   = Tc;
        Psb[ib]   = P[ic];
        rob[ib]   = Psb[ib]/(R*Tc);
        roUxb[ib] = rob[ib]*Ux_b;
        roUyb[ib] = rob[ib]*Uy_b;
        roUzb[ib] = rob[ib]*Uz_b;
    }
};

void wall_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    wall_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  
        var.p_d["fx"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}


__global__ 
void wall_isothermal_d 
( 

 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 geom_float* fx  , 

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
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 

)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        flow_float cv = cp/ga;
        flow_float R = cp-cv;
        flow_float inE = cv*Tsb[ib];

        Psb[ib] = P[ic];
        rob[ib] = Psb[ib]/(R*Tsb[ib]);

        ro[ig]   = rob[ib];
        Ux[ig]   = -roUx[ic]/ro[ic];
        Uy[ig]   = -roUy[ic]/ro[ic];
        Uz[ig]   = -roUz[ic]/ro[ic];
        roUx[ig] = -rob[ib]*Ux[ig];
        roUy[ig] = -rob[ib]*Uy[ig];
        roUz[ig] = -rob[ib]*Uz[ig];
        roe[ig]  = roe[ic];
//
//        flow_float ek = 0.5*(Ux[ig]*Ux[ig] +Uy[ig]*Uy[ig] +Uz[ig]*Uz[ig]);
//        P[ig] =(ga-1.0)*(roe[ig]-ro[ig]*ek);
//        Ht[ig]   = (roe[ig] + P[ig])/ro[ig];
//        sonic[ig]= sqrt(ga*P[ig]/ro[ig]);
//
//        rob[ib]   = ro[ic];
//        //Uxb[ib]   = roUx[ic]/ro[ic];
//        //Uyb[ib]   = roUy[ic]/ro[ic];
//        //Uzb[ib]   = roUz[ic]/ro[ic];
//        Uxb[ib]   = 0.0;
//        Uyb[ib]   = 0.0;
//        Uzb[ib]   = 0.0;
//        roUxb[ib] = 0.0;
//        roUyb[ib] = 0.0;
//        roUzb[ib] = 0.0;
//        roeb[ib]  = roe[ic];
//
//        //flow_float ek = 0.5*(Uxb[ib]*Uxb[ib] +Uyb[ib]*Uyb[ib] +Uzb[ib]*Uzb[ib]);
//        Psb[ib] =(ga-1.0)*(roeb[ib]-rob[ib]*ek);
//        //Tsb[ib] =Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
//
        flow_float Ux_b = Uxb[ib];
        flow_float Uy_b = Uyb[ib];
        flow_float Uz_b = Uzb[ib];
        flow_float ek = 0.5*(Ux_b*Ux_b +Uy_b*Uy_b +Uz_b*Uz_b);

        roUxb[ib] = rob[ib]*Ux_b;
        roUyb[ib] = rob[ib]*Uy_b;
        roUzb[ib] = rob[ib]*Uz_b;
        Psb[ib]   = P[ic];
        rob[ib]   = Psb[ib]/(R*Tsb[ib]);
        roeb[ib]  = Psb[ib]/(ga-1.0) + rob[ib]*ek;
    }
};

void wall_isothermal_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    wall_isothermal_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  
        var.p_d["fx"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}



__global__ 
void outlet_statPress_d 
( 
 // gas properties
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
 flow_float* sonic ,

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 


)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        flow_float sxx = sx[ip];
        flow_float syy = sy[ip];
        flow_float szz = sz[ip];
        flow_float sss = ss[ip];

        flow_float Pnew = Psb[ib];
        flow_float ronew = ro[ic];
        flow_float Uxnew = Ux[ic];
        flow_float Uynew = Uy[ic];
        flow_float Uznew = Uz[ic];
        flow_float Umagc = sqrt(Uxnew*Uxnew+Uynew*Uynew+Uznew*Uznew);
        flow_float roec = roe[ic];
        flow_float Umag_new;

        flow_float Un = (Uxnew*sxx+Uynew*syy+Uznew*szz)/sss;

        if (Un >= 0) {
            Umag_new = Umagc;
            if (Un/sonic[ig]>1.0) {
                Pnew = P[ic];
            }
        } else { // backflow
            Uxnew = -Uxnew*sxx/sss;
            Uynew = -Uynew*syy/sss;
            Uznew = -Uznew*szz/sss;

            Umag_new = sqrt(Uxnew*Uxnew +Uynew*Uynew +Uznew*Uznew);

            flow_float mach_new = Umag_new/sonic[ic];

            flow_float Pt_b = Ptb[ib];
            flow_float Tt_b = Ttb[ib];
            Pnew = Pt_b/pow(1.0+0.5*(ga-1.0)*mach_new*mach_new, ga/(ga-1.0));
            flow_float Tnew = Tt_b/(1.0+0.5*(ga-1.0)*mach_new*mach_new);
            ronew = ga*Pnew/((ga-1.0)*cp*Tnew);
        }

        ro[ig]    = ronew;
        P[ig]     = Pnew;
        Ux[ig]    = Uxnew*Umag_new/Umagc;
        Uy[ig]    = Uynew*Umag_new/Umagc;
        Uz[ig]    = Uznew*Umag_new/Umagc;
        roUx[ig]  = ronew*Uxnew*Umag_new/Umagc;
        roUy[ig]  = ronew*Uynew*Umag_new/Umagc;
        roUz[ig]  = ronew*Uznew*Umag_new/Umagc;
        roe[ig]   = roec;
        Ht[ig]    = roe[ig]/ronew + Pnew/ronew;
        sonic[ig] = sqrt(ga*Pnew/ronew);

        rob[ib]    = ronew;
        //Psb[ib]    = Pnew;
        Uxb[ib]    = Uxnew*Umag_new/Umagc;
        Uyb[ib]    = Uynew*Umag_new/Umagc;
        Uzb[ib]    = Uznew*Umag_new/Umagc;
        roUxb[ib]  = ronew*Uxnew*Umag_new/Umagc;
        roUyb[ib]  = ronew*Uynew*Umag_new/Umagc;
        roUzb[ib]  = ronew*Uznew*Umag_new/Umagc;
        roeb[ib]   = roec;
        Tsb[ib]    = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
    }
};

void outlet_statPress_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    outlet_statPress_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}

__global__ 
void inlet_uniformVelocity_d
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
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 


)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        //geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        ro[ig]    = rob[ib];
        roUx[ig]  = rob[ib]*Uxb[ib];
        roUy[ig]  = rob[ib]*Uyb[ib];
        roUz[ig]  = rob[ib]*Uzb[ib];
        roe[ig]   = Psb[ib]/(ga-1.0) + 0.5*rob[ib]*(Uxb[ib]*Uxb[ib] +Uyb[ib]*Uyb[ib] +Uzb[ib]*Uzb[ib]);
        P[ig]     = Psb[ib];
        Ux[ig]    = Uxb[ib];
        Uy[ig]    = Uyb[ib];
        Uz[ig]    = Uzb[ib];
        Ht[ig]    = roe[ig]/rob[ib] + Psb[ib]/rob[ib];
        sonic[ig] = sqrt(ga*Psb[ib]/rob[ib]);

        //rob[ib]    = rob[ib];
        roUxb[ib]  = rob[ib]*Uxb[ib];
        roUyb[ib]  = rob[ib]*Uyb[ib];
        roUzb[ib]  = rob[ib]*Uzb[ib];
        roeb[ib]   = Psb[ib]/(ga-1.0) + 0.5*rob[ib]*(Uxb[ib]*Uxb[ib] +Uyb[ib]*Uyb[ib] +Uzb[ib]*Uzb[ib]);
        //Psb[ib]    = Psb[ib];
        Tsb[ib]    = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
        //Uxb[ib]    = Uxb[ib];
        //Uyb[ib]    = Uyb[ib];
        //Uzb[ib]    = Uzb[ib];
 
    }
};

void inlet_uniformVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    inlet_uniformVelocity_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
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
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]
    ) ;
}


__global__ 
void inlet_Pressure_d
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
 flow_float* T   ,

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        flow_float Ux_c;
        flow_float Uy_c;
        flow_float Uz_c;
        flow_float Umag_c;

        flow_float mach_c;
        flow_float mach_b;
        flow_float mach_new;

        //flow_float T_c;
        //flow_float P_c;

        flow_float Pt_b;
        flow_float Tt_b;

        flow_float Ps_new;
        flow_float Ts_new;
        flow_float ro_new;
        flow_float sonic_new;

        flow_float rf=0.1;

        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        Pt_b = Ptb[ib];
        Tt_b = Ttb[ib];

        Ux_c = Ux[ic];
        Uy_c = Uy[ic];
        Uz_c = Uz[ic];

        flow_float sxx = sx[ip];
        flow_float syy = sy[ip];
        flow_float szz = sz[ip];
        flow_float sss = ss[ip];


        //Umag_c = sqrt(Ux_c*Ux_c + Uy_c*Uy_c + Uz_c*Uz_c);
        flow_float Un_c = (Ux_c*sxx + Uy_c*syy + Uz_c*szz)/sss;
        mach_c = Un_c/sonic[ic];

        flow_float Ux_new;
        flow_float Uy_new;
        flow_float Uz_new;

        //P_c = P[ic];
        //T_c = T[ic];


        if ( Un_c < 0) { 

            //mach_b = sqrt((pow((P_c/Pt_b),-(ga-1.0)/ga) -1.0)*2.0/(ga-1.0));

            //flow_float mach_new = rf*mach_b + (1.0-rf)*mach_c;

            //Ps_new = P_c;
            ////Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_b*mach_b);
            //Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_new*mach_new);
            //sonic_new = sqrt((ga-1.0)*cp*Ts_new);
            //ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);

            ////Ux_new = -mach_b*sonic_new*sxx/sss;
            ////Uy_new = -mach_b*sonic_new*syy/sss;
            ////Uz_new = -mach_b*sonic_new*szz/sss;
            //Ux_new = -mach_new*sonic_new*sxx/sss;
            //Uy_new = -mach_new*sonic_new*syy/sss;
            //Uz_new = -mach_new*sonic_new*szz/sss;

            Ux_new = Un_c*sxx/sss;
            Uy_new = Un_c*syy/sss;
            Uz_new = Un_c*szz/sss;

            Ps_new = Pt_b/pow(1.0+0.5*(ga-1.0)*mach_c*mach_c, ga/(ga-1.0));
            Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_c*mach_c);
            sonic_new = sqrt((ga-1.0)*cp*Ts_new);
            ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);

            //Ux_new = -mach_b*sonic_new*sxx/sss;
            //Uy_new = -mach_b*sonic_new*syy/sss;
            //Uz_new = -mach_b*sonic_new*szz/sss;



        } else { // reverse
            Ux_new = Un_c*sxx/sss;
            Uy_new = Un_c*syy/sss;
            Uz_new = Un_c*szz/sss;

            Ps_new = Pt_b/(pow(1.0+0.5*(ga-1.0)*mach_c*mach_c, ga/(ga-1.0)));
            Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_c*mach_c);
            sonic_new = sqrt((ga-1.0)*cp*Ts_new);
            ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);
        }

        ro[ig]    = ro_new;
        Ux[ig]    = Ux_new;
        Uy[ig]    = Uy_new;
        Uz[ig]    = Uz_new;
        roUx[ig]  = ro_new*Ux_new;
        roUy[ig]  = ro_new*Uy_new;
        roUz[ig]  = ro_new*Uz_new;
        roe[ig]   = Ps_new/(ga-1.0) + 0.5*ro_new*(Ux_new*Ux_new +Uy_new*Uy_new +Uz_new*Uz_new);
        P[ig]     = Ps_new;
        Ht[ig]    = roe[ig]/ro_new + Ps_new/ro_new;
        sonic[ig] = sqrt(ga*Ps_new/ro_new);

        rob[ib]    = ro_new;
        Uxb[ib]    = Ux_new;
        Uyb[ib]    = Uy_new;
        Uzb[ib]    = Uz_new;
        roUxb[ib]  = ro_new*Ux_new;
        roUyb[ib]  = ro_new*Uy_new;
        roUzb[ib]  = ro_new*Uz_new;
        flow_float ek = 0.5*(Ux_new*Ux_new +Uy_new*Uy_new +Uz_new*Uz_new);
        roeb[ib]   = Ps_new/(ga-1.0) + ro_new*ek;
        Psb[ib]    = Ps_new;
        Tsb[ib]    = Ts_new;

        //Htb[ib]    = roeb[ib]/ro_new + Ps_new/ro_new;
        //sonicb[ib] = sqrt(ga*Ps_new/ro_new);
 
    }
};

void inlet_Pressure_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    inlet_Pressure_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
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
        var.c_d["T"]  , 

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]
    ) ;
}

__global__ 
void inlet_Pressure_dir_d
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
 flow_float* T   ,

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 

)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        flow_float Ux_c;
        flow_float Uy_c;
        flow_float Uz_c;
        flow_float Umag_c;

        flow_float mach_c;
        flow_float mach_b;
        flow_float mach_new;

        flow_float T_c;
        flow_float P_c;

        flow_float Pt_b;
        flow_float Tt_b;

        flow_float Ps_new;
        flow_float Ts_new;
        flow_float ro_new;
        flow_float sonic_new;

        flow_float rf=0.5;

        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        Pt_b = Ptb[ib];
        Tt_b = Ttb[ib];

        flow_float Ubmag = sqrt(Uxb[ib]*Uxb[ib] + Uyb[ib]*Uyb[ib] + Uzb[ib]*Uzb[ib]);

        flow_float Unx_b = Uxb[ib]/Ubmag;
        flow_float Uny_b = Uyb[ib]/Ubmag;
        flow_float Unz_b = Uzb[ib]/Ubmag;

        Ux_c = Ux[ic];
        Uy_c = Uy[ic];
        Uz_c = Uz[ic];

        Umag_c = sqrt(Ux_c*Ux_c + Uy_c*Uy_c + Uz_c*Uz_c);
        mach_c = Umag_c/sonic[ic];

        T_c = T[ic];
        P_c = P[ic];

        mach_b = sqrt((pow((P_c/Pt_b),-(ga-1.0)/ga) -1.0)*2.0/(ga-1.0));

        mach_new = rf*mach_b + (1.0-rf)*mach_c;

        Ps_new = P_c;
        Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_new*mach_new);
        sonic_new = sqrt((ga-1.0)*cp*Ts_new);
        ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);

        flow_float Ux_new = mach_new*sonic_new*Unx_b;
        flow_float Uy_new = mach_new*sonic_new*Uny_b;
        flow_float Uz_new = mach_new*sonic_new*Unz_b;

        ro[ig]    = ro_new;
        Ux[ig]    = Ux_new;
        Uy[ig]    = Uy_new;
        Uz[ig]    = Uz_new;
        roUx[ig]  = ro_new*Ux_new;
        roUy[ig]  = ro_new*Uy_new;
        roUz[ig]  = ro_new*Uz_new;
        roe[ig]   = Ps_new/(ga-1.0) + 0.5*ro_new*(Ux_new*Ux_new +Uy_new*Uy_new +Uz_new*Uz_new);
        P[ig]     = Ps_new;
        Ht[ig]    = roe[ig]/ro_new + Ps_new/ro_new;
        sonic[ig] = sqrt(ga*Ps_new/ro_new);

        rob[ib]    = ro_new;
        Uxb[ib]    = Ux_new;
        Uyb[ib]    = Uy_new;
        Uzb[ib]    = Uz_new;
        roUxb[ib]  = ro_new*Ux_new;
        roUyb[ib]  = ro_new*Uy_new;
        roUzb[ib]  = ro_new*Uz_new;
        roeb[ib]   = Ps_new/(ga-1.0) + 0.5*ro_new*(Ux_new*Ux_new +Uy_new*Uy_new +Uz_new*Uz_new);
        Psb[ib]    = Ps_new;
        Tsb[ib]    = Ps_new*ga/(ro_new*(ga-1.0)*cp);
    }
};

void inlet_Pressure_dir_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    inlet_Pressure_dir_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
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
        var.c_d["T"]  , 

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}


__global__ 
void outflow_d 
( 
 // gas properties
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
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        flow_float Pnew ;
        flow_float ronew ;
        flow_float Uxnew = Ux[ic];
        flow_float Uynew = Uy[ic];
        flow_float Uznew = Uz[ic];
        flow_float Umagc = sqrt(Uxnew*Uxnew+Uynew*Uynew+Uznew*Uznew);
        flow_float roec = roe[ic];
        flow_float Umag_new;

        flow_float sxx = sx[ip];
        flow_float syy = sy[ip];
        flow_float szz = sz[ip];
        flow_float sss = ss[ip];

        flow_float Un = (Ux[ic]*sxx+Uy[ic]*syy+Uz[ic]*szz)/sss;

        if (Un <= 0) {
            Uxnew = -Uxnew*sxx/sss;
            Uynew = -Uynew*syy/sss;
            Uznew = -Uznew*szz/sss;

            Umag_new = sqrt(Uxnew*Uxnew +Uynew*Uynew +Uznew*Uznew);

            flow_float mach_new = Umag_new/sonic[ic];

            flow_float Pt_b = Ptb[ib];
            flow_float Tt_b = Ttb[ib];
            Pnew = Pt_b/pow(1.0+0.5*(ga-1.0)*mach_new*mach_new, ga/(ga-1.0));
            flow_float Tnew = Tt_b/(1.0+0.5*(ga-1.0)*mach_new*mach_new);
            ronew = ga*Pnew/((ga-1.0)*cp*Tnew);
        }


        ro[ig]   = ro[ic];
        P[ig]    = P[ic];
        Ux[ig]   = Ux[ic];
        Uy[ig]   = Uy[ic];
        Uz[ig]   = Uz[ic];
        roUx[ig] = roUx[ic];
        roUy[ig] = roUy[ic];
        roUz[ig] = roUz[ic];
        roe[ig]  = roe[ic];

        Ht[ig]   = (roe[ig] + P[ig])/ro[ig];
        sonic[ig]= sqrt(ga*P[ig]/ro[ig]);

        rob[ib]   = ro[ic];
        Psb[ib]   = P[ic];
        Uxb[ib]   = Ux[ic];
        Uyb[ib]   = Uy[ic];
        Uzb[ib]   = Uz[ic];
        roUxb[ib] = roUx[ic];
        roUyb[ib] = roUy[ic];
        roUzb[ib] = roUz[ic];
        roeb[ib]  = roe[ic];
        Tsb[ib]   = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
    }
};

void outflow_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    outflow_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]
    );
}

__global__ 
void periodic_d 
( 
 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_int* bplane_partnerCellID,  
 geom_float dtheta,

 geom_float* x   ,  geom_float* y   ,  geom_float* z ,
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
 flow_float* T   ,
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
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb 
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip         = bplane_plane[ib];
        geom_int  ic         = bplane_cell[ib];
        geom_int  ic_partner = bplane_partnerCellID[ib];
        geom_int  ig         = bplane_cell_ghst[ib];

        ro[ig]   = ro[ic];
        P[ig]    = P[ic];
        Ux[ig]   = Ux[ic];
        Uy[ig]   = Uy[ic]*cos(-dtheta) -Uz[ic]*sin(-dtheta);
        Uz[ig]   = Uz[ic]*sin(-dtheta) +Uz[ic]*cos(-dtheta);
        roUx[ig] = roUx[ic];
        roUy[ig] = roUy[ic]*cos(-dtheta) -roUz[ic]*sin(-dtheta);
        roUz[ig] = roUz[ic]*sin(-dtheta) +roUz[ic]*cos(-dtheta);
        roe[ig]  = roe[ic];
        T[ig]    = P[ic]*ga/(ro[ig]*(ga-1.0)*cp);

        Ht[ig]   = (roe[ig] + P[ig])/ro[ig];
        sonic[ig]= sqrt(ga*P[ig]/ro[ig]);

        rob[ib]   = ro[ic];
        Psb[ib]   = P[ic];
        Uxb[ib]   = Ux[ic];
        Uyb[ib]   = Uy[ic]*cos(-dtheta) -Uz[ic]*sin(-dtheta);
        Uzb[ib]   = Uz[ic]*sin(-dtheta) +Uz[ic]*cos(-dtheta);
        roUxb[ib] = roUx[ic];
        roUyb[ib] = roUy[ic]*cos(-dtheta) -roUz[ic]*sin(-dtheta);
        roUzb[ib] = roUz[ic]*sin(-dtheta) +roUz[ic]*cos(-dtheta);
        roeb[ib]  = roe[ic];
        Tsb[ib]   = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
    }
};

void periodic_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    periodic_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        bc.bint_d["partnerCellID"],
        bc.inputInts["dtheta"],

        var.p_d["x"],  
        var.p_d["y"],  
        var.p_d["z"],  
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["T"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]
    );
}



__global__ 
void copyBcondsGradient_d 
( 
 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 // variables
 flow_float* dTdx,
 flow_float* dTdy,
 flow_float* dTdz,

 flow_float* dHtdx,
 flow_float* dHtdy,
 flow_float* dHtdz,

 flow_float* drodx,
 flow_float* drody,
 flow_float* drodz,

 flow_float* dUxdx,
 flow_float* dUydx,
 flow_float* dUzdx,

 flow_float* dUxdy,
 flow_float* dUydy,
 flow_float* dUzdy,
 
 flow_float* dUxdz,
 flow_float* dUydz,
 flow_float* dUzdz,
 
 flow_float* dPdx,
 flow_float* dPdy,
 flow_float* dPdz

)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        dTdx[ig] = dTdx[ic];
        dTdy[ig] = dTdy[ic];
        dTdz[ig] = dTdz[ic];

        dHtdx[ig] = dHtdx[ic];
        dHtdy[ig] = dHtdy[ic];
        dHtdz[ig] = dHtdz[ic];

        drodx[ig] = drodx[ic];
        drody[ig] = drody[ic];
        drodz[ig] = drodz[ic];

        dUxdx[ig] = dUxdx[ic];
        dUydx[ig] = dUydx[ic];
        dUzdx[ig] = dUzdx[ic];

        dUxdy[ig] = dUxdy[ic];
        dUydy[ig] = dUydy[ic];
        dUzdy[ig] = dUzdy[ic];
 
        dUxdz[ig] = dUxdz[ic];
        dUydz[ig] = dUydz[ic];
        dUzdz[ig] = dUzdz[ic];

        dPdx[ig] = dPdx[ic];
        dPdy[ig] = dPdy[ic];
        dPdz[ig] = dPdz[ic];
    }
};

void copyBcondsGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    copyBcondsGradient_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.c_d["dTdx"] ,
        var.c_d["dTdy"] ,
        var.c_d["dTdz"] ,
        var.c_d["dHtdx"] ,
        var.c_d["dHtdy"] ,
        var.c_d["dHtdz"] ,
        var.c_d["drodx"] ,
        var.c_d["drody"] ,
        var.c_d["drodz"] ,
        var.c_d["dUxdx"] ,
        var.c_d["dUydx"] ,
        var.c_d["dUzdx"] ,
        var.c_d["dUxdy"] ,
        var.c_d["dUydy"] ,
        var.c_d["dUzdy"] ,
        var.c_d["dUxdz"] ,
        var.c_d["dUydz"] ,
        var.c_d["dUzdz"] ,
        var.c_d["dPdx"] ,
        var.c_d["dPdy"] ,
        var.c_d["dPdz"] 
    );
}



