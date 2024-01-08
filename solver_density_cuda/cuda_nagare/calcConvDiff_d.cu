#include "calcConvDiff_d.cuh"

__global__ void calcConvDiff_np_d
( 
 flow_float vis,
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* convx , flow_float* convy , flow_float* convz,
 flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 flow_float* USN ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    if (ip < nNormalPlanes) {
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float f = fx[ip];
        //flow_float Uxf = f*Ux[ic0] + (1.0-f)*Ux[ic1];
        //flow_float Uyf = f*Uy[ic0] + (1.0-f)*Uy[ic1];
        //flow_float Uzf = f*Uz[ic0] + (1.0-f)*Uz[ic1];

        flow_float Pf  = f*P [ic0] + (1.0-f)*P [ic1];
        flow_float Tf  = f*T [ic0] + (1.0-f)*T [ic1];
        flow_float rof = f*ro[ic0] + (1.0-f)*ro[ic1];

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        //flow_float US = Uxf*sxx + Uyf*syy + Uzf*szz;
        flow_float US = USN[ip];

        flow_float convx_temp  , convy_temp  , convz_temp;
        flow_float diffx_temp  , diffy_temp  , diffz_temp;
        flow_float temp_ndia_x , temp_ndia_y , temp_ndia_z;

        flow_float pcxx = pcx[ip];
        flow_float pcyy = pcy[ip];
        flow_float pczz = pcz[ip];

        flow_float ccx0= ccx[ic0];
        flow_float ccy0= ccy[ic0];
        flow_float ccz0= ccz[ic0];

        flow_float ccx1= ccx[ic1];
        flow_float ccy1= ccy[ic1];
        flow_float ccz1= ccz[ic1];

        flow_float dccvx = ccx1 - ccx0;
        flow_float dccvy = ccy1 - ccy0;
        flow_float dccvz = ccz1 - ccz0;
        flow_float dcc   = sqrt( dccvx*dccvx + dccvy*dccvy + dccvz*dccvz);
        flow_float dccv_dot_sv = dccvx*sxx + dccvy*syy + dccvz*szz;

        flow_float deltavx = dccvx*pow(sss, 2.0)/dccv_dot_sv;
        flow_float deltavy = dccvy*pow(sss, 2.0)/dccv_dot_sv;
        flow_float deltavz = dccvz*pow(sss, 2.0)/dccv_dot_sv;
        flow_float delta = dcc*pow(sss, 2.0)/dccv_dot_sv;

        flow_float sv_nodiax = sxx - deltavx;
        flow_float sv_nodiay = syy - deltavy;
        flow_float sv_nodiaz = szz - deltavz;

        // 1st upwind
        convx_temp = 0.5*rof*( (US+abs(US))*Ux[ic0] + (US-abs(US))*Ux[ic1] );
        convy_temp = 0.5*rof*( (US+abs(US))*Uy[ic0] + (US-abs(US))*Uy[ic1] );
        convz_temp = 0.5*rof*( (US+abs(US))*Uz[ic0] + (US-abs(US))*Uz[ic1] );

        temp_ndia_x = (f*dUxdx[ic0]+(1.0-f)*dUxdx[ic1])*sv_nodiax // non-diagonal term 
                    + (f*dUxdy[ic0]+(1.0-f)*dUxdy[ic1])*sv_nodiay 
                    + (f*dUxdz[ic0]+(1.0-f)*dUxdz[ic1])*sv_nodiaz;

        temp_ndia_y = (f*dUydx[ic0]+(1.0-f)*dUydx[ic1])*sv_nodiax // non-diagonal term 
                    + (f*dUydy[ic0]+(1.0-f)*dUydy[ic1])*sv_nodiay 
                    + (f*dUydz[ic0]+(1.0-f)*dUydz[ic1])*sv_nodiaz;

        temp_ndia_z = (f*dUzdx[ic0]+(1.0-f)*dUzdx[ic1])*sv_nodiax // non-diagonal term 
                    + (f*dUzdy[ic0]+(1.0-f)*dUzdy[ic1])*sv_nodiay 
                    + (f*dUzdz[ic0]+(1.0-f)*dUzdz[ic1])*sv_nodiaz;

        diffx_temp =  vis*( (Ux[ic1] - Ux[ic0])/dcc*delta + temp_ndia_x);
        diffy_temp =  vis*( (Uy[ic1] - Uy[ic0])/dcc*delta + temp_ndia_y);
        diffz_temp =  vis*( (Uz[ic1] - Uz[ic0])/dcc*delta + temp_ndia_z);

        __syncthreads();

        atomicAdd(&convx[ic0], +convx_temp);
        atomicAdd(&convy[ic0], +convy_temp);
        atomicAdd(&convz[ic0], +convz_temp);

        atomicAdd(&convx[ic1], -convx_temp);
        atomicAdd(&convy[ic1], -convy_temp);
        atomicAdd(&convz[ic1], -convz_temp);

        atomicAdd(&diffx[ic0], +diffx_temp);
        atomicAdd(&diffy[ic0], +diffy_temp);
        atomicAdd(&diffz[ic0], +diffz_temp);

        atomicAdd(&diffx[ic1], -diffx_temp);
        atomicAdd(&diffy[ic1], -diffy_temp);
        atomicAdd(&diffz[ic1], -diffz_temp);

        __syncthreads();
    }
}

__global__ void calcConvDiff_bp_d
( 
 flow_float vis,
 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* convx , flow_float* convy , flow_float* convz,
 flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 flow_float* USN ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        flow_float convx_temp  , convy_temp  , convz_temp;
        flow_float diffx_temp  , diffy_temp  , diffz_temp;

        geom_int  ip  = bplane_plane[ib];
        geom_int  ic0 = bplane_cell[ib];

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        //flow_float US = Uxp[ip]*sxx + Uyp[ip]*syy + Uzp[ip]*szz;
        flow_float US = USN[ip];

        flow_float pcxx = pcx[ip];
        flow_float pcyy = pcy[ip];
        flow_float pczz = pcz[ip];

        flow_float ccx0= ccx[ic0];
        flow_float ccy0= ccy[ic0];
        flow_float ccz0= ccz[ic0];

        // 1st upwind
        convx_temp = rop[ip]*Uxp[ip]*US;
        convy_temp = rop[ip]*Uyp[ip]*US;
        convz_temp = rop[ip]*Uzp[ip]*US;

        flow_float dn = ( (pcxx - ccx0)*sxx
                         +(pcyy - ccy0)*syy
                         +(pczz - ccz0)*szz ) /sss;

        diffx_temp =  vis*( (Uxp[ip] - Ux[ic0])/dn*sss);
        diffy_temp =  vis*( (Uyp[ip] - Uy[ic0])/dn*sss);
        diffz_temp =  vis*( (Uzp[ip] - Uz[ic0])/dn*sss);

        __syncthreads();
        atomicAdd(&convx[ic0], +convx_temp);
        atomicAdd(&convy[ic0], +convy_temp);
        atomicAdd(&convz[ic0], +convz_temp);

        atomicAdd(&diffx[ic0], +diffx_temp);
        atomicAdd(&diffy[ic0], +diffy_temp);
        atomicAdd(&diffz[ic0], +diffz_temp);
        __syncthreads();

    }
}

void calcConvDiff_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // initialize
    cudaMemset(var.c_d["convx"], 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["convy"], 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["convz"], 0.0, msh.nCells*sizeof(flow_float));

    cudaMemset(var.c_d["diffx"], 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["diffy"], 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["diffz"], 0.0, msh.nCells*sizeof(flow_float));

    // ------------------------------
    // *** sum over normal planes ***
    // ------------------------------
    calcConvDiff_np_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        cfg.visc, 
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_nplane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"]  , 
        var.p_d["USN"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"]
    ) ;

    // --------------------------------
    // *** sum over boudnary planes ***
    // --------------------------------
    // boundary conditions
    for (auto& bc : msh.bconds)
    {
        std::cout << bc.bcondKind << std::endl;
        std::cout << bc.iBPlanes[0] << " " << bc.iBPlanes.back() << std::endl;

        if (bc.valueTypes["Ux"] == -1) continue; // do nothing for slip boundaries.

        calcConvDiff_bp_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
            cfg.visc, 
            // mesh structure
            bc.iPlanes.size() , 
            bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
            msh.nCells,
            msh.nPlanes , msh.nNormalPlanes , msh.map_nplane_cells_d,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

            // basic variables
            var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
            var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
            var.p_d["T"]  , var.c_d["T"] ,
            var.p_d["Ux"] , var.c_d["Ux"] ,
            var.p_d["Uy"] , var.c_d["Uy"] ,
            var.p_d["Uz"] , var.c_d["Uz"] ,
            var.p_d["ro"] , var.c_d["ro"] ,
            var.p_d["P"]  , var.c_d["P"]  , 
            var.p_d["USN"] ,

            // gradient
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
            var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"]
        ) ;

    }

}