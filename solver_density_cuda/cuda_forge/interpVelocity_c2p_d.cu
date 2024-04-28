#include "interpVelocity_c2p_d.cuh"

#define CHECK_CUDA_ERROR(val) check((val), #val, __FILE__, __LINE__)
template <typename T>
void check(T err, const char* const func, const char* const file,
           const int line)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        // std::exit(EXIT_FAILURE);
    }
}


__global__ void interpVelocity_c2p_d
( 
 flow_float dt,
 // mesh structure
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* nplane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* Ux ,
 flow_float* Uy ,
 flow_float* Uz ,
 flow_float* ro ,
 flow_float* P  ,
 
 flow_float* divU_vol ,
 
 flow_float* USp,

 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    if (ip < nNormalPlanes) {
        geom_int  ic0 = nplane_cells[2*ip+0];
        geom_int  ic1 = nplane_cells[2*ip+1];

        geom_float f   = fx[ip];

        //flow_float Pf  = f*P [ic0] + (1.0-f)*P [ic1];
        //flow_float Tf  = f*T [ic0] + (1.0-f)*T [ic1];
        flow_float rof = f*ro[ic0] + (1.0-f)*ro[ic1];

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

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


        flow_float Uxf = f*(Ux[ic0] + dPdx[ic0]*dt/ro[ic0]) +(1.0-f)*(Ux[ic1] + dPdx[ic1]*dt/ro[ic1]);
        flow_float Uyf = f*(Uy[ic0] + dPdy[ic0]*dt/ro[ic0]) +(1.0-f)*(Uy[ic1] + dPdy[ic1]*dt/ro[ic1]);
        flow_float Uzf = f*(Uz[ic0] + dPdz[ic0]*dt/ro[ic0]) +(1.0-f)*(Uz[ic1] + dPdz[ic1]*dt/ro[ic1]);

        flow_float US = Uxf*sxx + Uyf*syy + Uzf*szz;

        flow_float temp_ndia = (f*dPdx[ic0]+(1.0-f)*dPdx[ic1])*sv_nodiax   // non-diagonal term 
                             + (f*dPdy[ic0]+(1.0-f)*dPdy[ic1])*sv_nodiay
                             + (f*dPdz[ic0]+(1.0-f)*dPdz[ic1])*sv_nodiaz ;

        US += (-(P[ic1] - P[ic0])/dcc*delta - temp_ndia)*dt/rof;

        USp[ip] = US;

        __syncthreads();
        atomicAdd(&divU_vol[ic0], US);
        atomicAdd(&divU_vol[ic1],-US);
        __syncthreads();
    }
}



__global__ void interpVelocity_c2p_outlet_d
( 
 flow_float dt,
 // mesh structure
 geom_int nb, 
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* nplane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* Ux ,
 flow_float* Uy ,
 flow_float* Uz ,
 flow_float* ro ,
 flow_float* P  ,
 
 flow_float* divU_vol ,
 
 flow_float* USp,
 flow_float* Pp,

 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz

)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        __syncthreads();
        geom_int  ip  = bplane_plane[ib];
        geom_int  ic0 = bplane_cell[ib];

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float pcxx = pcx[ip];
        flow_float pcyy = pcy[ip];
        flow_float pczz = pcz[ip];

        flow_float ccx0= ccx[ic0];
        flow_float ccy0= ccy[ic0];
        flow_float ccz0= ccz[ic0];

        flow_float dn = ( (pcxx - ccx0)*sxx
                         +(pcyy - ccy0)*syy 
                         +(pczz - ccz0)*szz )/sss;

        flow_float US = (Ux[ic0] + dPdx[ic0]*dt/ro[ic0])*sxx  
                       +(Uy[ic0] + dPdy[ic0]*dt/ro[ic0])*syy  
                       +(Uz[ic0] + dPdz[ic0]*dt/ro[ic0])*szz;

        US += (-(Pp[ip] - P[ic0])/dn*sss)*dt/ro[ic0];
        USp[ip] = US;

        atomicAdd(&divU_vol[ic0], +US);
        __syncthreads();
        //printf("ip=%d, US=%f, dPdx=%f, sx=%f, sy=%f, sz=%f, ss=%f, dn=%f, Pp=%f, P=%f, divU_vol=%f\n", ip, US, dPdx[ic0], sxx, syy, szz, sss, dn, Pp[ip], P[ic0], divU_vol[ic0]);
        //printf("ic0=%d, ip=%d, US=%f, dPdx=%f, sx=%f, sy=%f, sz=%f, ss=%f, dn=%f, Pp=%f, P=%f, divU_vol=%f\n", ic0, ip, US, dPdx[ic0], sxx, syy, szz, sss, dn, Pp[ip], P[ic0], divU_vol[ic0]);
    }
}

__global__ void interpVelocity_c2p_inlet_d
( 
 flow_float dt,
 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* nplane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* Uxp ,
 flow_float* Uyp ,
 flow_float* Uzp ,
 
 flow_float* divU_vol ,
 
 flow_float* USp
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip  = bplane_plane[ib];
        geom_int  ic0 = bplane_cell[ib];

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float US = Uxp[ip]*sxx + Uyp[ip]*syy + Uzp[ip]*szz;

        USp[ip] = US;

        __syncthreads();
        atomicAdd(&divU_vol[ic0], US);
        __syncthreads();
    }
}


void interpVelocity_c2p_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    cudaMemset(var.c_d["divU*vol"] , 0.0, msh.nCells *sizeof(flow_float));
    cudaMemset(var.p_d["US"]       , 0.0, msh.nPlanes*sizeof(flow_float));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
    
    // rhi-chow interpolation
    interpVelocity_c2p_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        cfg.dt, 
        // mesh structure
        msh.nPlanes , msh.nNormalPlanes , msh.map_nplane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["ro"] ,
        var.c_d["P"] ,

        var.c_d["divU*vol"]  , 

        var.p_d["US"] ,

        var.c_d["dPdx"] , var.c_d["dPdy"] , var.c_d["dPdz"] 
    ) ;

    // --------------------------------
    // *** sum over boudnary planes ***
    // --------------------------------
    // boundary conditions
    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "outlet_statPress")
        {
            interpVelocity_c2p_outlet_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
                cfg.dt, 
                // mesh structure
                bc.iPlanes.size() ,
                bc.map_bplane_plane_d,  
                bc.map_bplane_cell_d,
                msh.nPlanes , msh.nNormalPlanes , msh.map_nplane_cells_d,
                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

                // basic variables
                var.c_d["Ux"] ,
                var.c_d["Uy"] ,
                var.c_d["Uz"] ,
                var.c_d["ro"] ,
                var.c_d["P"]  , 

                var.c_d["divU*vol"],
                var.p_d["US"]  ,
                var.p_d["P"]  ,
                var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"]
            ) ;

            //forprint = (geom_float*) malloc(msh.nCells*sizeof(geom_float));

            //cudaMemcpy(var.c_d["divU*vol"] , forprint , msh.nCells*sizeof(geom_float) , cudaMemcpyDeviceToHost);
            //printf("after outlet_statPress\n");
            //printf("32087 divU=%f\n", forprint[32087]);

            //free(forprint);


        } else if (bc.bcondKind == "inlet_uniformVelocity")
        {
            interpVelocity_c2p_inlet_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
                cfg.dt, 
                // mesh structure
                bc.iPlanes.size() ,
                bc.map_bplane_plane_d,  
                bc.map_bplane_cell_d,
                msh.nPlanes , msh.nNormalPlanes , msh.map_nplane_cells_d,
                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

                // basic variables
                var.p_d["Ux"] ,
                var.p_d["Uy"] ,
                var.p_d["Uz"] ,

                var.c_d["divU*vol"],
                var.p_d["US"]
            ) ;
        }
    }
}