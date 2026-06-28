#include "calcGradient_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

__global__ void calcGradient_1_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* roe  ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx , flow_float* droUxdy, flow_float* droUxdz,
 //flow_float* droUydx , flow_float* droUydy, flow_float* droUydz,
 //flow_float* droUzdx , flow_float* droUzdy, flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    __syncthreads();

    //if (ip < nNormalPlanes) {
    if (ip < nPlanes) {
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float f   = fx[ip];
        flow_float Uxf = f*Ux[ic0] + (1.0-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0-f)*Uz[ic1];

        flow_float Pf  = f*P[ic0]  + (1.0-f)*P[ic1];
        flow_float Tf  = f*T[ic0]  + (1.0-f)*T[ic1];
        flow_float rof = f*ro[ic0] + (1.0-f)*ro[ic1];
        flow_float roef= f*roe[ic0]+ (1.0-f)*roe[ic1];
        flow_float Htf= f*Ht[ic0]+ (1.0-f)*Ht[ic1];
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];

        flow_float US  = Uxf*sxx +Uyf*syy +Uzf*szz;

//        if (ic0 == 0 or ic1 == 0) {
//            printf("ic0 = %d, ic1 = %d, f = %f, Uxf = %f, Uyf = %f, Uzf = %f\n", ic0, ic1, f, Uxf, Uyf, Uzf);
//            printf("sxx = %f, syy = %f, szz = %f\n", sxx, syy, szz);
//            printf("Pf  = %f\n", Pf);    
//            printf("US  = %f\n", US);    
//        }
//
//        if (ic0 == 1 or ic1 == 1) {
//            printf("ic0 = %d, ic1 = %d, f = %f, Uxf = %f, Uyf = %f, Uzf = %f\n", ic0, ic1, f, Uxf, Uyf, Uzf);
//        }



        atomicAdd(&dUxdx[ic0], sxx*Uxf);
        atomicAdd(&dUxdy[ic0], syy*Uxf);
        atomicAdd(&dUxdz[ic0], szz*Uxf);

        atomicAdd(&dUydx[ic0], sxx*Uyf);
        atomicAdd(&dUydy[ic0], syy*Uyf);
        atomicAdd(&dUydz[ic0], szz*Uyf);

        atomicAdd(&dUzdx[ic0], sxx*Uzf);
        atomicAdd(&dUzdy[ic0], syy*Uzf);
        atomicAdd(&dUzdz[ic0], szz*Uzf);

        atomicAdd(&dTdx[ic0], sxx*Tf);
        atomicAdd(&dTdy[ic0], syy*Tf);
        atomicAdd(&dTdz[ic0], szz*Tf);

        atomicAdd(&dPdx[ic0], sxx*Pf);
        atomicAdd(&dPdy[ic0], syy*Pf);
        atomicAdd(&dPdz[ic0], szz*Pf);

        //atomicAdd(&dHtdx[ic0], sxx*Htf);
        //atomicAdd(&dHtdy[ic0], syy*Htf);
        //atomicAdd(&dHtdz[ic0], szz*Htf);



        atomicAdd(&drodx[ic0], sxx*rof);
        atomicAdd(&drody[ic0], syy*rof);
        atomicAdd(&drodz[ic0], szz*rof);


        //atomicAdd(&droUxdx[ic0], sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic0], syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic0], szz*rof*Uxf);

        //atomicAdd(&droUydx[ic0], sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic0], syy*rof*Uyf);
        //atomicAdd(&droUydz[ic0], szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic0], sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic0], syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic0], szz*rof*Uzf);

        //atomicAdd(&droedx[ic0], sxx*roef);
        //atomicAdd(&droedy[ic0], syy*roef);
        //atomicAdd(&droedz[ic0], szz*roef);


        atomicAdd(&divU[ic0], US);


        atomicAdd(&dUxdx[ic1], -sxx*Uxf);
        atomicAdd(&dUxdy[ic1], -syy*Uxf);
        atomicAdd(&dUxdz[ic1], -szz*Uxf);

        atomicAdd(&dUydx[ic1], -sxx*Uyf);
        atomicAdd(&dUydy[ic1], -syy*Uyf);
        atomicAdd(&dUydz[ic1], -szz*Uyf);

        atomicAdd(&dUzdx[ic1], -sxx*Uzf);
        atomicAdd(&dUzdy[ic1], -syy*Uzf);
        atomicAdd(&dUzdz[ic1], -szz*Uzf);

        atomicAdd(&dTdx[ic1], -sxx*Tf);
        atomicAdd(&dTdy[ic1], -syy*Tf);
        atomicAdd(&dTdz[ic1], -szz*Tf);

        atomicAdd(&dPdx[ic1], -sxx*Pf);
        atomicAdd(&dPdy[ic1], -syy*Pf);
        atomicAdd(&dPdz[ic1], -szz*Pf);

        //atomicAdd(&dHtdx[ic1], -sxx*Pf);
        //atomicAdd(&dHtdy[ic1], -syy*Pf);
        //atomicAdd(&dHtdz[ic1], -szz*Pf);

        atomicAdd(&drodx[ic1], -sxx*rof);
        atomicAdd(&drody[ic1], -syy*rof);
        atomicAdd(&drodz[ic1], -szz*rof);

        //atomicAdd(&droUxdx[ic1], -sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic1], -syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic1], -szz*rof*Uxf);

        //atomicAdd(&droUydx[ic1], -sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic1], -syy*rof*Uyf);
        //atomicAdd(&droUydz[ic1], -szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic1], -sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic1], -syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic1], -szz*rof*Uzf);

        //atomicAdd(&droedx[ic1], -sxx*roef);
        //atomicAdd(&droedy[ic1], -syy*roef);
        //atomicAdd(&droedz[ic1], -szz*roef);

        atomicAdd(&divU[ic1], -US);
    }
    __syncthreads();
}

__global__ void calcGradient_b_d
( 
  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

// variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Tt   ,
 flow_float* Pt   ,
 flow_float* Ts   ,
 flow_float* Ps   ,
 flow_float* Ht   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx  , flow_float* droUxdy , flow_float* droUxdz,
 //flow_float* droUydx  , flow_float* droUydy , flow_float* droUydz,
 //flow_float* droUzdx  , flow_float* droUzdy , flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];
       
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float rof   = ro[ib];
        flow_float roUxf = roUx[ib];
        flow_float roUyf = roUy[ib];
        flow_float roUzf = roUz[ib];
        flow_float roef  = roe[ib];
        flow_float Uxf = roUxf/rof;
        flow_float Uyf = roUyf/rof;
        flow_float Uzf = roUzf/rof;
        flow_float ek = 0.5*(Uxf*Uxf + Uyf*Uyf + Uzf*Uzf);
        flow_float Pf = Ps[ib];
        flow_float Tf = Ts[ib];
        flow_float Htf= roef/rof + Pf/rof;

        flow_float US  = Uxf*sxx +Uyf*syy +Uzf*szz;

        atomicAdd(&dUxdx[ic], sxx*Uxf);
        atomicAdd(&dUxdy[ic], syy*Uxf);
        atomicAdd(&dUxdz[ic], szz*Uxf);

        atomicAdd(&dUydx[ic], sxx*Uyf);
        atomicAdd(&dUydy[ic], syy*Uyf);
        atomicAdd(&dUydz[ic], szz*Uyf);

        atomicAdd(&dUzdx[ic], sxx*Uzf);
        atomicAdd(&dUzdy[ic], syy*Uzf);
        atomicAdd(&dUzdz[ic], szz*Uzf);

        atomicAdd(&dTdx[ic], sxx*Tf);
        atomicAdd(&dTdy[ic], syy*Tf);
        atomicAdd(&dTdz[ic], szz*Tf);

        atomicAdd(&dPdx[ic], sxx*Pf);
        atomicAdd(&dPdy[ic], syy*Pf);
        atomicAdd(&dPdz[ic], szz*Pf);

        //atomicAdd(&dHtdx[ic], sxx*Htf);
        //atomicAdd(&dHtdy[ic], syy*Htf);
        //atomicAdd(&dHtdz[ic], szz*Htf);

        atomicAdd(&drodx[ic], sxx*rof);
        atomicAdd(&drody[ic], syy*rof);
        atomicAdd(&drodz[ic], szz*rof);

        //atomicAdd(&droUxdx[ic], sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic], syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic], szz*rof*Uxf);

        //atomicAdd(&droUydx[ic], sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic], syy*rof*Uyf);
        //atomicAdd(&droUydz[ic], szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic], sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic], syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic], szz*rof*Uzf);

        //atomicAdd(&droedx[ic], sxx*roef);
        //atomicAdd(&droedy[ic], syy*roef);
        //atomicAdd(&droedz[ic], szz*roef);

        atomicAdd(&divU[ic], US);
    }
    __syncthreads();
}


// K5: atomicAdd を使わない cell-gather 版（face並列+atomic を cell並列に置換）。
// 各セルが自分の面を走査し、面値 φf=f*φ[ic0]+(1-f)*φ[ic1] と外向き符号で寄与を集約。
// 内部面は両側セルが各自で再計算（face値補間は安価）。raw sum を書き、正規化は calcGradient_2_d が行う。
__global__ void calcGradient_cellgather_d
(
 geom_int nCells,
 geom_int* plane_cells,
 geom_int* cell_planes_index, geom_int* cell_planes,
 geom_float* sx, geom_float* sy, geom_float* sz,
 geom_float* fx,
 flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* P, flow_float* T,
 flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
 flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
 flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
 flow_float* drodx, flow_float* drody, flow_float* drodz,
 flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
 flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz,
 flow_float* divU,
 geom_int gradPlaneSkipFrom   // node 弱形式: ip>=この値 (=nNormalPlanes) の境界面はスキップし
                              // calcGradient_b_d (bvar) に委ねる。cell では nPlanes 以上を渡し従来どおり全面処理。
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    flow_float gUxx=0,gUxy=0,gUxz=0, gUyx=0,gUyy=0,gUyz=0, gUzx=0,gUzy=0,gUzz=0;
    flow_float grox=0,groy=0,groz=0, gPx=0,gPy=0,gPz=0, gTx=0,gTy=0,gTz=0, gdiv=0;

    const geom_int st = cell_planes_index[ic];
    const geom_int en = cell_planes_index[ic+1];
    for (geom_int ilp=st; ilp<en; ilp++) {
        geom_int ip  = cell_planes[ilp];
        if (ip >= gradPlaneSkipFrom) continue;   // node: 境界半割面は calcGradient_b_d で処理
        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];
        flow_float sgn = (ic0==ic) ? 1.0 : -1.0;   // 外向き法線符号（owner ic0 が +）
        geom_float f = fx[ip];
        flow_float g = 1.0 - f;

        flow_float Uxf=f*Ux[ic0]+g*Ux[ic1];
        flow_float Uyf=f*Uy[ic0]+g*Uy[ic1];
        flow_float Uzf=f*Uz[ic0]+g*Uz[ic1];
        flow_float Pf =f*P[ic0] +g*P[ic1];
        flow_float Tf =f*T[ic0] +g*T[ic1];
        flow_float rof=f*ro[ic0]+g*ro[ic1];

        flow_float sxx=sgn*sx[ip], syy=sgn*sy[ip], szz=sgn*sz[ip];

        gUxx+=sxx*Uxf; gUxy+=syy*Uxf; gUxz+=szz*Uxf;
        gUyx+=sxx*Uyf; gUyy+=syy*Uyf; gUyz+=szz*Uyf;
        gUzx+=sxx*Uzf; gUzy+=syy*Uzf; gUzz+=szz*Uzf;
        gTx +=sxx*Tf;  gTy +=syy*Tf;  gTz +=szz*Tf;
        gPx +=sxx*Pf;  gPy +=syy*Pf;  gPz +=szz*Pf;
        grox+=sxx*rof; groy+=syy*rof; groz+=szz*rof;
        gdiv+= Uxf*sxx + Uyf*syy + Uzf*szz;
    }

    dUxdx[ic]=gUxx; dUxdy[ic]=gUxy; dUxdz[ic]=gUxz;
    dUydx[ic]=gUyx; dUydy[ic]=gUyy; dUydz[ic]=gUyz;
    dUzdx[ic]=gUzx; dUzdy[ic]=gUzy; dUzdz[ic]=gUzz;
    drodx[ic]=grox; drody[ic]=groy; drodz[ic]=groz;
    dPdx[ic]=gPx;   dPdy[ic]=gPy;   dPdz[ic]=gPz;
    dTdx[ic]=gTx;   dTdy[ic]=gTy;   dTdz[ic]=gTz;
    divU[ic]=gdiv;
}

__global__ void calcGradient_2_d
(
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* roe  ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx  , flow_float* droUxdy , flow_float* droUxdz,
 //flow_float* droUydx  , flow_float* droUydy , flow_float* droUydz,
 //flow_float* droUzdx  , flow_float* droUzdy , flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells ) {

        geom_float volume = vol[ic];

        dUxdx[ic] = dUxdx[ic]/volume;
        dUxdy[ic] = dUxdy[ic]/volume;
        dUxdz[ic] = dUxdz[ic]/volume;
                                    
        dUydx[ic] = dUydx[ic]/volume;
        dUydy[ic] = dUydy[ic]/volume;
        dUydz[ic] = dUydz[ic]/volume;

        dUzdx[ic] = dUzdx[ic]/volume;
        dUzdy[ic] = dUzdy[ic]/volume;
        dUzdz[ic] = dUzdz[ic]/volume;

        dPdx[ic] = dPdx[ic]/volume;
        dPdy[ic] = dPdy[ic]/volume;
        dPdz[ic] = dPdz[ic]/volume;

        drodx[ic] = drodx[ic]/volume;
        drody[ic] = drody[ic]/volume;
        drodz[ic] = drodz[ic]/volume;

        dTdx[ic] = dTdx[ic]/volume;
        dTdy[ic] = dTdy[ic]/volume;
        dTdz[ic] = dTdz[ic]/volume;

        //droUxdx[ic] = droUxdx[ic]/volume;
        //droUxdy[ic] = droUxdy[ic]/volume;
        //droUxdz[ic] = droUxdz[ic]/volume;

        //droUydx[ic] = droUydx[ic]/volume;
        //droUydy[ic] = droUydy[ic]/volume;
        //droUydz[ic] = droUydz[ic]/volume;

        //droUzdx[ic] = droUzdx[ic]/volume;
        //droUzdy[ic] = droUzdy[ic]/volume;
        //droUzdz[ic] = droUzdz[ic]/volume;

        //droedx[ic] = droedx[ic]/volume;
        //droedy[ic] = droedy[ic]/volume;
        //droedz[ic] = droedz[ic]/volume;

        //dHtdx[ic] = dHtdx[ic]/volume;
        //dHtdy[ic] = dHtdy[ic]/volume;
        //dHtdz[ic] = dHtdz[ic]/volume;

        divU[ic] = divU[ic]/volume;
    }
}


// bndFirstOrder 診断/ロバスト化: 境界隣接 CV の再構成勾配 (ro,Ux,Uy,Uz,P) を 0 にし、
// その CV からの 2 次 MUSCL 外挿を 1 次に落とす (壁近傍の過大再構成切り分け/抑制)。
// divU・dT は粘性で使うため残す (本オプションは対流再構成のみを対象とする)。
__global__ void zeroBndNodeGradient_d
(
 geom_int nCells, geom_int* bnode_flag,
 flow_float* drodx, flow_float* drody, flow_float* drodz,
 flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
 flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
 flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
 flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells || bnode_flag[ic] == 0) return;
    drodx[ic]=0; drody[ic]=0; drodz[ic]=0;
    dUxdx[ic]=0; dUxdy[ic]=0; dUxdz[ic]=0;
    dUydx[ic]=0; dUydy[ic]=0; dUydz[ic]=0;
    dUzdx[ic]=0; dUzdy[ic]=0; dUzdz[ic]=0;
    dPdx[ic]=0;  dPdy[ic]=0;  dPdz[ic]=0;
}

// ===================== 最小二乗 (LSQ) 勾配 (node-centered 用) =====================
// median-dual の近壁で Green-Gauss 面勾配が checkerboard を持ち粘性 BL を振動させるため、
// LSQ 勾配 (近傍 CV 中心の重み付き最小二乗) で checkerboard-free に求める。
// b 配列は既存の勾配配列を流用 (accum で b を貯め、solve で in-place に勾配へ上書き)。
// M (対称 3x3) は scratch。境界は bvar 境界値 (壁=速度0) を半割面重心 pc の LSQ 点として含め勾配を閉じる。

// 対称 3x3 を解く。2D (mzz≈0) は xy の 2x2 を解き gz=0。
__device__ __forceinline__ void lsq_solve_sym3(
    double mxx,double mxy,double mxz,double myy,double myz,double mzz,
    double bx,double by,double bz,
    flow_float& gx, flow_float& gy, flow_float& gz)
{
    const double tol = 1.0e-12 * (mxx + myy + 1.0e-300);
    if (mzz <= tol) { // 2D (z 方向に格子が無い)
        double det = mxx*myy - mxy*mxy;
        double idet = (fabs(det) > 1.0e-300) ? 1.0/det : 0.0;
        gx = (flow_float)((myy*bx - mxy*by)*idet);
        gy = (flow_float)((mxx*by - mxy*bx)*idet);
        gz = (flow_float)0.0;
    } else {
        double c00=myy*mzz-myz*myz, c01=mxz*myz-mxy*mzz, c02=mxy*myz-mxz*myy;
        double det=mxx*c00+mxy*c01+mxz*c02;
        double idet=(fabs(det)>1.0e-300)?1.0/det:0.0;
        double c11=mxx*mzz-mxz*mxz, c12=mxz*mxy-mxx*myz, c22=mxx*myy-mxy*mxy;
        gx=(flow_float)((c00*bx+c01*by+c02*bz)*idet);
        gy=(flow_float)((c01*bx+c11*by+c12*bz)*idet);
        gz=(flow_float)((c02*bx+c12*by+c22*bz)*idet);
    }
}

// 内部面 (ip<nNormalPlanes) について per-cell gather で M と b を貯める。
__global__ void lsqGrad_accumInternal_d(
    geom_int nCells, geom_int* plane_cells,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int nNormalPlanes,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* P, flow_float* T,
    flow_float* Mxx, flow_float* Mxy, flow_float* Mxz, flow_float* Myy, flow_float* Myz, flow_float* Mzz,
    flow_float* drodx, flow_float* drody, flow_float* drodz,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
    flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    double mxx=0,mxy=0,mxz=0,myy=0,myz=0,mzz=0;
    double bro[3]={0,0,0}, bux[3]={0,0,0}, buy[3]={0,0,0}, buz[3]={0,0,0}, bp[3]={0,0,0}, bt[3]={0,0,0};
    const double cx=ccx[ic], cy=ccy[ic], cz=ccz[ic];
    const double r0=ro[ic], ux0=Ux[ic], uy0=Uy[ic], uz0=Uz[ic], p0=P[ic], t0=T[ic];
    const geom_int st=cell_planes_index[ic], en=cell_planes_index[ic+1];
    for (geom_int ilp=st; ilp<en; ++ilp) {
        const geom_int ip=cell_planes[ilp];
        if (ip>=nNormalPlanes) continue;          // 境界面は別途 bvar で
        const geom_int ic0=plane_cells[2*ip+0], ic1=plane_cells[2*ip+1];
        const geom_int jc=(ic0==ic)?ic1:ic0;
        const double dx=ccx[jc]-cx, dy=ccy[jc]-cy, dz=ccz[jc]-cz;
        const double w=1.0/fmax(dx*dx+dy*dy+dz*dz,1.0e-300);
        mxx+=w*dx*dx; mxy+=w*dx*dy; mxz+=w*dx*dz; myy+=w*dy*dy; myz+=w*dy*dz; mzz+=w*dz*dz;
        const double dr=ro[jc]-r0, dux=Ux[jc]-ux0, duy=Uy[jc]-uy0, duz=Uz[jc]-uz0, dp=P[jc]-p0, dt=T[jc]-t0;
        bro[0]+=w*dx*dr; bro[1]+=w*dy*dr; bro[2]+=w*dz*dr;
        bux[0]+=w*dx*dux; bux[1]+=w*dy*dux; bux[2]+=w*dz*dux;
        buy[0]+=w*dx*duy; buy[1]+=w*dy*duy; buy[2]+=w*dz*duy;
        buz[0]+=w*dx*duz; buz[1]+=w*dy*duz; buz[2]+=w*dz*duz;
        bp[0]+=w*dx*dp;  bp[1]+=w*dy*dp;  bp[2]+=w*dz*dp;
        bt[0]+=w*dx*dt;  bt[1]+=w*dy*dt;  bt[2]+=w*dz*dt;
    }
    Mxx[ic]=mxx; Mxy[ic]=mxy; Mxz[ic]=mxz; Myy[ic]=myy; Myz[ic]=myz; Mzz[ic]=mzz;
    drodx[ic]=bro[0]; drody[ic]=bro[1]; drodz[ic]=bro[2];
    dUxdx[ic]=bux[0]; dUxdy[ic]=bux[1]; dUxdz[ic]=bux[2];
    dUydx[ic]=buy[0]; dUydy[ic]=buy[1]; dUydz[ic]=buy[2];
    dUzdx[ic]=buz[0]; dUzdy[ic]=buz[1]; dUzdz[ic]=buz[2];
    dPdx[ic]=bp[0];   dPdy[ic]=bp[1];   dPdz[ic]=bp[2];
    dTdx[ic]=bt[0];   dTdy[ic]=bt[1];   dTdz[ic]=bt[2];
}

// 境界半割面の bvar 境界値を LSQ 点 (半割面重心 pc) として M,b へ atomicAdd で加える。
__global__ void lsqGrad_accumBoundary_d(
    geom_int nb, geom_int* bplane_plane, geom_int* bplane_cell,
    geom_float* pcx, geom_float* pcy, geom_float* pcz,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    flow_float* rob, flow_float* Uxb, flow_float* Uyb, flow_float* Uzb, flow_float* Psb, flow_float* Tsb,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* P, flow_float* T,
    flow_float* Mxx, flow_float* Mxy, flow_float* Mxz, flow_float* Myy, flow_float* Myz, flow_float* Mzz,
    flow_float* drodx, flow_float* drody, flow_float* drodz,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
    flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    const geom_int ip=bplane_plane[ib], ic=bplane_cell[ib];
    const double dx=pcx[ip]-ccx[ic], dy=pcy[ip]-ccy[ic], dz=pcz[ip]-ccz[ic];
    const double w=1.0/fmax(dx*dx+dy*dy+dz*dz,1.0e-300);
    atomicAdd(&Mxx[ic],(flow_float)(w*dx*dx)); atomicAdd(&Mxy[ic],(flow_float)(w*dx*dy));
    atomicAdd(&Mxz[ic],(flow_float)(w*dx*dz)); atomicAdd(&Myy[ic],(flow_float)(w*dy*dy));
    atomicAdd(&Myz[ic],(flow_float)(w*dy*dz)); atomicAdd(&Mzz[ic],(flow_float)(w*dz*dz));
    const double dr=rob[ib]-ro[ic], dux=Uxb[ib]-Ux[ic], duy=Uyb[ib]-Uy[ic],
                 duz=Uzb[ib]-Uz[ic], dp=Psb[ib]-P[ic], dt=Tsb[ib]-T[ic];
    atomicAdd(&drodx[ic],(flow_float)(w*dx*dr)); atomicAdd(&drody[ic],(flow_float)(w*dy*dr)); atomicAdd(&drodz[ic],(flow_float)(w*dz*dr));
    atomicAdd(&dUxdx[ic],(flow_float)(w*dx*dux)); atomicAdd(&dUxdy[ic],(flow_float)(w*dy*dux)); atomicAdd(&dUxdz[ic],(flow_float)(w*dz*dux));
    atomicAdd(&dUydx[ic],(flow_float)(w*dx*duy)); atomicAdd(&dUydy[ic],(flow_float)(w*dy*duy)); atomicAdd(&dUydz[ic],(flow_float)(w*dz*duy));
    atomicAdd(&dUzdx[ic],(flow_float)(w*dx*duz)); atomicAdd(&dUzdy[ic],(flow_float)(w*dy*duz)); atomicAdd(&dUzdz[ic],(flow_float)(w*dz*duz));
    atomicAdd(&dPdx[ic],(flow_float)(w*dx*dp)); atomicAdd(&dPdy[ic],(flow_float)(w*dy*dp)); atomicAdd(&dPdz[ic],(flow_float)(w*dz*dp));
    atomicAdd(&dTdx[ic],(flow_float)(w*dx*dt)); atomicAdd(&dTdy[ic],(flow_float)(w*dy*dt)); atomicAdd(&dTdz[ic],(flow_float)(w*dz*dt));
}

// M g = b を変数ごとに解いて勾配を in-place に上書きし、divU を更新する。
__global__ void lsqGrad_solve_d(
    geom_int nCells,
    flow_float* Mxx, flow_float* Mxy, flow_float* Mxz, flow_float* Myy, flow_float* Myz, flow_float* Mzz,
    flow_float* drodx, flow_float* drody, flow_float* drodz,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
    flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz,
    flow_float* divU)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double mxx=Mxx[ic], mxy=Mxy[ic], mxz=Mxz[ic], myy=Myy[ic], myz=Myz[ic], mzz=Mzz[ic];
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, drodx[ic],drody[ic],drodz[ic], drodx[ic],drody[ic],drodz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dUxdx[ic],dUxdy[ic],dUxdz[ic], dUxdx[ic],dUxdy[ic],dUxdz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dUydx[ic],dUydy[ic],dUydz[ic], dUydx[ic],dUydy[ic],dUydz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dUzdx[ic],dUzdy[ic],dUzdz[ic], dUzdx[ic],dUzdy[ic],dUzdz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dPdx[ic],dPdy[ic],dPdz[ic], dPdx[ic],dPdy[ic],dPdz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dTdx[ic],dTdy[ic],dTdz[ic], dTdx[ic],dTdy[ic],dTdz[ic]);
    divU[ic] = dUxdx[ic] + dUydy[ic] + dUzdz[ic];
}

void calcGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    flow_float* grad_volume = (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"];
    flow_float* grad_sx = (cfg.isAxisymmetric == 1) ? var.p_d["sx_planar"] : var.p_d["sx"];
    flow_float* grad_sy = (cfg.isAxisymmetric == 1) ? var.p_d["sy_planar"] : var.p_d["sy"];
    flow_float* grad_sz = (cfg.isAxisymmetric == 1) ? var.p_d["sz_planar"] : var.p_d["sz"];
    flow_float* grad_ss = (cfg.isAxisymmetric == 1) ? var.p_d["ss_planar"] : var.p_d["ss"];

    // initialize
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drody"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["divU"] , 0.0, msh.nCells*sizeof(flow_float)));

    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedz"], 0.0, msh.nCells*sizeof(flow_float)));


    // ---- LSQ 勾配 (node-centered, gradLSQ=1): Green-Gauss の checkerboard を避ける ----
    if (cfg.gradLSQ && cfg.discretization == "node") {
        static flow_float *Mxx=nullptr,*Mxy=nullptr,*Mxz=nullptr,*Myy=nullptr,*Myz=nullptr,*Mzz=nullptr;
        static geom_int M_n=0;
        if (Mxx==nullptr || M_n!=msh.nCells) {
            if (Mxx){cudaFree(Mxx);cudaFree(Mxy);cudaFree(Mxz);cudaFree(Myy);cudaFree(Myz);cudaFree(Mzz);}
            size_t sz=sizeof(flow_float)*msh.nCells;
            gpuErrchk(cudaMalloc(&Mxx,sz)); gpuErrchk(cudaMalloc(&Mxy,sz)); gpuErrchk(cudaMalloc(&Mxz,sz));
            gpuErrchk(cudaMalloc(&Myy,sz)); gpuErrchk(cudaMalloc(&Myz,sz)); gpuErrchk(cudaMalloc(&Mzz,sz));
            M_n=msh.nCells;
        }
        lsqGrad_accumInternal_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.nNormalPlanes,
            var.c_d["ccx"],var.c_d["ccy"],var.c_d["ccz"],
            var.c_d["ro"],var.c_d["Ux"],var.c_d["Uy"],var.c_d["Uz"],var.c_d["P"],var.c_d["T"],
            Mxx,Mxy,Mxz,Myy,Myz,Mzz,
            var.c_d["drodx"],var.c_d["drody"],var.c_d["drodz"],
            var.c_d["dUxdx"],var.c_d["dUxdy"],var.c_d["dUxdz"],
            var.c_d["dUydx"],var.c_d["dUydy"],var.c_d["dUydz"],
            var.c_d["dUzdx"],var.c_d["dUzdy"],var.c_d["dUzdz"],
            var.c_d["dPdx"],var.c_d["dPdy"],var.c_d["dPdz"],
            var.c_d["dTdx"],var.c_d["dTdy"],var.c_d["dTdz"]);
        for (auto& bc : msh.bconds) {
            if (bc.iPlanes.size()==0) continue;
            lsqGrad_accumBoundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
                (geom_int)bc.iPlanes.size(), bc.map_bplane_plane_d, bc.map_bplane_cell_d,
                var.p_d["pcx"],var.p_d["pcy"],var.p_d["pcz"],
                var.c_d["ccx"],var.c_d["ccy"],var.c_d["ccz"],
                bc.bvar_d["ro"],bc.bvar_d["Ux"],bc.bvar_d["Uy"],bc.bvar_d["Uz"],bc.bvar_d["Ps"],bc.bvar_d["Ts"],
                var.c_d["ro"],var.c_d["Ux"],var.c_d["Uy"],var.c_d["Uz"],var.c_d["P"],var.c_d["T"],
                Mxx,Mxy,Mxz,Myy,Myz,Mzz,
                var.c_d["drodx"],var.c_d["drody"],var.c_d["drodz"],
                var.c_d["dUxdx"],var.c_d["dUxdy"],var.c_d["dUxdz"],
                var.c_d["dUydx"],var.c_d["dUydy"],var.c_d["dUydz"],
                var.c_d["dUzdx"],var.c_d["dUzdy"],var.c_d["dUzdz"],
                var.c_d["dPdx"],var.c_d["dPdy"],var.c_d["dPdz"],
                var.c_d["dTdx"],var.c_d["dTdy"],var.c_d["dTdz"]);
        }
        lsqGrad_solve_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, Mxx,Mxy,Mxz,Myy,Myz,Mzz,
            var.c_d["drodx"],var.c_d["drody"],var.c_d["drodz"],
            var.c_d["dUxdx"],var.c_d["dUxdy"],var.c_d["dUxdz"],
            var.c_d["dUydx"],var.c_d["dUydy"],var.c_d["dUydz"],
            var.c_d["dUzdx"],var.c_d["dUzdy"],var.c_d["dUzdz"],
            var.c_d["dPdx"],var.c_d["dPdy"],var.c_d["dPdz"],
            var.c_d["dTdx"],var.c_d["dTdy"],var.c_d["dTdz"],
            var.c_d["divU"]);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        return;   // LSQ は正規化不要 (直接勾配を解く)。GG 経路はスキップ。
    }

    // sum over planes
    // K5: atomicAdd 版 calcGradient_1_d を cell-gather 版へ置換（raw sum を書く。正規化は下の _2_d）。
    calcGradient_cellgather_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
        msh.nCells,
        msh.map_plane_cells_d,
        msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
        grad_sx , grad_sy , grad_sz ,
        var.p_d["fx"],
        var.c_d["ro"] , var.c_d["Ux"] , var.c_d["Uy"] , var.c_d["Uz"] , var.c_d["P"] , var.c_d["T"] ,
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
        var.c_d["divU"],
        // node 弱形式: 境界半割面 (ip>=nNormalPlanes) を cellgather でスキップし calcGradient_b_d (bvar+真の
        // 半割面) に委ねる。壁ノードは壁上に乗りミラーゴーストが退化するため、ミラー経由の cellgather 勾配が
        // 不正確 → 粘性 BL 振動の原因。弱形式の境界寄与で正しく閉じる。cell は nPlanes (全面処理, 従来)。
        (cfg.discretization == "node") ? msh.nNormalPlanes : msh.nPlanes
    ) ;

    // node 弱形式 境界勾配: 各境界半割面の φ_b·S_b を bvar から加える (cellgather は内部のみ計算済み)。
    // 軸対称は cellgather と同じ planar 面積 (grad_sx=sx_planar) を使い r 重みを一致させる。
    if (cfg.discretization == "node") {
        for (auto& bc : msh.bconds)
        {
            // periodic は DOF 同一視 (§4.5) で扱うため境界半割面の勾配寄与を加えない。
            // 勾配は内部双対面のみ (片側) で計算し、後段 periodicGradientGather で両側を厳密合併する。
            if (bc.bcondKind == "periodic") continue;
            calcGradient_b_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
                bc.iPlanes.size(),
                bc.map_bplane_plane_d, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                grad_sx          , grad_sy       , grad_sz       , grad_ss,
                bc.bvar_d["ro"], bc.bvar_d["roUx"], bc.bvar_d["roUy"], bc.bvar_d["roUz"], bc.bvar_d["roe"],
                bc.bvar_d["Ux"], bc.bvar_d["Uy"], bc.bvar_d["Uz"],
                bc.bvar_d["Tt"], bc.bvar_d["Pt"], bc.bvar_d["Ts"], bc.bvar_d["Ps"], bc.bvar_d["Ht"],
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
                var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
                var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
                var.c_d["divU"]
            ) ;
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
    }


    calcGradient_2_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        grad_volume, var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        grad_sx , grad_sy , grad_sz , grad_ss,

        // basic variables
        var.c_d["ro"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"]  , 
        var.c_d["T"]  ,
        var.c_d["roe"] ,
        var.c_d["Ht"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

        //var.c_d["droUxdx"] , var.c_d["droUxdy"] , var.c_d["droUxdz"],
        //var.c_d["droUydx"] , var.c_d["droUydy"] , var.c_d["droUydz"],
        //var.c_d["droUzdx"] , var.c_d["droUzdy"] , var.c_d["droUzdz"],
        //var.c_d["droedx"]  , var.c_d["droedy"]  , var.c_d["droedz"],
        //var.c_d["dHtdx"]  , var.c_d["dHtdy"]  , var.c_d["dHtdz"],
 
        var.c_d["divU"]
    ) ;

    // 境界隣接 CV の再構成勾配を 0 に (境界側 1 次化)。bndFirstOrder=1 のときのみ。
    if (cfg.bndFirstOrder == 1 && msh.bnode_flag_d != nullptr) {
        zeroBndNodeGradient_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, msh.bnode_flag_d,
            var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"]
        );
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}