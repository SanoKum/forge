#include "convectiveFlux_d.cuh"
#include <cstdlib>
#include <cstdio>
#include <vector>

__global__ void viscousFlux_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 flow_float mu ,  flow_float Prt , flow_float* cp , flow_float* thermCond,
 flow_float* vis_lam   , flow_float* vis_turb  ,

 // variables
//flow_float* convx , flow_float* convy , flow_float* convz,
// flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* ro   ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Ps  ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* Ts  ,

 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 // 軸対称: 完全な発散 ∇·u = ∂xUx+∂yUy+u_r/r (axisym_divU)。-2/3 μ(∇·u) の体積粘性項を
 // τθθ ソース (axisymmetricSource_d) と整合させるため planar 面でもフープ項込みの divu を使う。
 int isAxisymmetric, flow_float* axisym_divU
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    //if (ip < nPlanes) { 
    if (ip < nNormalPlanes) {

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float f = fx[ip];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float ccx_0 = ccx[ic0];
        flow_float ccy_0 = ccy[ic0];
        flow_float ccz_0 = ccz[ic0];

        flow_float ccx_1 = ccx[ic1];
        flow_float ccy_1 = ccy[ic1];
        flow_float ccz_1 = ccz[ic1];

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float Uxf = f*Ux[ic0] + (1.0-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0-f)*Uz[ic1];

        flow_float dUxdxf = f*dUxdx[ic0] + (1.0-f)*dUxdx[ic1];
        flow_float dUxdyf = f*dUxdy[ic0] + (1.0-f)*dUxdy[ic1];
        flow_float dUxdzf = f*dUxdz[ic0] + (1.0-f)*dUxdz[ic1];

        flow_float dUydxf = f*dUydx[ic0] + (1.0-f)*dUydx[ic1];
        flow_float dUydyf = f*dUydy[ic0] + (1.0-f)*dUydy[ic1];
        flow_float dUydzf = f*dUydz[ic0] + (1.0-f)*dUydz[ic1];

        flow_float dUzdxf = f*dUzdx[ic0] + (1.0-f)*dUzdx[ic1];
        flow_float dUzdyf = f*dUzdy[ic0] + (1.0-f)*dUzdy[ic1];
        flow_float dUzdzf = f*dUzdz[ic0] + (1.0-f)*dUzdz[ic1];

        flow_float dTdxf = f*dTdx[ic0] + (1.0-f)*dTdx[ic1];
        flow_float dTdyf = f*dTdy[ic0] + (1.0-f)*dTdy[ic1];
        flow_float dTdzf = f*dTdz[ic0] + (1.0-f)*dTdz[ic1];
#ifdef VISC_DIAG_NOGRAD
        // 診断: 面勾配 (Green-Gauss) 由来の項 (転置・非直交補正) を 0 にし、over-relaxed 法線 Laplacian のみ残す。
        dUxdxf=dUxdyf=dUxdzf=0; dUydxf=dUydyf=dUydzf=0; dUzdxf=dUzdyf=dUzdzf=0; dTdxf=dTdyf=dTdzf=0;
#endif

        // implicit 側（timeIntegration_d.cu）と同じ fabs()+下限で分母を保護する
        flow_float D_safe  = max(fabs(dcc_x*sxx +dcc_y*syy +dcc_z*szz), (flow_float)1.0e-30);
        flow_float delta   = dcc*sss*sss/D_safe; // over-relaxed
        flow_float delta_x = dcc_x*sss*sss/D_safe;
        flow_float delta_y = dcc_y*sss*sss/D_safe;
        flow_float delta_z = dcc_z*sss*sss/D_safe;
        flow_float k_x = sxx - delta_x;
        flow_float k_y = syy - delta_y;
        flow_float k_z = szz - delta_z;
        // 軸対称は完全発散 (u_r/r 込み) を面補間。デカルトは従来どおり ∂xUx+∂yUy+∂zUz。
        flow_float divu = (isAxisymmetric == 1)
            ? (f*axisym_divU[ic0] + (1.0-f)*axisym_divU[ic1])
            : (dUxdxf+dUydyf+dUzdzf);

        flow_float v_lam  = f*vis_lam [ic0] + (1.0-f)*vis_lam [ic1] ;
        flow_float v_turb = f*vis_turb[ic0] + (1.0-f)*vis_turb[ic1] ;
        flow_float mu_total = v_lam + v_turb;

        // 完全な Newton 応力 tau_ij S_j = mu(du_i/dx_j + du_j/dx_i)S_j - (2/3)mu divu S_i。
        // 第1項 (Laplacian, mu grad(u_i).S) は over-relaxed: 法線スカラー delta + 同成分勾配.k。
        // 第2項 (転置, mu du_j/dx_i S_j) は面平均勾配にフル S を内積。第3項 (発散) は成分 s**。
        flow_float tau_x = mu_total*((Ux[ic1] -Ux[ic0])/dcc)*delta;
        tau_x += mu_total*(dUxdxf*k_x +dUxdyf*k_y +dUxdzf*k_z);
        tau_x += mu_total*(dUxdxf*sxx +dUydxf*syy +dUzdxf*szz);
        tau_x += -mu_total*2.0/3.0*(divu)*sxx;

        flow_float tau_y = mu_total*((Uy[ic1] -Uy[ic0])/dcc)*delta;
        tau_y += mu_total*(dUydxf*k_x +dUydyf*k_y +dUydzf*k_z);
        tau_y += mu_total*(dUxdyf*sxx +dUydyf*syy +dUzdyf*szz);
        tau_y += -mu_total*2.0/3.0*(divu)*syy;

        flow_float tau_z = mu_total*((Uz[ic1] -Uz[ic0])/dcc)*delta;
        tau_z += mu_total*(dUzdxf*k_x +dUzdyf*k_y +dUzdzf*k_z);
        tau_z += mu_total*(dUxdzf*sxx +dUydzf*syy +dUzdzf*szz);
        tau_z += -mu_total*2.0/3.0*(divu)*szz;

        // 乱流熱伝導: 有効熱伝導率 k_eff = k_lam + cp*mu_turb/Pr_t。
        // 応力(摩擦発熱)が mu_total=v_lam+v_turb を使うのに対し従来は層流 k のみで
        // 熱を運んでいたため、乱流境界層で散逸熱が逃げ場を失い T が全温を超えて
        // overshoot していた。RANS のエネルギー保存には乱流熱伝導が必須。
        // cp は thermally-perfect の温度依存を反映するためセル配列を面平均で使う。
        flow_float cp_face = f*cp[ic0] + (1.0-f)*cp[ic1];
        flow_float tc_face = f*thermCond[ic0] + (1.0-f)*thermCond[ic1];
        tc_face += cp_face*v_turb/Prt;
        flow_float heatflux = tc_face*((Ts[ic1] -Ts[ic0])/dcc)*delta;
        heatflux += tc_face*(dTdxf*k_x +dTdyf*k_y +dTdzf*k_z);

        flow_float res_ro_temp   = 0.0;
        flow_float res_roUx_temp = tau_x;
        flow_float res_roUy_temp = tau_y;
        flow_float res_roUz_temp = tau_z;
        flow_float res_roe_temp  = tau_x*Uxf +tau_y*Uyf +tau_z*Uzf; 
        res_roe_temp += heatflux;

        atomicAdd(&res_ro[ic0]  , res_ro_temp);
        atomicAdd(&res_roUx[ic0], res_roUx_temp);
        atomicAdd(&res_roUy[ic0], res_roUy_temp);
        atomicAdd(&res_roUz[ic0], res_roUz_temp);
        atomicAdd(&res_roe[ic0] , res_roe_temp);

        atomicAdd(&res_ro[ic1]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic1], -res_roUx_temp);
        atomicAdd(&res_roUy[ic1], -res_roUy_temp);
        atomicAdd(&res_roUz[ic1], -res_roUz_temp);
        atomicAdd(&res_roe[ic1] , -res_roe_temp);
    }

    __syncthreads();
}


__global__ void viscousFlux_wall_d
( 
  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 flow_float mu ,  flow_float Prt , flow_float* cp , flow_float* thermCond,
 flow_float* vis_lam   , flow_float* vis_turb  ,

 // variables
//flow_float* convx , flow_float* convy , flow_float* convz,
// flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* ro   ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Ps  ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* Ts  ,

 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 // bvar
 flow_float* ypls_b ,
 flow_float* twall_x_b ,
 flow_float* twall_y_b ,
 flow_float* twall_z_b ,
 flow_float* Ux_b ,
 flow_float* Uy_b ,
 flow_float* Uz_b ,
 flow_float* Ts_b ,
//flow_float* sx_b ,
 //flow_float* sy_b ,
 //flow_float* sz_b

 // 軸対称: 完全発散 (u_r/r 込み) を体積粘性項に使う (内部面と同趣旨)
 int isAxisymmetric, flow_float* axisym_divU,

 // SST automatic wall treatment (docs/turbulence §6.5 (c)):
 //   wallTreatment==1 で接線せん断を modeled τ_w=ρu_τ² に置換。utau_b は事前算出済み u_τ。
 int wallTreatment, flow_float* utau_b,
 // node-centered (1): 壁法線 Laplacian/熱流束を ghost+dcc でなくセル勾配 ∇φ·S (bvar 壁閉包) で評価。
 // 壁ノードが壁面に乗り dcc≈0 に退化するのを回避する。cell (0) は従来の (φ[ig]-φ[ic])/dcc。
 int isNode
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

        flow_float ccx_0 = ccx[ic];
        flow_float ccy_0 = ccy[ic];
        flow_float ccz_0 = ccz[ic];

        flow_float pcx_1 = pcx[ip];
        flow_float pcy_1 = pcy[ip];
        flow_float pcz_1 = pcz[ip];

        flow_float ccx_1 = ccx[ig];
        flow_float ccy_1 = ccy[ig];
        flow_float ccz_1 = ccz[ig];

        //flow_float dcc_x = pcx_1 - ccx_0;
        //flow_float dcc_y = pcy_1 - ccy_0;
        //flow_float dcc_z = pcz_1 - ccz_0;

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;
 
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float dUxdxf = dUxdx[ic] ;
        flow_float dUxdyf = dUxdy[ic] ;
        flow_float dUxdzf = dUxdz[ic] ;

        flow_float dUydxf = dUydx[ic] ;
        flow_float dUydyf = dUydy[ic] ;
        flow_float dUydzf = dUydz[ic] ;

        flow_float dUzdxf = dUzdx[ic] ;
        flow_float dUzdyf = dUzdy[ic] ;
        flow_float dUzdzf = dUzdz[ic] ;

        flow_float v_turb = vis_turb[ic] ;
        flow_float mu_total = vis_lam[ic] + v_turb;

        // 軸対称は完全発散 (u_r/r 込み, セル値)。壁面はセル中心勾配を直接使うため cell の axisym_divU。
        flow_float divu = (isAxisymmetric == 1) ? axisym_divU[ic] : (dUxdxf+dUydyf+dUzdzf);

        // ミラーゴースト (u^g=-u^c) では d∥S なので over-relaxed の delta=sss、k=0。
        // 法線項: 法線勾配 (u^g-u^c)/dcc に面積 sss。転置項: セル中心勾配にフル S。発散項: 成分 s**。
        // 注: dcc フロア等で物理 flux をマスクしない方針 (退化 dcc はメッシュ/境界クロージャの欠陥であり、
        // 直すべきはメッシュとゴーストレス弱形式。診断は viscousWallDiag_d で別途行う)。
        // 法線 Laplacian 項: node は ∇U·S (セル勾配, bvar 壁閉包) で退化 dcc を回避、cell は従来の法線差分。
        flow_float tau_x = (isNode != 0) ? mu_total*(dUxdxf*sxx +dUxdyf*syy +dUxdzf*szz)
                                         : mu_total*((Ux[ig] - Ux[ic])/dcc)*sss;
        tau_x += mu_total*(dUxdxf*sxx +dUydxf*syy +dUzdxf*szz);
        tau_x += -mu_total*2.0/3.0*(divu)*sxx;

        flow_float tau_y = (isNode != 0) ? mu_total*(dUydxf*sxx +dUydyf*syy +dUydzf*szz)
                                         : mu_total*((Uy[ig] - Uy[ic])/dcc)*sss;
        tau_y += mu_total*(dUxdyf*sxx +dUydyf*syy +dUzdyf*szz);
        tau_y += -mu_total*2.0/3.0*(divu)*syy;

        flow_float tau_z = (isNode != 0) ? mu_total*(dUzdxf*sxx +dUzdyf*syy +dUzdzf*szz)
                                         : mu_total*((Uz[ig] - Uz[ic])/dcc)*sss;
        tau_z += mu_total*(dUxdzf*sxx +dUydzf*syy +dUzdzf*szz);
        tau_z += -mu_total*2.0/3.0*(divu)*szz;

        // SST automatic wall treatment (docs/turbulence §6.5 (c)): 粗メッシュで分子勾配が τ_w を
        // 過小評価するため、接線せん断を modeled τ_w=ρu_τ² (有効壁粘性) に置換する。u_τ は壁関数で
        // 事前算出済み (ransWallFunction)。法線は単位法線 n、接線方向 ê_t はセル接線速度から取る。
        // y⁺→0 では u_τ²=νU_t/y より分子勾配に縮退、対数層では壁関数化し、全層で連続。
        if (wallTreatment == 1) {
            const flow_float n_x = sxx / sss;
            const flow_float n_y = syy / sss;
            const flow_float n_z = szz / sss;
            const flow_float uc = Ux[ic];
            const flow_float vc = Uy[ic];
            const flow_float wc = Uz[ic];
            const flow_float un = uc*n_x + vc*n_y + wc*n_z;
            const flow_float utx = uc - un*n_x;
            const flow_float uty = vc - un*n_y;
            const flow_float utz = wc - un*n_z;
            const flow_float Ut = sqrt(utx*utx + uty*uty + utz*utz);
            const flow_float utau_w = utau_b[ib];
            const flow_float tauw = ro[ic]*utau_w*utau_w;   // modeled 壁せん断応力の大きさ
            flow_float etx = 0.0, ety = 0.0, etz = 0.0;
            if (Ut > 1.0e-12) { etx = utx/Ut; ety = uty/Ut; etz = utz/Ut; }
            // 流れを減速させる向き (-ê_t) に τ_w·面積 を課す。法線粘性・体積項は落とす (壁関数)。
            tau_x = -tauw*etx*sss;
            tau_y = -tauw*ety*sss;
            tau_z = -tauw*etz*sss;
        }

        // 乱流熱伝導 (内部面と同じ k_eff = k_lam + cp*mu_turb/Pr_t)。断熱壁ではミラーゴーストで
        // Ts[ig]=Ts[ic] のため寄与 0、isothermal 壁では有効。
        flow_float tc_w = thermCond[ic] + cp[ic]*v_turb/Prt;
        // 熱流束法線項: node は ∇T·S (セル勾配) で退化 dcc を回避。cell は従来の (Ts[ig]-Ts[ic])/dcc。
        flow_float heatflux = (isNode != 0) ? tc_w*(dTdx[ic]*sxx +dTdy[ic]*syy +dTdz[ic]*szz)
                                            : tc_w*((Ts[ig]- Ts[ic])/dcc)*sss;

        flow_float res_ro_temp   = 0.0;
        flow_float res_roUx_temp = tau_x;
        flow_float res_roUy_temp = tau_y;
        flow_float res_roUz_temp = tau_z;
        flow_float res_roe_temp  = tau_x*Ux_b[ib] +tau_y*Uy_b[ib] +tau_z*Uz_b[ib]; 

        //printf("res_roe_temp=%e\n", res_roe_temp);
        //printf("Uyb[ib]=%e\n", Uy_b[ib]);
        //printf("Uzb[ib]=%e\n", Uz_b[ib]);

        res_roe_temp += heatflux;

        //res_roe_temp = 0.0;

        atomicAdd(&res_ro[ic]  , res_ro_temp);
        atomicAdd(&res_roUx[ic], res_roUx_temp);
        atomicAdd(&res_roUy[ic], res_roUy_temp);
        atomicAdd(&res_roUz[ic], res_roUz_temp);
        atomicAdd(&res_roe[ic] , res_roe_temp);

        twall_x_b[ib] = tau_x/sss;
        twall_y_b[ib] = tau_y/sss;
        twall_z_b[ib] = tau_z/sss;

        flow_float twall = sqrt(tau_x*tau_x + tau_y*tau_y + tau_z*tau_z)/sss;
        flow_float utau = sqrt(twall/ro[ic]);

        // mode 1 では壁関数カーネル (ransWallFunction) が y⁺=u_τ y/ν_lam を既に格納済み。
        // ここで dcc/mu_total ベースの値で上書きすると定義が不整合になるため mode 0 のみ更新する。
        if (wallTreatment != 1) {
            ypls_b[ib] = ro[ic]*utau*dcc/mu_total;
        }
        //if (ib == 1) {
        //    printf("ib = %d\n", ib);
        //    printf("ip = %d\n", ip);
        //    printf("tau_x = %e, tau_y = %e, tau_z = %e\n", tau_x, tau_y, tau_z);
        //    printf("c Ux = %e, Uy = %e, Uz = %e\n", Ux[ic], Uy[ic], Uz[ic]);
        //    printf("g Ux = %e, Uy = %e, Uz = %e\n", Ux[ig], Uy[ig], Uz[ig]);
        //    printf("b Ux = %e, Uy = %e, Uz = %e\n", Ux_b[ib], Uy_b[ib], Uz_b[ib]);
        //}
    }

    __syncthreads();
}

// ---- 壁粘性 flux の距離診断 (env FORGE_VISC_WALL_DIAG=1 で 1 回だけ実行) ----
// 壁半割面ごとに、signed wall distance dn=(pc-cc)·n、ミラーゴースト距離 dcc=|cc_ghost-cc|、
// 接線オフセット |(pc-cc)-dn n| を書き出す。node-centered で dcc≈0 (壁ノードが壁面上) に
// 退化しているか、cell-centered で dcc≈2dn が成立しているかを切り分けるための診断。
__global__ void viscousWallDiag_d
(
 geom_int nb,
 geom_int* bplane_plane, geom_int* bplane_cell, geom_int* bplane_cell_ghst,
 geom_float* ccx, geom_float* ccy, geom_float* ccz,
 geom_float* pcx, geom_float* pcy, geom_float* pcz,
 geom_float* sx , geom_float* sy , geom_float* sz, geom_float* ss,
 flow_float* dn_o, flow_float* dcc_o, flow_float* tang_o
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;
    if (ib < nb) {
        geom_int ip = bplane_plane[ib];
        geom_int ic = bplane_cell[ib];
        geom_int ig = bplane_cell_ghst[ib];
        flow_float nx = sx[ip]/ss[ip], ny = sy[ip]/ss[ip], nz = sz[ip]/ss[ip];
        flow_float dx = pcx[ip]-ccx[ic], dy = pcy[ip]-ccy[ic], dz = pcz[ip]-ccz[ic];
        flow_float dn = dx*nx + dy*ny + dz*nz;
        flow_float tx = dx - dn*nx, ty = dy - dn*ny, tz = dz - dn*nz;
        flow_float gx = ccx[ig]-ccx[ic], gy = ccy[ig]-ccy[ic], gz = ccz[ig]-ccz[ic];
        dn_o[ib]   = dn;
        dcc_o[ib]  = sqrt(gx*gx + gy*gy + gz*gz);
        tang_o[ib] = sqrt(tx*tx + ty*ty + tz*tz);
    }
}

void viscousFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // 距離診断 (1 回限り)。FORGE_VISC_WALL_DIAG=1 のとき壁半割面の dn/dcc/tangential を集計表示。
    static bool diag_done = false;
    if (!diag_done && getenv("FORGE_VISC_WALL_DIAG") && atoi(getenv("FORGE_VISC_WALL_DIAG")) != 0) {
        diag_done = true;
        for (auto& bc : msh.bconds) {
            if (!(bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")) continue;
            geom_int nbf = (geom_int)bc.iPlanes.size();
            if (nbf == 0) { printf("[viscWallDiag] bcond '%s' (physID %d): node 壁 plane なし\n",
                                    bc.bcondKind.c_str(), bc.physID); continue; }
            flow_float *dn_d, *dcc_d, *tang_d;
            gpuErrchk(cudaMalloc(&dn_d , sizeof(flow_float)*nbf));
            gpuErrchk(cudaMalloc(&dcc_d, sizeof(flow_float)*nbf));
            gpuErrchk(cudaMalloc(&tang_d,sizeof(flow_float)*nbf));
            viscousWallDiag_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nbf, bc.map_bplane_plane_d, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"], var.p_d["pcy"], var.p_d["pcz"],
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
                dn_d, dcc_d, tang_d);
            gpuErrchkKernelSync();
            std::vector<flow_float> dn(nbf), dcc(nbf), tang(nbf);
            gpuErrchk(cudaMemcpy(dn.data() , dn_d , sizeof(flow_float)*nbf, cudaMemcpyDeviceToHost));
            gpuErrchk(cudaMemcpy(dcc.data(), dcc_d, sizeof(flow_float)*nbf, cudaMemcpyDeviceToHost));
            gpuErrchk(cudaMemcpy(tang.data(),tang_d,sizeof(flow_float)*nbf, cudaMemcpyDeviceToHost));
            // 統計
            double dn_mn=1e30,dn_mx=-1e30,dn_amn=1e30, dcc_mn=1e30,dcc_mx=-1e30;
            double tang_mx=-1e30, ratio_mn=1e30, ratio_mx=-1e30; int nFlip=0;
            for (geom_int i=0;i<nbf;i++){
                double adn=fabs(dn[i]);
                dn_mn=fmin(dn_mn,dn[i]); dn_mx=fmax(dn_mx,dn[i]); dn_amn=fmin(dn_amn,adn);
                dcc_mn=fmin(dcc_mn,dcc[i]); dcc_mx=fmax(dcc_mx,dcc[i]);
                tang_mx=fmax(tang_mx,tang[i]);
                double r = dcc[i]/fmax(2.0*adn,1e-30);
                ratio_mn=fmin(ratio_mn,r); ratio_mx=fmax(ratio_mx,r);
                if (dn[i]<0) nFlip++;
            }
            printf("[viscWallDiag] bcond '%s' physID %d nFaces=%ld disc=%s\n",
                   bc.bcondKind.c_str(), bc.physID, (long)nbf, cfg.discretization.c_str());
            printf("  dn (signed): min=%.3e max=%.3e | |dn| min=%.3e | dn<0 count=%d\n",
                   dn_mn, dn_mx, dn_amn, nFlip);
            printf("  dcc=|cc_ghost-cc|: min=%.3e max=%.3e\n", dcc_mn, dcc_mx);
            printf("  dcc/(2|dn|): min=%.3e max=%.3e (cell では ~1, node 退化なら ~0/大ばらつき)\n",
                   ratio_mn, ratio_mx);
            printf("  tangential offset |(pc-cc)-dn n|: max=%.3e (node 壁ノードでは支配的のはず)\n", tang_mx);
            cudaFree(dn_d); cudaFree(dcc_d); cudaFree(tang_d);
        }
    }

    // ------------------------------
    // *** sum over normal planes ***
    // ------------------------------
    viscousFlux_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        cfg.visc , cfg.turbulentPrandtl , var.c_d["cp"] , var.c_d["thermCond"],
        var.c_d["vis_lam"], var.c_d["vis_turb"],

        // basic variables
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

        var.c_d["res_ro"] ,
        var.c_d["res_roUx"] ,
        var.c_d["res_roUy"] ,
        var.c_d["res_roUz"] ,
        var.c_d["res_roe"]  ,
       
        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["dTdx"] , var.c_d["dTdy"] , var.c_d["dTdz"],
        cfg.isAxisymmetric, var.c_d["axisym_divU"]
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "wall" or bc.bcondKind == "wall_isothermal") {
            viscousFlux_wall_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
                // mesh structure
                bc.iPlanes.size(),
                bc.map_bplane_plane_d,  
                bc.map_bplane_cell_d,  
                bc.map_bplane_cell_ghst_d,

                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

                cfg.visc , cfg.turbulentPrandtl , var.c_d["cp"] , var.c_d["thermCond"],
                var.c_d["vis_lam"], var.c_d["vis_turb"],

                // basic variables
                //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
                //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
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

                var.c_d["res_ro"] ,
                var.c_d["res_roUx"] ,
                var.c_d["res_roUy"] ,
                var.c_d["res_roUz"] ,
                var.c_d["res_roe"]  ,
       
                // gradient
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

                bc.bvar_d["ypls"],
                bc.bvar_d["twall_x"],
                bc.bvar_d["twall_y"],
                bc.bvar_d["twall_z"],
                bc.bvar_d["Ux"],
                bc.bvar_d["Uy"],
                bc.bvar_d["Uz"],
                bc.bvar_d["Ts"],

                //bc.bvar_d["sx"],
                //bc.bvar_d["sy"],
                //bc.bvar_d["sz"]
                cfg.isAxisymmetric, var.c_d["axisym_divU"],

                // SST automatic wall treatment: modeled τ_w 用の flag と u_τ (ransWallFunction で算出済み)
                (cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? cfg.wallTreatmentSST : 0,
                bc.bvar_d["utau"],
                (cfg.discretization == "node") ? 1 : 0   // node: 壁法線/熱流束を ∇φ·S で評価 (ghostless)
            ) ;
        }
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}