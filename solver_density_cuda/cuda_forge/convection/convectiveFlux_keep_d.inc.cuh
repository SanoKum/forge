// =============================================================================
// 対流フラックス KEEP_d の実装断片。
//   - KEEP (Kinetic Energy & Entropy Preserving) 中心流束 + Roe 特性行列散逸。
//     legacy/convectiveFlux_ausm_keep_d.inc.cuh の KEEP_FVS_d を modern bundled API
//     (FaceGeom/PrimState/...) へ移植し、無効化 (if(false)) されていた KEEP 中心流束を
//     有効化したもの。中心流束は隣接セル/ノード対 (ic0,ic1) の生値で構成し (KE/エントロピー
//     保存)、散逸は MUSCL 再構成した L/R 状態の Roe 行列で与える (低マッハ/LES 用に後で
//     Ducros 等で散逸を絞る余地あり: plans/active/convection-keep-revive-node.md)。
//   - convectiveFlux_d.cu に textual include される断片であり **スタンドアロン TU ではない**。
//     CMake の source には足さないこと (interp_dispatch / apply_ducros_limiter は include 元上方)。
//   - 周回面は geom.nLoopPlanes (= convPlaneBound)。node 弱形式では内部双対面のみ (境界は
//     convectiveFlux_boundary_d)。cell では内部+境界 ghost を周回し、境界 plane は 1 次化する。
// =============================================================================

__global__ void KEEP_d
(
 int conv_scheme, int limit_scheme,
 flow_float ga,
 CondArgs      cnd,    // (CPG KEEP では未使用。dispatch 統一のため受け取る)
 FaceGeom      geom,
 PrimState     st,
 ResidualOut   reso,
 LimiterFields lim,
 GradFields    grd
)
{
    // --- struct 引数をローカルへ展開 ---
    const geom_int nCells       = geom.nCells;
    geom_int*      plane_cells  = geom.plane_cells;
    const geom_int nNormal_halo_Planes = geom.nLoopPlanes;
    geom_int*      normal_halo_planes_d = geom.loop_planes;
    geom_float *vol=geom.vol, *ccx=geom.ccx, *ccy=geom.ccy, *ccz=geom.ccz;
    geom_float *pcx=geom.pcx, *pcy=geom.pcy, *pcz=geom.pcz, *fx=geom.fx;
    geom_float *sx=geom.sx, *sy=geom.sy, *sz=geom.sz, *ss=geom.ss;
    flow_float* massflux=geom.massflux;

    flow_float *ro=st.ro, *Ux=st.Ux, *Uy=st.Uy, *Uz=st.Uz, *Ps=st.Ps;

    flow_float *res_ro=reso.res_ro, *res_roUx=reso.res_roUx, *res_roUy=reso.res_roUy, *res_roUz=reso.res_roUz, *res_roe=reso.res_roe;

    flow_float *limiter_ro=lim.limiter_ro, *limiter_Ux=lim.limiter_Ux, *limiter_Uy=lim.limiter_Uy, *limiter_Uz=lim.limiter_Uz, *limiter_P=lim.limiter_P;
    flow_float* ducros=lim.ducros;

    flow_float *drodx=grd.drodx, *drody=grd.drody, *drodz=grd.drodz;
    flow_float *dUxdx=grd.dUxdx, *dUxdy=grd.dUxdy, *dUxdz=grd.dUxdz;
    flow_float *dUydx=grd.dUydx, *dUydy=grd.dUydy, *dUydz=grd.dUydz;
    flow_float *dUzdx=grd.dUzdx, *dUzdy=grd.dUzdy, *dUzdz=grd.dUzdz;
    flow_float *dPdx=grd.dPdx, *dPdy=grd.dPdy, *dPdz=grd.dPdz;
    // --- ローカル展開ここまで ---

    geom_int ip_orig = blockDim.x*blockIdx.x + threadIdx.x;

    if (ip_orig < nNormal_halo_Planes) {

        geom_int ip = normal_halo_planes_d[ip_orig];

        flow_float lam[5];
        flow_float R[5][5];
        flow_float L[5][5];
        flow_float dW[5];
        flow_float difQ1[5];
        flow_float difQ2[5];
        flow_float Q[5];
        int i, j;

        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];

        // 境界 ghost (ic1>=nCells): ghost には勾配/リミタが無いので 1 次化。
        if (ic1 >= nCells) conv_scheme = -1;

        geom_float f = fx[ip];
        geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];
        geom_float nx = sxx/sss, ny = syy/sss, nz = szz/sss;

        flow_float ccx_0 = ccx[ic0], ccy_0 = ccy[ic0], ccz_0 = ccz[ic0];
        flow_float ccx_1 = ccx[ic1], ccy_1 = ccy[ic1], ccz_1 = ccz[ic1];
        flow_float dcc_x = ccx_1 - ccx_0, dcc_y = ccy_1 - ccy_0, dcc_z = ccz_1 - ccz_0;
        flow_float dc0p_x = pcx[ip] - ccx_0, dc0p_y = pcy[ip] - ccy_0, dc0p_z = pcz[ip] - ccz_0;
        flow_float dc1p_x = pcx[ip] - ccx_1, dc1p_y = pcy[ip] - ccy_1, dc1p_z = pcz[ip] - ccz_1;

        const flow_float small_rho = 1.0e-4f;
        const flow_float small_p   = 1.0f;
        const flow_float small_a2  = 1.0f;

        // Ducros 連動リミタ (ROE と同じ扱い)
        flow_float duc = max(min(max(ducros[ic0], ducros[ic1]), 1.0), 0.0);
        flow_float limiter_ro_L = apply_ducros_limiter(limiter_ro[ic0], duc);
        flow_float limiter_Ux_L = apply_ducros_limiter(limiter_Ux[ic0], duc);
        flow_float limiter_Uy_L = apply_ducros_limiter(limiter_Uy[ic0], duc);
        flow_float limiter_Uz_L = apply_ducros_limiter(limiter_Uz[ic0], duc);
        flow_float limiter_P_L  = apply_ducros_limiter(limiter_P[ic0], duc);
        flow_float limiter_ro_R = apply_ducros_limiter(limiter_ro[ic1], duc);
        flow_float limiter_Ux_R = apply_ducros_limiter(limiter_Ux[ic1], duc);
        flow_float limiter_Uy_R = apply_ducros_limiter(limiter_Uy[ic1], duc);
        flow_float limiter_Uz_R = apply_ducros_limiter(limiter_Uz[ic1], duc);
        flow_float limiter_P_R  = apply_ducros_limiter(limiter_P[ic1], duc);

        // ---- MUSCL 再構成した L/R 状態 (散逸用) ----
        flow_float ro_L = interp_dispatch(conv_scheme, limit_scheme, ro[ic0], ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro_L);
        flow_float Ux_L = interp_dispatch(conv_scheme, limit_scheme, Ux[ic0], Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux_L);
        flow_float Uy_L = interp_dispatch(conv_scheme, limit_scheme, Uy[ic0], Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy_L);
        flow_float Uz_L = interp_dispatch(conv_scheme, limit_scheme, Uz[ic0], Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz_L);
        flow_float P_L  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic0], Ps[ic1], dPdx[ic0], dPdy[ic0], dPdz[ic0], dPdx[ic1], dPdy[ic1], dPdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P_L);
        ro_L = max(ro_L, small_rho); P_L = max(P_L, small_p);
        flow_float roe_L = P_L/(ga-1.0) + 0.5*ro_L*(Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L);
        flow_float Ht_L  = roe_L/ro_L + P_L/ro_L;

        flow_float ro_R = interp_dispatch(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro_R);
        flow_float Ux_R = interp_dispatch(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux_R);
        flow_float Uy_R = interp_dispatch(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy_R);
        flow_float Uz_R = interp_dispatch(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz_R);
        flow_float P_R  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1], dPdy[ic1], dPdz[ic1], dPdx[ic0], dPdy[ic0], dPdz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P_R);
        ro_R = max(ro_R, small_rho); P_R = max(P_R, small_p);
        flow_float roe_R = P_R/(ga-1.0) + 0.5*ro_R*(Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R);
        flow_float Ht_R  = roe_R/ro_R + P_R/ro_R;

        flow_float U_L = (Ux_L*sxx + Uy_L*syy + Uz_L*szz)/sss;
        flow_float U_R = (Ux_R*sxx + Uy_R*syy + Uz_R*szz)/sss;

        // ---- Roe 特性行列散逸 (L 状態評価) ----
        flow_float ua = Ux_L, va = Uy_L, wa = Uz_L, Ha = Ht_L;
        flow_float ca = sqrt(max((ga-1.0)*(Ha-0.5*(ua*ua+va*va+wa*wa)), small_a2));
        flow_float Ua = U_L;

        flow_float P_ro  = 0.5*(ga-1.0)*(ua*ua + va*va + wa*wa);
        flow_float P_roe = ga-1.0;
        flow_float P_rou = -ua*(ga-1.0);
        flow_float P_rov = -va*(ga-1.0);
        flow_float P_row = -wa*(ga-1.0);
        flow_float z = ua*ua+va*va+wa*wa - P_ro/P_roe;
        flow_float fai = P_ro-ca*ca;

        R[0][0]=1.0;      R[0][1]=nx;               R[0][2]=ny;               R[0][3]=nz;               R[0][4]=1.0;
        R[1][0]=Ha-ca*Ua; R[1][1]=z*nx+va*nz-wa*ny; R[1][2]=z*ny+wa*nx-ua*nz; R[1][3]=z*nz+ua*ny-va*nx; R[1][4]=Ha+ca*Ua;
        R[2][0]=ua-ca*nx; R[2][1]=ua*nx;            R[2][2]=ua*ny-nz;         R[2][3]=ua*nz+ny;         R[2][4]=ua+ca*nx;
        R[3][0]=va-ca*ny; R[3][1]=va*nx+nz;         R[3][2]=va*ny;            R[3][3]=va*nz-nx;         R[3][4]=va+ca*ny;
        R[4][0]=wa-ca*nz; R[4][1]=wa*nx-ny;         R[4][2]=wa*ny+nx;         R[4][3]=wa*nz;            R[4][4]=wa+ca*nz;

        L[0][0]=0.5*(P_ro+ca*Ua);            L[0][1]=P_roe*0.5; L[0][2]=(P_rou-ca*nx)*0.5; L[0][3]=(P_rov-ca*ny)*0.5; L[0][4]=(P_row-ca*nz)*0.5;
        L[1][0]=-fai*nx-ca*ca*(va*nz-wa*ny); L[1][1]=-P_roe*nx; L[1][2]=-P_rou*nx;         L[1][3]=-P_rov*nx+ca*ca*nz; L[1][4]=-P_row*nx-ca*ca*ny;
        L[2][0]=-fai*ny-ca*ca*(wa*nx-ua*nz); L[2][1]=-P_roe*ny; L[2][2]=-P_rou*ny-ca*ca*nz; L[2][3]=-P_rov*ny;         L[2][4]=-P_row*ny+ca*ca*nx;
        L[3][0]=-fai*nz-ca*ca*(ua*ny-va*nx); L[3][1]=-P_roe*nz; L[3][2]=-P_rou*nz+ca*ca*ny; L[3][3]=-P_rov*nz-ca*ca*nx; L[3][4]=-P_row*nz;
        L[4][0]=0.5*(P_ro-ca*Ua); L[4][1]=0.5*P_roe; L[4][2]=0.5*(P_rou+ca*nx); L[4][3]=0.5*(P_rov+ca*ny); L[4][4]=0.5*(P_row+ca*nz);

        Q[0]=ro_L; Q[1]=ro_L*Ht_L-P_L; Q[2]=ro_L*Ux_L; Q[3]=ro_L*Uy_L; Q[4]=ro_L*Uz_L;
        for (i=0;i<5;i++) for (j=0;j<5;j++) L[i][j] = L[i][j]/(ca*ca);
        lam[0]=abs(Ua-ca); lam[1]=abs(Ua); lam[2]=abs(Ua); lam[3]=abs(Ua); lam[4]=abs(Ua+ca);
        for (i=0;i<5;i++) { dW[i]=0.0; for (j=0;j<5;j++) dW[i]+=L[i][j]*Q[j]; }
        for (i=0;i<5;i++) difQ1[i]=lam[i]*dW[i];
        for (i=0;i<5;i++) { difQ2[i]=0.0; for (j=0;j<5;j++) difQ2[i]+=R[i][j]*difQ1[j]; }

        flow_float res_ro_temp   = 0.5*difQ2[0]*sss;
        flow_float res_roe_temp  = 0.5*difQ2[1]*sss;
        flow_float res_roUx_temp = 0.5*difQ2[2]*sss;
        flow_float res_roUy_temp = 0.5*difQ2[3]*sss;
        flow_float res_roUz_temp = 0.5*difQ2[4]*sss;

        // ---- Roe 特性行列散逸 (R 状態評価) ----
        ua=Ux_R; va=Uy_R; wa=Uz_R; Ha=Ht_R;
        ca=sqrt(max((ga-1.0)*(Ha-0.5*(ua*ua+va*va+wa*wa)), small_a2));
        Ua=U_R;
        P_ro=0.5*(ga-1.0)*(ua*ua+va*va+wa*wa); P_roe=ga-1.0;
        P_rou=-ua*(ga-1.0); P_rov=-va*(ga-1.0); P_row=-wa*(ga-1.0);
        z=ua*ua+va*va+wa*wa - P_ro/P_roe; fai=P_ro-ca*ca;

        R[0][0]=1.0;      R[0][1]=nx;               R[0][2]=ny;               R[0][3]=nz;               R[0][4]=1.0;
        R[1][0]=Ha-ca*Ua; R[1][1]=z*nx+va*nz-wa*ny; R[1][2]=z*ny+wa*nx-ua*nz; R[1][3]=z*nz+ua*ny-va*nx; R[1][4]=Ha+ca*Ua;
        R[2][0]=ua-ca*nx; R[2][1]=ua*nx;            R[2][2]=ua*ny-nz;         R[2][3]=ua*nz+ny;         R[2][4]=ua+ca*nx;
        R[3][0]=va-ca*ny; R[3][1]=va*nx+nz;         R[3][2]=va*ny;            R[3][3]=va*nz-nx;         R[3][4]=va+ca*ny;
        R[4][0]=wa-ca*nz; R[4][1]=wa*nx-ny;         R[4][2]=wa*ny+nx;         R[4][3]=wa*nz;            R[4][4]=wa+ca*nz;

        L[0][0]=0.5*(P_ro+ca*Ua);            L[0][1]=P_roe*0.5; L[0][2]=(P_rou-ca*nx)*0.5; L[0][3]=(P_rov-ca*ny)*0.5; L[0][4]=(P_row-ca*nz)*0.5;
        L[1][0]=-fai*nx-ca*ca*(va*nz-wa*ny); L[1][1]=-P_roe*nx; L[1][2]=-P_rou*nx;         L[1][3]=-P_rov*nx+ca*ca*nz; L[1][4]=-P_row*nx-ca*ca*ny;
        L[2][0]=-fai*ny-ca*ca*(wa*nx-ua*nz); L[2][1]=-P_roe*ny; L[2][2]=-P_rou*ny-ca*ca*nz; L[2][3]=-P_rov*ny;         L[2][4]=-P_row*ny+ca*ca*nx;
        L[3][0]=-fai*nz-ca*ca*(ua*ny-va*nx); L[3][1]=-P_roe*nz; L[3][2]=-P_rou*nz+ca*ca*ny; L[3][3]=-P_rov*nz-ca*ca*nx; L[3][4]=-P_row*nz;
        L[4][0]=0.5*(P_ro-ca*Ua); L[4][1]=0.5*P_roe; L[4][2]=0.5*(P_rou+ca*nx); L[4][3]=0.5*(P_rov+ca*ny); L[4][4]=0.5*(P_row+ca*nz);

        Q[0]=ro_R; Q[1]=ro_R*Ht_R-P_R; Q[2]=ro_R*Ux_R; Q[3]=ro_R*Uy_R; Q[4]=ro_R*Uz_R;
        for (i=0;i<5;i++) for (j=0;j<5;j++) L[i][j] = L[i][j]/(ca*ca);
        lam[0]=abs(Ua-ca); lam[1]=abs(Ua); lam[2]=abs(Ua); lam[3]=abs(Ua); lam[4]=abs(Ua+ca);
        for (i=0;i<5;i++) { dW[i]=0.0; for (j=0;j<5;j++) dW[i]+=L[i][j]*Q[j]; }
        for (i=0;i<5;i++) difQ1[i]=lam[i]*dW[i];
        for (i=0;i<5;i++) { difQ2[i]=0.0; for (j=0;j<5;j++) difQ2[i]+=R[i][j]*difQ1[j]; }

        res_ro_temp   -= 0.5*difQ2[0]*sss;
        res_roe_temp  -= 0.5*difQ2[1]*sss;
        res_roUx_temp -= 0.5*difQ2[2]*sss;
        res_roUy_temp -= 0.5*difQ2[3]*sss;
        res_roUz_temp -= 0.5*difQ2[4]*sss;

        // ---- KEEP 中心流束 (隣接対の生値; KE/エントロピー保存) ----
        flow_float Ctilde  = 0.5*(ro[ic0]+ro[ic1])*0.5*( (Ux[ic0]+Ux[ic1])*nx
                                                        +(Uy[ic0]+Uy[ic1])*ny
                                                        +(Uz[ic0]+Uz[ic1])*nz );
        flow_float Mtildex = Ctilde*(Ux[ic0]+Ux[ic1])*0.5;
        flow_float Mtildey = Ctilde*(Uy[ic0]+Uy[ic1])*0.5;
        flow_float Mtildez = Ctilde*(Uz[ic0]+Uz[ic1])*0.5;
        flow_float Ktilde  = Ctilde*0.5*(Ux[ic0]*Ux[ic1] +Uy[ic0]*Uy[ic1] +Uz[ic0]*Uz[ic1]);
        flow_float Itilde  = Ctilde*0.5*(Ps[ic0]/ro[ic0] +Ps[ic1]/ro[ic1])/(ga-1.0);
        flow_float Gtildex = 0.5*(Ps[ic0]+Ps[ic1])*nx;
        flow_float Gtildey = 0.5*(Ps[ic0]+Ps[ic1])*ny;
        flow_float Gtildez = 0.5*(Ps[ic0]+Ps[ic1])*nz;
        flow_float Ptilde  = 0.5*( (Ux[ic0]*Ps[ic1] + Ux[ic1]*Ps[ic0])*nx
                                  +(Uy[ic0]*Ps[ic1] + Uy[ic1]*Ps[ic0])*ny
                                  +(Uz[ic0]*Ps[ic1] + Uz[ic1]*Ps[ic0])*nz );

        res_ro_temp   = Ctilde*sss                     + res_ro_temp;
        res_roe_temp  = (Ktilde + Itilde + Ptilde)*sss + res_roe_temp;
        res_roUx_temp = (Mtildex + Gtildex)*sss        + res_roUx_temp;
        res_roUy_temp = (Mtildey + Gtildey)*sss        + res_roUy_temp;
        res_roUz_temp = (Mtildez + Gtildez)*sss        + res_roUz_temp;

        // スカラー輸送用の面質量流束 (連続式と整合: 散逸込みの総流束 ic0->ic1)
        massflux[ip] = res_ro_temp;

        atomicAdd(&res_ro[ic0]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic0], -res_roUx_temp);
        atomicAdd(&res_roUy[ic0], -res_roUy_temp);
        atomicAdd(&res_roUz[ic0], -res_roUz_temp);
        atomicAdd(&res_roe[ic0] , -res_roe_temp);

        atomicAdd(&res_ro[ic1]  , res_ro_temp);
        atomicAdd(&res_roUx[ic1], res_roUx_temp);
        atomicAdd(&res_roUy[ic1], res_roUy_temp);
        atomicAdd(&res_roUz[ic1], res_roUz_temp);
        atomicAdd(&res_roe[ic1] , res_roe_temp);
    }

    __syncthreads();
}
