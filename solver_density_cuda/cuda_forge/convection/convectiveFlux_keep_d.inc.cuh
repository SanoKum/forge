// =============================================================================
// 対流フラックス KEEP_d の実装断片。
//   - 純粋 KEEP (Kinetic Energy & Entropy Preserving) 中心流束。**中心部は散逸項なし**。
//     隣接セル/ノード対 (ic0,ic1) の生値で中心流束を構成する (KE/エントロピー保存)。
//     Roe 行列散逸・MUSCL 再構成・リミタ・Ducros は持たない (低散逸 LES/ILES 用; SGS は WALE が担う)。
//   - opt-in 散逸レイヤ (keepDissType, plans/active/convection-keep-es-dissipation.md):
//     keepDissType==0 で従来どおり散逸ゼロ (ビット不変)。==1 で scalar entropy-stable 散逸
//     F -= 0.5*σ*λ'*ΔU を全 5 式に加算 (λ'=|Un|+c'; lowMachPrecond>=1 なら c'=lowMachCprime 低マッハ
//     スケール、else c)。LLF 型は Δw·ΔU>=0 (エントロピー凸性) で ES。純 KEEP は低マッハ圧力 2Δ 市松を
//     原理的に減衰できない (中心平均の null-mode) ため、その減衰はこのレイヤが担う。
//   - convectiveFlux_d.cu に textual include される断片であり **スタンドアロン TU ではない**。
//     CMake の source には足さないこと。
//   - 周回面は geom.nLoopPlanes (= convPlaneBound)。cell では内部+境界 ghost を周回し、境界面は
//     ic1=ghost セルの生値で中心流束を組む。node 弱形式では内部双対面のみ (境界は別カーネル/gather)。
// =============================================================================

__global__ void KEEP_d
(
 flow_float ga,
 int keepDissType, flow_float keepDissCoeff,      // opt-in ES 散逸レイヤ (0: off=ビット不変)
 int lowMachPrecond, flow_float precondEps,       // λ' の音速スケール (>=1 で c'=lowMachCprime)
 FaceGeom      geom,
 PrimState     st,
 ResidualOut   reso
)
{
    // --- struct 引数をローカルへ展開 ---
    const geom_int nCells       = geom.nCells;
    geom_int*      plane_cells  = geom.plane_cells;
    const geom_int nLoopPlanes  = geom.nLoopPlanes;
    geom_int*      loop_planes  = geom.loop_planes;
    geom_float *sx=geom.sx, *sy=geom.sy, *sz=geom.sz, *ss=geom.ss;
    flow_float* massflux=geom.massflux;

    flow_float *ro=st.ro, *Ux=st.Ux, *Uy=st.Uy, *Uz=st.Uz, *Ps=st.Ps;

    flow_float *res_ro=reso.res_ro, *res_roUx=reso.res_roUx, *res_roUy=reso.res_roUy, *res_roUz=reso.res_roUz, *res_roe=reso.res_roe;
    // --- ローカル展開ここまで ---

    geom_int ip_orig = blockDim.x*blockIdx.x + threadIdx.x;

    if (ip_orig < nLoopPlanes) {

        geom_int ip = loop_planes[ip_orig];

        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];
        (void)nCells;

        geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];
        geom_float nx = sxx/sss, ny = syy/sss, nz = szz/sss;

        // ---- 純粋 KEEP 中心流束 (隣接対の生値; KE/エントロピー保存・散逸なし) ----
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

        flow_float res_ro_temp   = Ctilde*sss;
        flow_float res_roe_temp  = (Ktilde + Itilde + Ptilde)*sss;
        flow_float res_roUx_temp = (Mtildex + Gtildex)*sss;
        flow_float res_roUy_temp = (Mtildey + Gtildey)*sss;
        flow_float res_roUz_temp = (Mtildez + Gtildez)*sss;

        // ---- opt-in scalar entropy-stable 散逸レイヤ (keepDissType==1) ----
        //   F -= 0.5*σ*λ'*ΔU (ΔU=保存量ジャンプ ic0→ic1)。LLF 型は Δw·ΔU>=0 で ES。
        //   λ'=|Un|+c': lowMachPrecond>=1 なら c'=lowMachCprime (低マッハで c'→O(|u|)、
        //   M>=1 で c'=c に復帰)。過大 σ は解像乱流を殺すため TGV (L2) で較正する。
        if (keepDissType == 1) {
            // 保存量ジャンプ (primitives から構成; CPG: roe = P/(γ-1) + ½ρ|u|²)
            const flow_float dro   = ro[ic1] - ro[ic0];
            const flow_float droUx = ro[ic1]*Ux[ic1] - ro[ic0]*Ux[ic0];
            const flow_float droUy = ro[ic1]*Uy[ic1] - ro[ic0]*Uy[ic0];
            const flow_float droUz = ro[ic1]*Uz[ic1] - ro[ic0]*Uz[ic0];
            const flow_float roe0  = Ps[ic0]/(ga-1.0)
                                   + 0.5*ro[ic0]*(Ux[ic0]*Ux[ic0]+Uy[ic0]*Uy[ic0]+Uz[ic0]*Uz[ic0]);
            const flow_float roe1  = Ps[ic1]/(ga-1.0)
                                   + 0.5*ro[ic1]*(Ux[ic1]*Ux[ic1]+Uy[ic1]*Uy[ic1]+Uz[ic1]*Uz[ic1]);
            const flow_float droe  = roe1 - roe0;

            // 面代表の波速 λ' = |Un| + c(') (面平均状態で評価)
            const flow_float roF = 0.5*(ro[ic0]+ro[ic1]);
            const flow_float PsF = 0.5*(Ps[ic0]+Ps[ic1]);
            const flow_float uxF = 0.5*(Ux[ic0]+Ux[ic1]);
            const flow_float uyF = 0.5*(Uy[ic0]+Uy[ic1]);
            const flow_float uzF = 0.5*(Uz[ic0]+Uz[ic1]);
            const flow_float Un  = uxF*nx + uyF*ny + uzF*nz;
            const flow_float c   = sqrt(ga*PsF/roF);
            const flow_float cd  = (lowMachPrecond >= 1)
                                 ? lowMachCprime(c, sqrt(uxF*uxF+uyF*uyF+uzF*uzF), Un, precondEps)
                                 : c;
            const flow_float coef = 0.5*keepDissCoeff*(fabs(Un) + cd)*sss;

            res_ro_temp   -= coef*dro;
            res_roUx_temp -= coef*droUx;
            res_roUy_temp -= coef*droUy;
            res_roUz_temp -= coef*droUz;
            res_roe_temp  -= coef*droe;
        }

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
