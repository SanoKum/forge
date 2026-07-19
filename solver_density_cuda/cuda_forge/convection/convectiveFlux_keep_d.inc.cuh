// =============================================================================
// 対流フラックス KEEP_d の実装断片。
//   - 純粋 KEEP (Kinetic Energy & Entropy Preserving) 中心流束。**中心部は散逸項なし**。
//     隣接セル/ノード対 (ic0,ic1) の生値で中心流束を構成する (KE/エントロピー保存)。
//     Roe 行列散逸・MUSCL 再構成・リミタ・Ducros は持たない (低散逸 LES/ILES 用; SGS は WALE が担う)。
//   - opt-in 散逸レイヤ (keepDissType, plans/accepted/convection-keep-es-dissipation.md):
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
 int thermalMethod, const SpeciesThermo* sp,      // Step 3/4: TP (thermalMethod==2)。CPG では未使用
 int nSpecies, flow_float** roY,                  // Step 4: 多成分の組成 (nSpecies<=1 なら未参照)
 int keepDissType, flow_float keepDissCoeff,      // opt-in ES 散逸レイヤ (0: off=ビット不変)
 int keepDissCprime, flow_float precondEps,       // 散逸波速: 1 で音響 c'=lowMachCprime (lowMachPrecond から独立)
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

    flow_float *ro=st.ro, *Ux=st.Ux, *Uy=st.Uy, *Uz=st.Uz, *Ps=st.Ps, *roe_c=st.roe;

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
        flow_float Itilde;
        if (thermalMethod == 2) {
            // TP (Step 3): 内部エネルギーを保存量 roe から取る (EOS 非依存・保存整合)。
            // CPG の P/ρ/(γ-1) は TP では e(T) と不一致になるため使わない。
            const flow_float e0c = roe_c[ic0]/ro[ic0]
                                 - 0.5*(Ux[ic0]*Ux[ic0]+Uy[ic0]*Uy[ic0]+Uz[ic0]*Uz[ic0]);
            const flow_float e1c = roe_c[ic1]/ro[ic1]
                                 - 0.5*(Ux[ic1]*Ux[ic1]+Uy[ic1]*Uy[ic1]+Uz[ic1]*Uz[ic1]);
            Itilde = Ctilde*0.5*(e0c + e1c);
        } else {
            Itilde = Ctilde*0.5*(Ps[ic0]/ro[ic0] +Ps[ic1]/ro[ic1])/(ga-1.0);   // CPG (ビット不変)
        }
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

        // ---- opt-in entropy-stable 散逸レイヤ (keepDissType 1: scalar / 2: matrix) ----
        //   1 (scalar): F -= 0.5*σ*λ'*ΔU。LLF 型は Δw·ΔU>=0 で ES。全成分同一 λ' なので
        //     市松と一緒に渦も食う (TGV KE cost 2.7%@σ=0.05)。
        //   2 (matrix): F -= 0.5*σ*R|Λ'|S Rᵀ Δw (entropy-scaled Roe 型, Chandrashekar/Barth)。
        //     w=∂η/∂U (η=-ρs/(γ-1)), H=∂U/∂w=RSRᵀ (S: 音響 ρ/(2γ)・エントロピー ρ(γ-1)/γ・
        //     せん断 p — tools/verify_entropy_scaling.py で数値検証済)。
        //     |Λ'| は音響のみ |Un|+c'、せん断/エントロピーは |Un| → 渦を食わず市松 (音響) を減衰。
        //   λ'/c': lowMachPrecond>=1 なら c'=lowMachCprime (低マッハで c'→O(|u|)、M>=1 で c に復帰)。
        //   過大 σ は解像乱流を殺すため TGV (L2) で較正する。
        if (keepDissType == 1) {
            // 保存量ジャンプ (primitives から構成; CPG: roe = P/(γ-1) + ½ρ|u|²)
            const flow_float dro   = ro[ic1] - ro[ic0];
            const flow_float droUx = ro[ic1]*Ux[ic1] - ro[ic0]*Ux[ic0];
            const flow_float droUy = ro[ic1]*Uy[ic1] - ro[ic0]*Uy[ic0];
            const flow_float droUz = ro[ic1]*Uz[ic1] - ro[ic0]*Uz[ic0];
            // TP では保存量 roe のジャンプそのもの (EOS 非依存)、CPG は従来式 (ビット不変)
            flow_float droe;
            if (thermalMethod == 2) {
                droe = roe_c[ic1] - roe_c[ic0];
            } else {
                const flow_float roe0  = Ps[ic0]/(ga-1.0)
                                       + 0.5*ro[ic0]*(Ux[ic0]*Ux[ic0]+Uy[ic0]*Uy[ic0]+Uz[ic0]*Uz[ic0]);
                const flow_float roe1  = Ps[ic1]/(ga-1.0)
                                       + 0.5*ro[ic1]*(Ux[ic1]*Ux[ic1]+Uy[ic1]*Uy[ic1]+Uz[ic1]*Uz[ic1]);
                droe = roe1 - roe0;
            }

            // 面代表の波速 λ' = |Un| + c(') (面平均状態で評価)
            const flow_float roF = 0.5*(ro[ic0]+ro[ic1]);
            const flow_float PsF = 0.5*(Ps[ic0]+Ps[ic1]);
            const flow_float uxF = 0.5*(Ux[ic0]+Ux[ic1]);
            const flow_float uyF = 0.5*(Uy[ic0]+Uy[ic1]);
            const flow_float uzF = 0.5*(Uz[ic0]+Uz[ic1]);
            const flow_float Un  = uxF*nx + uyF*ny + uzF*nz;
            const flow_float c   = (thermalMethod == 2)
                                 ? 0.5*(st.sonic[ic0]+st.sonic[ic1])   // TP: dependent vars の凍結音速
                                 : sqrt(ga*PsF/roF);                   // CPG (ビット不変)
            const flow_float cd  = (keepDissCprime == 1)
                                 ? lowMachCprime(c, sqrt(uxF*uxF+uyF*uyF+uzF*uzF), Un, precondEps)
                                 : c;
            const flow_float coef = 0.5*keepDissCoeff*(fabs(Un) + cd)*sss;

            res_ro_temp   -= coef*dro;
            res_roUx_temp -= coef*droUx;
            res_roUy_temp -= coef*droUy;
            res_roUz_temp -= coef*droUz;
            res_roe_temp  -= coef*droe;
        }
        else if (keepDissType == 2 && thermalMethod != 2) {
            // ---- matrix ES 散逸 (CPG): D = R|Λ'|S Rᵀ Δw ----
            // Δw (エントロピー変数ジャンプ): w=[(γ-s)/(γ-1)-β|u|², 2βu, -2β], β=ρ/(2p), s=ln p-γ ln ρ
            const flow_float b0 = ro[ic0]/(2.0*Ps[ic0]);
            const flow_float b1 = ro[ic1]/(2.0*Ps[ic1]);
            const flow_float q0 = Ux[ic0]*Ux[ic0]+Uy[ic0]*Uy[ic0]+Uz[ic0]*Uz[ic0];
            const flow_float q1 = Ux[ic1]*Ux[ic1]+Uy[ic1]*Uy[ic1]+Uz[ic1]*Uz[ic1];
            const flow_float ds = log(Ps[ic1]/Ps[ic0]) - ga*log(ro[ic1]/ro[ic0]);
            const flow_float dw0 = -ds/(ga-1.0) - (b1*q1 - b0*q0);
            const flow_float dw1 = 2.0*(b1*Ux[ic1] - b0*Ux[ic0]);
            const flow_float dw2 = 2.0*(b1*Uy[ic1] - b0*Uy[ic0]);
            const flow_float dw3 = 2.0*(b1*Uz[ic1] - b0*Uz[ic0]);
            const flow_float dw4 = -2.0*(b1 - b0);

            // 面平均状態と固有系 (python 検証済の正規化: r=[1,u∓cn,Ht∓c·Un] 等)
            const flow_float roF = 0.5*(ro[ic0]+ro[ic1]);
            const flow_float PsF = 0.5*(Ps[ic0]+Ps[ic1]);
            const flow_float uxF = 0.5*(Ux[ic0]+Ux[ic1]);
            const flow_float uyF = 0.5*(Uy[ic0]+Uy[ic1]);
            const flow_float uzF = 0.5*(Uz[ic0]+Uz[ic1]);
            const flow_float qF  = uxF*uxF+uyF*uyF+uzF*uzF;
            const flow_float Un  = uxF*nx + uyF*ny + uzF*nz;
            const flow_float c   = sqrt(ga*PsF/roF);
            const flow_float Ht  = ga*PsF/((ga-1.0)*roF) + 0.5*qF;

            // 接線基底 (n に対しロバストに構築)
            flow_float ax, ay, az;
            if (fabs(nx) < 0.9) { ax=1.0; ay=0.0; az=0.0; } else { ax=0.0; ay=1.0; az=0.0; }
            flow_float t1x = ny*az - nz*ay, t1y = nz*ax - nx*az, t1z = nx*ay - ny*ax;
            const flow_float t1n = sqrt(t1x*t1x+t1y*t1y+t1z*t1z);
            t1x/=t1n; t1y/=t1n; t1z/=t1n;
            const flow_float t2x = ny*t1z - nz*t1y, t2y = nz*t1x - nx*t1z, t2z = nx*t1y - ny*t1x;

            // 固有値: 音響のみ c(') 込み、せん断/エントロピーは |Un| (渦を守る)
            const flow_float cd  = (keepDissCprime == 1)
                                 ? lowMachCprime(c, sqrt(qF), Un, precondEps)
                                 : c;
            const flow_float lamA = fabs(Un) + cd;   // 音響 (±両波共通の安全上界)
            const flow_float lamU = fabs(Un);        // エントロピー・せん断

            // スケーリング S (Barth): 音響 ρ/(2γ), エントロピー ρ(γ-1)/γ, せん断 p
            const flow_float sA = roF/(2.0*ga);
            const flow_float sE = roF*(ga-1.0)/ga;
            const flow_float sS = PsF;

            // z_k = S_k |λ_k| (r_k · Δw)
            const flow_float ut1 = uxF*t1x + uyF*t1y + uzF*t1z;
            const flow_float ut2 = uxF*t2x + uyF*t2y + uzF*t2z;
            const flow_float rd1 = dw0 + (uxF-c*nx)*dw1 + (uyF-c*ny)*dw2 + (uzF-c*nz)*dw3 + (Ht-c*Un)*dw4; // r(un-c)·Δw
            const flow_float rd2 = dw0 + uxF*dw1 + uyF*dw2 + uzF*dw3 + 0.5*qF*dw4;                          // r(entropy)·Δw
            const flow_float rd3 = t1x*dw1 + t1y*dw2 + t1z*dw3 + ut1*dw4;                                   // r(shear1)·Δw
            const flow_float rd4 = t2x*dw1 + t2y*dw2 + t2z*dw3 + ut2*dw4;                                   // r(shear2)·Δw
            const flow_float rd5 = dw0 + (uxF+c*nx)*dw1 + (uyF+c*ny)*dw2 + (uzF+c*nz)*dw3 + (Ht+c*Un)*dw4; // r(un+c)·Δw
            const flow_float z1 = sA*lamA*rd1;
            const flow_float z2 = sE*lamU*rd2;
            const flow_float z3 = sS*lamU*rd3;
            const flow_float z4 = sS*lamU*rd4;
            const flow_float z5 = sA*lamA*rd5;

            // D = Σ z_k r_k、flux -= 0.5σ D S
            const flow_float coef = 0.5*keepDissCoeff*sss;
            const flow_float Dro   = z1 + z2 + z5;
            const flow_float DroUx = z1*(uxF-c*nx) + z2*uxF + z3*t1x + z4*t2x + z5*(uxF+c*nx);
            const flow_float DroUy = z1*(uyF-c*ny) + z2*uyF + z3*t1y + z4*t2y + z5*(uyF+c*ny);
            const flow_float DroUz = z1*(uzF-c*nz) + z2*uzF + z3*t1z + z4*t2z + z5*(uzF+c*nz);
            const flow_float Droe  = z1*(Ht-c*Un) + z2*0.5*qF + z3*ut1 + z4*ut2 + z5*(Ht+c*Un);

            res_ro_temp   -= coef*Dro;
            res_roUx_temp -= coef*DroUx;
            res_roUy_temp -= coef*DroUy;
            res_roUz_temp -= coef*DroUz;
            res_roe_temp  -= coef*Droe;
        }
        else if (keepDissType == 2) {
            // ---- matrix ES 散逸 (TP, Step 3/4): D = R|Λ'|S Rᵀ Δw ----
            // 次元付きエントロピー η=-ρs_mix, w=[(g-½|u|²)/T, u/T, -1/T] (Chalot-Hughes-Shakib),
            // S: 音響 ρ/(2γR_mix)・エントロピー ρ/cp_mix(T,Y)・せん断 ρT、エントロピー波 r_E=½q+e-cv·T
            // (tools/verify_entropy_scaling_tp.py + 凍結組成混合の数値検証で H=RSRᵀ・Ar=λr・SPD 確認済)。
            // 多成分 (Step 4): 混合物性 (R_mix/cp_mix/s0_mix) + 混合エントロピー -ΣY_kR_k ln X_k
            //   (Y ln X → 0 でゼロ濃度種でも正則)。組成一様なら凍結組成の TP 1 気体に厳密縮退。
            //   組成ジャンプ面では全系の厳密 ES は主張しない (Q は SPD なので無条件に散逸的=安定)。
            //   種輸送は散逸込み massflux 経由で連続式と整合 (種方程式自体への行列散逸は将来課題)。
            // thermo 評価は double (thermoHrefTemp datum 焼き込み係数をそのまま使用 → w も自動整合)。
            const int ns = (nSpecies > 1) ? nSpecies : 1;
            double Y0[THERMO_MAX_SPECIES], Y1[THERMO_MAX_SPECIES], YF[THERMO_MAX_SPECIES];
            if (nSpecies > 1) {
                // ★ Y はカーネル内で正規化 (Σ_k Y_k = 1 を強制)。ρY と ρ は別カーネル/atomicAdd で
                //   更新されるため Σ ρY_k ≠ ρ の共通モードノイズ (~1e-7) が乗り、これが s⁰ (×~10) と
                //   ln X (対数微分特異) で増幅されて w ノイズ→弱減衰モードへ注入され、市松減衰が
                //   プラトー化する (Y=[1,0] 縮退 bisect で確定)。正規化で共通モードを除去する。
                double sum0 = 0.0, sum1 = 0.0;
                for (int k = 0; k < ns; k++) {
                    Y0[k] = fmax((double)roY[k][ic0], 0.0);
                    Y1[k] = fmax((double)roY[k][ic1], 0.0);
                    sum0 += Y0[k]; sum1 += Y1[k];
                }
                const double inv0 = (sum0 > 1e-300) ? 1.0/sum0 : 0.0;
                const double inv1 = (sum1 > 1e-300) ? 1.0/sum1 : 0.0;
                for (int k = 0; k < ns; k++) {
                    Y0[k] *= inv0; Y1[k] *= inv1;
                    YF[k] = 0.5*(Y0[k]+Y1[k]);
                }
                if (sum0 <= 1e-300) { Y0[0] = 1.0; YF[0] = 0.5*(Y0[0]+Y1[0]); }
                if (sum1 <= 1e-300) { Y1[0] = 1.0; YF[0] = 0.5*(Y0[0]+Y1[0]); }
            } else { Y0[0] = 1.0; Y1[0] = 1.0; YF[0] = 1.0; }
            const double R0g = thermo_R_mix(sp, ns, Y0);
            const double R1g = thermo_R_mix(sp, ns, Y1);
            const double q0d = (double)Ux[ic0]*Ux[ic0]+(double)Uy[ic0]*Uy[ic0]+(double)Uz[ic0]*Uz[ic0];
            const double q1d = (double)Ux[ic1]*Ux[ic1]+(double)Uy[ic1]*Uy[ic1]+(double)Uz[ic1]*Uz[ic1];
            const double T0d = (double)Ps[ic0]/((double)ro[ic0]*R0g);
            const double T1d = (double)Ps[ic1]/((double)ro[ic1]*R1g);
            const double PsFd = 0.5*((double)Ps[ic0]+(double)Ps[ic1]);
            // 混合エントロピー項 -Σ Y_k R_k ln X_k (単成分では 0)
            double mixent0 = 0.0, mixent1 = 0.0;
            if (nSpecies > 1) {
                double mol0 = 0.0, mol1 = 0.0;
                for (int k = 0; k < ns; k++) { mol0 += Y0[k]/sp[k].MW; mol1 += Y1[k]/sp[k].MW; }
                for (int k = 0; k < ns; k++) {
                    const double Rk = THERMO_RU/sp[k].MW;
                    if (Y0[k] > 0.0) mixent0 -= Y0[k]*Rk*log(fmax(Y0[k]/sp[k].MW/mol0, 1e-300));
                    if (Y1[k] > 0.0) mixent1 -= Y1[k]*Rk*log(fmax(Y1[k]/sp[k].MW/mol1, 1e-300));
                }
            }
            // s_i (基準圧 pref=PsF: log 引数 ~1 で桁落ち回避。pref は Δw で相殺)
            const double s0d = thermo_s0_mix(sp, ns, Y0, T0d) + mixent0 - R0g*log((double)Ps[ic0]/PsFd);
            const double s1d = thermo_s0_mix(sp, ns, Y1, T1d) + mixent1 - R1g*log((double)Ps[ic1]/PsFd);
            const double g0d = thermo_h_mix(sp, ns, Y0, T0d) - T0d*s0d;
            const double g1d = thermo_h_mix(sp, ns, Y1, T1d) - T1d*s1d;
            const flow_float dw0 = (flow_float)((g1d-0.5*q1d)/T1d - (g0d-0.5*q0d)/T0d);
            const flow_float dw1 = (flow_float)((double)Ux[ic1]/T1d - (double)Ux[ic0]/T0d);
            const flow_float dw2 = (flow_float)((double)Uy[ic1]/T1d - (double)Uy[ic0]/T0d);
            const flow_float dw3 = (flow_float)((double)Uz[ic1]/T1d - (double)Uz[ic0]/T0d);
            const flow_float dw4 = (flow_float)(-(1.0/T1d - 1.0/T0d));

            // 面平均状態 (γ,cp,cv は面平均温度で凍結評価)
            const flow_float roF = 0.5*(ro[ic0]+ro[ic1]);
            const flow_float uxF = 0.5*(Ux[ic0]+Ux[ic1]);
            const flow_float uyF = 0.5*(Uy[ic0]+Uy[ic1]);
            const flow_float uzF = 0.5*(Uz[ic0]+Uz[ic1]);
            const flow_float qF  = uxF*uxF+uyF*uyF+uzF*uzF;
            const flow_float Un  = uxF*nx + uyF*ny + uzF*nz;
            const double TFd  = 0.5*(T0d+T1d);
            const double RFg  = thermo_R_mix(sp, ns, YF);
            const double cpFd = thermo_cp_mix(sp, ns, YF, TFd);
            const double cvFd = cpFd - RFg;
            const double gaFd = cpFd/cvFd;
            const flow_float c  = (flow_float)sqrt(gaFd*RFg*TFd);     // 凍結音速 (面平均 T, 面平均組成)
            // Ht・e は保存量から (EOS 厳密・datum 整合)
            const flow_float Ht = 0.5*( (roe_c[ic0]+Ps[ic0])/ro[ic0] + (roe_c[ic1]+Ps[ic1])/ro[ic1] );
            const flow_float eF = 0.5*( (flow_float)(roe_c[ic0]/ro[ic0]-0.5*q0d)
                                       +(flow_float)(roe_c[ic1]/ro[ic1]-0.5*q1d) );
            const flow_float r2E = 0.5*qF + eF - (flow_float)(cvFd*TFd);  // エントロピー波エネルギー成分

            // 接線基底
            flow_float ax, ay, az;
            if (fabs(nx) < 0.9) { ax=1.0; ay=0.0; az=0.0; } else { ax=0.0; ay=1.0; az=0.0; }
            flow_float t1x = ny*az - nz*ay, t1y = nz*ax - nx*az, t1z = nx*ay - ny*ax;
            const flow_float t1n = sqrt(t1x*t1x+t1y*t1y+t1z*t1z);
            t1x/=t1n; t1y/=t1n; t1z/=t1n;
            const flow_float t2x = ny*t1z - nz*t1y, t2y = nz*t1x - nx*t1z, t2z = nx*t1y - ny*t1x;

            // 固有値 (音響のみ c(') 込み) と TP スケーリング
            const flow_float cd  = (keepDissCprime == 1)
                                 ? lowMachCprime(c, sqrt(qF), Un, precondEps)
                                 : c;
            const flow_float lamA = fabs(Un) + cd;
            const flow_float lamU = fabs(Un);
            const flow_float sA = (flow_float)((double)roF/(2.0*gaFd*RFg));
            const flow_float sE = (flow_float)((double)roF/cpFd);
            const flow_float sS = (flow_float)((double)roF*TFd);

            // z_k = S_k |λ_k| (r_k · Δw) → D = Σ z_k r_k (エントロピー波のみ r_E=r2E)
            const flow_float ut1 = uxF*t1x + uyF*t1y + uzF*t1z;
            const flow_float ut2 = uxF*t2x + uyF*t2y + uzF*t2z;
            const flow_float rd1 = dw0 + (uxF-c*nx)*dw1 + (uyF-c*ny)*dw2 + (uzF-c*nz)*dw3 + (Ht-c*Un)*dw4;
            const flow_float rd2 = dw0 + uxF*dw1 + uyF*dw2 + uzF*dw3 + r2E*dw4;
            const flow_float rd3 = t1x*dw1 + t1y*dw2 + t1z*dw3 + ut1*dw4;
            const flow_float rd4 = t2x*dw1 + t2y*dw2 + t2z*dw3 + ut2*dw4;
            const flow_float rd5 = dw0 + (uxF+c*nx)*dw1 + (uyF+c*ny)*dw2 + (uzF+c*nz)*dw3 + (Ht+c*Un)*dw4;
            const flow_float z1 = sA*lamA*rd1;
            const flow_float z2 = sE*lamU*rd2;
            const flow_float z3 = sS*lamU*rd3;
            const flow_float z4 = sS*lamU*rd4;
            const flow_float z5 = sA*lamA*rd5;

            const flow_float coef = 0.5*keepDissCoeff*sss;
            const flow_float Dro   = z1 + z2 + z5;
            const flow_float DroUx = z1*(uxF-c*nx) + z2*uxF + z3*t1x + z4*t2x + z5*(uxF+c*nx);
            const flow_float DroUy = z1*(uyF-c*ny) + z2*uyF + z3*t1y + z4*t2y + z5*(uyF+c*ny);
            const flow_float DroUz = z1*(uzF-c*nz) + z2*uzF + z3*t1z + z4*t2z + z5*(uzF+c*nz);
            const flow_float Droe  = z1*(Ht-c*Un) + z2*r2E + z3*ut1 + z4*ut2 + z5*(Ht+c*Un);

            res_ro_temp   -= coef*Dro;
            res_roUx_temp -= coef*DroUx;
            res_roUy_temp -= coef*DroUy;
            res_roUz_temp -= coef*DroUz;
            res_roe_temp  -= coef*Droe;
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
