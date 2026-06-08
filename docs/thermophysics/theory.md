# 多成分 thermally-perfect gas — 理論

## 1. 目的と適用範囲

本章は forge に導入する**多成分 thermally-perfect gas (熱的完全気体)** の理論的背景をまとめる。対象は**非反応**の多成分混合 (混合・希釈・多成分境界層など) であり、化学反応源項は扱わない (将来拡張点)。

従来の forge は単一成分・**熱量的完全気体 (calorically perfect gas, CPG)** のみを扱っていた。すなわち比熱 $c_p,\,c_v$ と比熱比 $\gamma$ は定数で、温度は内部エネルギーから線形に陽に求まる。本拡張では:

- 比熱 $c_p(T)$ を温度依存とする (**NASA-9 多項式**, CEA 準拠)。
- 複数化学種の質量分率 $Y_s$ を輸送し、混合物性を **ideal-gas mixing** で評価する。
- 粘性・熱伝導率・質量拡散係数を **kinetic theory** (Chapman-Enskog) で評価する。

`thermalMethod==0` (CPG) は完全に保持し、`thermalMethod==2` で本モデルを有効化する。

## 2. 状態方程式と熱力学

### 2.1 理想気体の状態方程式

各化学種・混合物とも理想気体とする。混合気体定数は質量分率 $Y_s$ から

$$ R_{\mathrm{mix}} = R_u \sum_s \frac{Y_s}{W_s} = \frac{R_u}{W_{\mathrm{mix}}}, \qquad \frac{1}{W_{\mathrm{mix}}} = \sum_s \frac{Y_s}{W_s} $$

ここで $R_u=8.314462618\ \mathrm{J/(mol\,K)}$、$W_s$ は化学種 $s$ の分子量 [kg/mol]。圧力は

$$ P = \rho R_{\mathrm{mix}} T. $$

### 2.2 NASA-9 多項式による比熱・エンタルピー・エントロピー

CEA (McBride & Gordon 2002) のネイティブ形式である 9 係数多項式を用いる。各化学種・温度域ごとに係数 $a_0\dots a_8$ をもち、モル基準で

$$ \frac{C_{p,s}}{R_u} = a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T + a_4 T^2 + a_5 T^3 + a_6 T^4 $$

$$ \frac{H_s}{R_u T} = -a_0 T^{-2} + a_1 \frac{\ln T}{T} + a_2 + \frac{a_3}{2} T + \frac{a_4}{3} T^2 + \frac{a_5}{4} T^3 + \frac{a_6}{5} T^4 + \frac{a_7}{T} $$

$$ \frac{S^{\circ}_s}{R_u} = -\frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln T + a_3 T + \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3 + \frac{a_6}{4} T^4 + a_8 $$

積分定数 $a_7,a_8$ により、エンタルピーは**標準生成エンタルピーを含む絶対エンタルピー基準**となる。係数は 2 温度域 ($T_{\mathrm{lo}}\le T<T_{\mathrm{mid}}$ と $T_{\mathrm{mid}}\le T\le T_{\mathrm{hi}}$、標準は 200/1000/6000 K) で切り替える。範囲外は端でクランプし、エンタルピーは $h(T)\approx h(T_c)+c_p(T_c)(T-T_c)$ と線形外挿して衝撃波での暴走を防ぐ。

質量基準は $c_{p,s}=C_{p,s}/W_s$、$h_s=H_s/W_s$。

### 2.3 ideal-gas mixing による混合物性

混合比熱・絶対エンタルピーは質量分率重み:

$$ c_{p,\mathrm{mix}}(T) = \sum_s Y_s\,c_{p,s}(T), \qquad h_{\mathrm{mix}}(T) = \sum_s Y_s\,h_s(T). $$

定容比熱・比熱比は

$$ c_{v,\mathrm{mix}} = c_{p,\mathrm{mix}} - R_{\mathrm{mix}}, \qquad \gamma_{\mathrm{mix}} = \frac{c_{p,\mathrm{mix}}}{c_{v,\mathrm{mix}}}. $$

比内部エネルギーと凍結音速は

$$ e_{\mathrm{mix}}(T) = h_{\mathrm{mix}}(T) - R_{\mathrm{mix}} T, \qquad a = \sqrt{\gamma_{\mathrm{mix}} R_{\mathrm{mix}} T}. $$

### 2.4 温度反転 (Newton 法)

保存変数からは比内部エネルギー $e = \rho e_t/\rho - \tfrac12|\mathbf{u}|^2$ が得られる。CPG では $T=e/c_v$ と陽に解けるが、TP では $e_{\mathrm{mix}}(T)$ が非線形なので Newton 反転で解く:

$$ f(T) = e_{\mathrm{mix}}(T) - e = 0, \qquad f'(T) = c_{v,\mathrm{mix}}(T) = c_{p,\mathrm{mix}}(T) - R_{\mathrm{mix}}. $$

$$ T^{(k+1)} = T^{(k)} - \frac{e_{\mathrm{mix}}(T^{(k)}) - e}{c_{v,\mathrm{mix}}(T^{(k)})}. $$

厳密微分 $c_{v,\mathrm{mix}}$ により 2 次収束する。初期値は前ステップの $T$ をウォームスタートに使い (定常で約 2 反復)、各反復で $T\in[T_{\min},T_{\max}]$ にクランプする。FP32 では多項式和・Newton が破綻するため、**反転と多項式評価はカーネル内で double** で行う。

## 3. 化学種輸送方程式

非反応・質量分率 $Y_s$ の保存形 ($\rho Y_s$):

$$ \frac{\partial (\rho Y_s)}{\partial t} + \nabla\!\cdot(\rho \mathbf{u} Y_s) = -\nabla\!\cdot \mathbf{J}_s, \qquad s=1,\dots,N $$

$\mathbf{J}_s$ は化学種拡散流束 (§5)。移流は**連続の式と同じ質量流束** $\dot m$ を再利用して上流化し、$\sum_s$ が連続の式と整合する (質量保存)。実現可能性として $Y_s\ge 0$ をフロアし、各ステージで $\sum_s Y_s = 1$ に再正規化する。

## 4. kinetic theory による輸送係数

### 4.1 純成分 (Chapman-Enskog)

換算温度 $T^*=k_B T/\varepsilon_s$ と Lennard-Jones パラメータ ($\sigma_s$ [Å], $\varepsilon_s/k_B$ [K]) を用い、衝突積分 $\Omega^{(2,2)*},\Omega^{(1,1)*}$ は Neufeld ら (1972) の曲線近似で評価する。純成分粘性・熱伝導率 (modified Eucken) は

$$ \mu_s = 2.6693\times10^{-6}\,\frac{\sqrt{W_s[\mathrm{g/mol}]\,T}}{\sigma_s^2\,\Omega^{(2,2)*}(T^*)}\ [\mathrm{Pa\,s}], \qquad k_s = \mu_s\,(1.32\,c_{v,s} + 1.77\,R_s). $$

二成分拡散係数は

$$ D_{ij} = 1.8583\times10^{-7}\,\frac{\sqrt{T^3\,(1/W_i+1/W_j)}}{P\,\sigma_{ij}^2\,\Omega^{(1,1)*}(T^*_{ij})}, \qquad \sigma_{ij}=\tfrac12(\sigma_i+\sigma_j),\ \varepsilon_{ij}=\sqrt{\varepsilon_i\varepsilon_j}. $$

### 4.2 混合則

粘性は Wilke、熱伝導率は Wassiljewa/Mason-Saxena:

$$ \mu_{\mathrm{mix}} = \sum_i \frac{X_i \mu_i}{\sum_j X_j \phi_{ij}}, \quad k_{\mathrm{mix}} = \sum_i \frac{X_i k_i}{\sum_j X_j \phi_{ij}}, \quad \phi_{ij} = \frac{[1+\sqrt{\mu_i/\mu_j}\,(W_j/W_i)^{1/4}]^2}{\sqrt{8(1+W_i/W_j)}}. $$

ここで $X_i$ はモル分率。質量拡散は混合平均 (Curtiss-Hirschfelder):

$$ D_{i,\mathrm{mix}} = \frac{1-X_i}{\sum_{j\ne i} X_j/D_{ij}}. $$

## 5. 化学種拡散とエネルギー結合

化学種拡散流束は $\mathbf{J}_i = -\rho D_{i,\mathrm{mix}} \nabla Y_i$。混合平均は $\sum_i \mathbf{J}_i\ne 0$ となるため、補正速度で $\sum_i \mathbf{J}_i=0$ を担保する (定数 Schmidt 数フォールバックも用意)。エネルギー方程式には**化学種拡散によるエンタルピー輸送** $\sum_i h_i \mathbf{J}_i$ を熱伝導 $-k\nabla T$ に加える:

$$ \mathbf{q} = -k\,\nabla T + \sum_i h_i(T)\,\mathbf{J}_i. $$

## 6. 数値スキームとの整合 (対流流束)

密度ベース流束 (SLAU, Roe 等) は CPG では被移流全エンタルピーを $\gamma/(\gamma-1)\,P/\rho + \tfrac12|\mathbf u|^2$ と再構成するが、これは $c_p T$ に等しく、TP の絶対エンタルピー $h(T)=\int c_p\,dT'$ とは**一致しない**。TP では面の全エンタルピーを NASA の $h(T_{\mathrm{face}})$ で再構成する (全エネルギー保存を決定する最重要項)。音速は混合整合な $a=\sqrt{\gamma_{\mathrm{mix}}R_{\mathrm{mix}}T}$ を用いる。詳細は [`docs/convection/`](../convection/) および実装解説を参照。

## 参考文献

- B. J. McBride, M. J. Gordon, *NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species*, NASA/TP-2002-211556 (2002).
- R. A. Svehla, *Estimated Viscosities and Thermal Conductivities of Gases at High Temperatures*, NASA TR R-132 (1962).
- P. D. Neufeld, A. R. Janzen, R. A. Aziz, *J. Chem. Phys.* 57, 1100 (1972).
- C. R. Wilke, *J. Chem. Phys.* 18, 517 (1950).
- M. Vinokur, J.-P. Montagné, *J. Comput. Phys.* 89, 276 (1990) — 一般 EOS の Roe 平均。
