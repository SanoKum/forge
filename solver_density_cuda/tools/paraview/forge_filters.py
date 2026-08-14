"""forge の res_*.xmf / res_*.h5 用 ParaView Python プラグイン。

提供フィルタ:
  Forge Derived Quantities  … マッハ数 / シュリーレン (無次元密度勾配) /
                              Q 値 / ヘリシティ を一括で計算して出力に追加する。

読み込み方法 (ParaView GUI):
  Tools > Manage Plugins > Load New... で本ファイルを選び、
  必要なら "Auto Load" を有効にする。以後 Filters > Alphabetical に
  "Forge Derived Quantities" が現れる。
  "Failed to import paraview.detail.pythonalgorithm" で失敗する環境
  (ParaView 5.11 + Python 3.12 の組み合わせ。ParaView 側が Python 3.11 で削除された
  inspect.getargspec を import するため) では、代わりに同ディレクトリの
  macro_load_forge_filters.py を Macros > Add new macro... で登録し、そのボタンから読み込む。

pvpython / pvbatch から使う場合:
  from paraview.simple import *
  LoadPlugin("<repo>/solver_density_cuda/tools/paraview/forge_filters.py", ns=globals())
  src = XDMFReader(FileNames=["res_12000.xmf"])
  drv = ForgeDerivedQuantities(Input=src)

入力に期待する配列 (cell モードは Cell Data、node モードは Point Data):
  必須: ro, Ux, Uy, Uz          (速度ベクトル U と各種微分量の元)
  任意: P, sonic                (マッハ数)
  任意: dUxdx…dUzdz, drodx…drodz (ソルバ自身の勾配。あれば既定で優先使用)
勾配配列が無い場合は vtkGradientFilter でメッシュから勾配を計算する。

出力配列:
  U, Mach, sound_speed
  grad_ro, schlieren_mag, schlieren_dir, schlieren
  vorticity, vorticity_mag, Q, Q_norm, (lambda2)
  helicity, helicity_norm
"""

import inspect

# ParaView 5.11 系の paraview.detail.pythonalgorithm は Python 3.11+ で
# 削除された inspect.getargspec を import するため、ここで補っておく
# (使われ方は `.args` の参照のみなので getfullargspec で代替できる)。
if not hasattr(inspect, "getargspec"):
    inspect.getargspec = inspect.getfullargspec

import numpy as np

from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa  # noqa: F401  (ParaView 側の初期化用)
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkDataObject
from vtkmodules.vtkFiltersGeneral import vtkGradientFilter

POINTS = vtkDataObject.FIELD_ASSOCIATION_POINTS
CELLS = vtkDataObject.FIELD_ASSOCIATION_CELLS

VEL_GRAD_NAMES = (
    "dUxdx", "dUxdy", "dUxdz",
    "dUydx", "dUydy", "dUydz",
    "dUzdx", "dUzdy", "dUzdz",
)
RO_GRAD_NAMES = ("drodx", "drody", "drodz")


# ---------------------------------------------------------------- 補助関数

def _attr(ds, assoc):
    return ds.GetPointData() if assoc == POINTS else ds.GetCellData()


def _has(ds, assoc, name):
    return _attr(ds, assoc).GetArray(name) is not None


def _get(ds, assoc, name):
    arr = _attr(ds, assoc).GetArray(name)
    if arr is None:
        return None
    return vtk_to_numpy(arr).astype(np.float64)


def _add(ds, assoc, name, values):
    values = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    varr = numpy_to_vtk(values, deep=1)
    varr.SetName(name)
    _attr(ds, assoc).AddArray(varr)


def _detect_assoc(ds):
    """ro / Ux が乗っている側を自動判定する (cell モード=Cell, node モード=Point)。"""
    for assoc in (POINTS, CELLS):
        if _has(ds, assoc, "ro") or _has(ds, assoc, "Ux"):
            return assoc
    # 判定不能ならデータ数の多い方
    return POINTS if ds.GetPointData().GetNumberOfArrays() >= ds.GetCellData().GetNumberOfArrays() else CELLS


def _velocity(ds, assoc):
    ux, uy, uz = (_get(ds, assoc, n) for n in ("Ux", "Uy", "Uz"))
    if ux is None:
        vec = _get(ds, assoc, "U")
        if vec is not None and vec.ndim == 2 and vec.shape[1] == 3:
            return vec
        return None
    n = ux.shape[0]
    uy = np.zeros(n) if uy is None else uy
    uz = np.zeros(n) if uz is None else uz
    return np.column_stack((ux, uy, uz))


def _mesh_gradient(ds, assoc, values, name="__forge_tmp"):
    """vtkGradientFilter でメッシュから勾配を計算する。

    スカラー入力なら (N,3)、3 成分ベクトル入力なら (N,9) を返す。
    (N,9) の並びは行優先 = (du/dx, du/dy, du/dz, dv/dx, …) で
    forge の dUxdx… と同じ規約。
    """
    tmp = ds.NewInstance()
    tmp.ShallowCopy(ds)
    _add(tmp, assoc, name, values)

    gf = vtkGradientFilter()
    gf.SetInputData(tmp)
    gf.SetInputArrayToProcess(0, 0, 0, assoc, name)
    gf.SetResultArrayName("__forge_grad")
    gf.Update()

    out = gf.GetOutput()
    grad = vtk_to_numpy(_attr(out, assoc).GetArray("__forge_grad")).astype(np.float64)
    return grad


def _velocity_gradient(ds, assoc, vel, use_solver):
    """速度勾配テンソル g[:, i, j] = dU_i/dx_j を返す。"""
    if use_solver and all(_has(ds, assoc, n) for n in VEL_GRAD_NAMES):
        cols = [_get(ds, assoc, n) for n in VEL_GRAD_NAMES]
        return np.column_stack(cols).reshape(-1, 3, 3), "solver"
    return _mesh_gradient(ds, assoc, vel).reshape(-1, 3, 3), "mesh"


def _density_gradient(ds, assoc, use_solver):
    if use_solver and all(_has(ds, assoc, n) for n in RO_GRAD_NAMES):
        return np.column_stack([_get(ds, assoc, n) for n in RO_GRAD_NAMES]), "solver"
    ro = _get(ds, assoc, "ro")
    if ro is None:
        return None, None
    return _mesh_gradient(ds, assoc, ro), "mesh"


def _leaf_pairs(inp, out):
    """(入力 leaf, 出力 leaf) の組を返す。composite / 単一データセット両対応。"""
    if inp.IsA("vtkCompositeDataSet"):
        out.CopyStructure(inp)
        pairs = []
        it = inp.NewIterator()
        it.InitTraversal()
        while not it.IsDoneWithTraversal():
            blk = it.GetCurrentDataObject()
            if blk is not None:
                new_blk = blk.NewInstance()
                new_blk.ShallowCopy(blk)
                out.SetDataSet(it, new_blk)
                pairs.append((blk, new_blk))
            it.GoToNextItem()
        return pairs
    out.ShallowCopy(inp)
    return [(inp, out)]


# ---------------------------------------------------------------- フィルタ本体

@smproxy.filter(label="Forge Derived Quantities")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet", "vtkCompositeDataSet"], composite_data_supported=True)
class ForgeDerivedQuantities(VTKPythonAlgorithmBase):
    """forge の結果からマッハ数・シュリーレン・Q 値・ヘリシティを計算する。"""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=1,
            inputType="vtkDataObject", outputType="vtkDataObject")
        self._gamma = 1.4
        self._use_solver_grad = True
        self._use_solver_sonic = True
        self._mach = True
        self._schlieren = True
        self._qcrit = True
        self._helicity = True
        self._lambda2 = False
        self._sch_dir = [1.0, 0.0, 0.0]
        self._sch_exp = 15.0

    # -- properties -------------------------------------------------

    @smproperty.doublevector(name="Gamma", default_values=1.4)
    @smdomain.doublerange(min=1.0, max=2.0)
    def SetGamma(self, value):
        """比熱比 (音速 a=sqrt(gamma*P/ro) 用)。"""
        self._gamma = float(value)
        self.Modified()

    @smproperty.intvector(name="UseSolverGradients", label="Use Solver Gradients", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetUseSolverGradients(self, value):
        """ON: ソルバが出力した勾配 (dUxdx…, drodx…) をそのまま使う。
        OFF (または配列が無い場合): vtkGradientFilter でメッシュから勾配を計算する。"""
        self._use_solver_grad = bool(value)
        self.Modified()

    @smproperty.intvector(name="UseSolverSonic", label="Use Solver Sound Speed", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetUseSolverSonic(self, value):
        """ON: 配列 sonic をそのまま音速に使う。OFF/無い場合は a=sqrt(gamma*P/ro)。"""
        self._use_solver_sonic = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeMach", label="Compute Mach", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeMach(self, value):
        self._mach = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeSchlieren", label="Compute Schlieren", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeSchlieren(self, value):
        self._schlieren = bool(value)
        self.Modified()

    @smproperty.doublevector(name="SchlierenDirection", label="Schlieren Direction",
                             default_values=[1.0, 0.0, 0.0], number_of_elements=3)
    def SetSchlierenDirection(self, x, y, z):
        """方向シュリーレン schlieren_dir = (n^*grad ro)/max|n^*grad ro| の方向 n。"""
        self._sch_dir = [float(x), float(y), float(z)]
        self.Modified()

    @smproperty.doublevector(name="SchlierenExponent", label="Schlieren Exponent", default_values=15.0)
    @smdomain.doublerange(min=0.0, max=200.0)
    def SetSchlierenExponent(self, value):
        """schlieren = exp(-k*|grad ro|/max|grad ro|) の k。大きいほど弱い勾配を強調。"""
        self._sch_exp = float(value)
        self.Modified()

    @smproperty.intvector(name="ComputeQ", label="Compute Q", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeQ(self, value):
        self._qcrit = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeLambda2", label="Compute Lambda2", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeLambda2(self, value):
        self._lambda2 = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeHelicity", label="Compute Helicity", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeHelicity(self, value):
        self._helicity = bool(value)
        self.Modified()

    # -- pipeline ---------------------------------------------------

    def RequestDataObject(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0])
        if inp is None:
            return 0
        out = vtkDataObject.GetData(outInfo)
        if out is None or not out.IsA(inp.GetClassName()):
            out = inp.NewInstance()
            outInfo.GetInformationObject(0).Set(vtkDataObject.DATA_OBJECT(), out)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0])
        out = vtkDataObject.GetData(outInfo)

        pairs = _leaf_pairs(inp, out)

        # 1 パス目: 各 leaf の派生量を計算 (シュリーレン正規化用に全体最大も集める)
        results = []
        gmax_mag = 0.0
        gmax_dir = 0.0
        for _, leaf in pairs:
            if not leaf.IsA("vtkDataSet"):
                results.append(None)
                continue
            res = self._compute_leaf(leaf)
            results.append(res)
            if res is not None and "grad_ro_mag" in res:
                gmax_mag = max(gmax_mag, float(res["grad_ro_mag"].max(initial=0.0)))
                gmax_dir = max(gmax_dir, float(np.abs(res["grad_ro_dir"]).max(initial=0.0)))

        # 2 パス目: 正規化してから書き出す
        for (_, leaf), res in zip(pairs, results):
            if res is None:
                continue
            assoc = res["assoc"]
            for name, values in res["arrays"].items():
                _add(leaf, assoc, name, values)
            if "grad_ro_mag" in res:
                mag = res["grad_ro_mag"]
                dir_ = res["grad_ro_dir"]
                norm_mag = mag / gmax_mag if gmax_mag > 0.0 else np.zeros_like(mag)
                norm_dir = dir_ / gmax_dir if gmax_dir > 0.0 else np.zeros_like(dir_)
                _add(leaf, assoc, "schlieren_mag", norm_mag)
                _add(leaf, assoc, "schlieren_dir", norm_dir)
                _add(leaf, assoc, "schlieren", np.exp(-self._sch_exp * norm_mag))
        return 1

    # -- 中身 --------------------------------------------------------

    def _compute_leaf(self, ds):
        n = ds.GetNumberOfPoints()
        if n == 0:
            return None
        assoc = _detect_assoc(ds)
        if assoc == CELLS and ds.GetNumberOfCells() == 0:
            return None

        res = {"assoc": assoc, "arrays": {}}
        arrays = res["arrays"]

        vel = _velocity(ds, assoc)
        if vel is not None:
            arrays["U"] = vel

        ro = _get(ds, assoc, "ro")

        # --- マッハ数
        if self._mach and vel is not None:
            sonic = _get(ds, assoc, "sonic") if self._use_solver_sonic else None
            if sonic is None:
                pres = _get(ds, assoc, "P")
                if pres is not None and ro is not None:
                    sonic = np.sqrt(np.maximum(self._gamma * pres / np.maximum(ro, 1e-30), 0.0))
            if sonic is not None:
                umag = np.linalg.norm(vel, axis=1)
                arrays["sound_speed"] = sonic
                arrays["Mach"] = umag / np.maximum(sonic, 1e-30)

        # --- シュリーレン (正規化は RequestData 側の全体最大で行う)
        if self._schlieren and ro is not None:
            grad_ro, _ = _density_gradient(ds, assoc, self._use_solver_grad)
            if grad_ro is not None:
                arrays["grad_ro"] = grad_ro
                d = np.asarray(self._sch_dir, dtype=np.float64)
                dnorm = np.linalg.norm(d)
                d = d / dnorm if dnorm > 0.0 else np.array([1.0, 0.0, 0.0])
                res["grad_ro_mag"] = np.linalg.norm(grad_ro, axis=1)
                res["grad_ro_dir"] = grad_ro @ d

        # --- 速度勾配ベースの量 (Q / lambda2 / ヘリシティ)
        need_grad = (self._qcrit or self._helicity or self._lambda2) and vel is not None
        if need_grad:
            g, _src = _velocity_gradient(ds, assoc, vel, self._use_solver_grad)

            # 渦度 omega = rot U
            omega = np.column_stack((
                g[:, 2, 1] - g[:, 1, 2],
                g[:, 0, 2] - g[:, 2, 0],
                g[:, 1, 0] - g[:, 0, 1],
            ))
            omega_mag = np.linalg.norm(omega, axis=1)
            arrays["vorticity"] = omega
            arrays["vorticity_mag"] = omega_mag

            if self._qcrit or self._lambda2:
                gT = np.transpose(g, (0, 2, 1))
                S = 0.5 * (g + gT)
                W = 0.5 * (g - gT)
                s2 = np.einsum("nij,nij->n", S, S)
                w2 = np.einsum("nij,nij->n", W, W)
                if self._qcrit:
                    arrays["Q"] = 0.5 * (w2 - s2)
                    denom = 0.5 * (w2 + s2)
                    arrays["Q_norm"] = np.where(denom > 0.0, 0.5 * (w2 - s2) / np.maximum(denom, 1e-30), 0.0)
                if self._lambda2:
                    M = np.matmul(S, S) + np.matmul(W, W)
                    M = 0.5 * (M + np.transpose(M, (0, 2, 1)))  # 数値的な非対称を除去
                    # 発散した run では場に NaN が混ざる。固有値計算が落ちないよう
                    # 有限な要素だけ解き、残りは NaN のままにする。
                    ok = np.isfinite(M).all(axis=(1, 2))
                    lam2 = np.full(M.shape[0], np.nan)
                    if ok.any():
                        lam2[ok] = np.linalg.eigvalsh(M[ok])[:, 1]  # 昇順の中間固有値
                    arrays["lambda2"] = lam2

            if self._helicity:
                hel = np.einsum("ni,ni->n", vel, omega)
                arrays["helicity"] = hel
                denom = np.linalg.norm(vel, axis=1) * omega_mag
                arrays["helicity_norm"] = np.where(denom > 0.0, hel / np.maximum(denom, 1e-30), 0.0)

        return res
