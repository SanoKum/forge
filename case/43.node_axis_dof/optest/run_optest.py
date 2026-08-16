#!/usr/bin/env python3
"""node 軸対称の離散連続式 (質量) 演算子テスト。

場: ρ=const, u = U0 - 2 a x, v = a r, p = const  →  ∂x u + (1/r)∂r(r v) = -2a + 2a = 0。
   軸対称では連続式残差が恒等的に 0。planar 扱いなら残差 = -ρ a (0 にならない)。
測定: 陽解法 Euler 1 step (固定 dt) で  res_ro/V = (ρ¹-ρ⁰)/dt  を各 CV で取り、
      e = (res_ro/V)/(ρ a) を無次元残差とする。軸 CV / 第一内点行 / 内部 CV を分けて h 収束を見る。
可変: メッシュ幅 h, convMethod, gradLSQ, nodeMidpointFx, 軸ノード IC 規約 (node: v_A=0 / centroid: v_A=a·r̄)。
"""
import argparse, itertools, json, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np, h5py

W = Path("/home/sano/work/forge-axisnode")
sys.path.insert(0, str(W / "design"))
from forge_design.meshing.mesh2d import write_msh41_2d  # noqa: E402

BUILD = W / "solver_density_cuda/build"
TOOLS = W / "solver_density_cuda/tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
HERE = Path(__file__).resolve().parent

RHO0, U0, A, P0 = 1.2, 1200.0, 500.0, 101325.0   # u∈[200,1200] on x∈[0,1], v=a r
GAMMA, CP = 1.4, 1004.5
L, R = 1.0, 0.2
DT = 1.0e-3   # 陽解法 1 step は dt 任意 (線形評価)。float32 の Δρ 分解能を稼ぐため大きく取る


def duct_mesh(h):
    ni, nj = int(round(L / h)) + 1, int(round(R / h)) + 1
    xs, rs = np.linspace(0, L, ni), np.linspace(0, R, nj)
    X, Rr = np.meshgrid(xs, rs, indexing="ij")
    coords = np.zeros((ni * nj, 3)); coords[:, 0] = X.ravel(); coords[:, 1] = Rr.ravel()
    nid = lambda i, j: i * nj + j
    quads = np.array([(nid(i, j), nid(i+1, j), nid(i+1, j+1), nid(i, j+1))
                      for i in range(ni-1) for j in range(nj-1)], dtype=np.int64)
    bedges = {"axis": [(nid(i, 0), nid(i+1, 0)) for i in range(ni-1)],
              "wall": [(nid(i, nj-1), nid(i+1, nj-1)) for i in range(ni-1)],
              "inlet": [(nid(0, j), nid(0, j+1)) for j in range(nj-1)],
              "outlet": [(nid(ni-1, j), nid(ni-1, j+1)) for j in range(nj-1)]}
    return coords, quads, bedges


def solver_cfg(conv, lsq, midfx, extra_mesh=""):
    return f"""mesh: {{meshFormat: "hdf5", discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, nodeAxisDirichlet: 0, nodeMidpointFx: {midfx}, gradLSQ: {lsq}{extra_mesh}, meshFileName: "duct.h5", valueFileName: "duct.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: {CP}, gamma: {GAMMA}}}
time:
  unsteady: 1
  dualTime: 0
  last: {{control: 0, nStepOuter: 1}}
  deltaT: {{control: 0, dt: {DT}, cfl: 1.0, cfl_pseudo: 1.0, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 0, lowMachPrecond: 0, detectNaN: 0}}
  outStepStart: 0
  outStepInterval: 1
  timeIntegration: 1
  nStepInner: 1
space: {{convMethod: {conv}, limiter: 0}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
"""
# limiter: 0 は意図的 (線形場の再構成を厳密にするため; Venkat は線形場でも ψ=3/4 に落ちる)。

BCOND = """inlet:  {physID: 1, kind: slip, outputHDFflg: 0, ints: , floats: }
outlet: {physID: 2, kind: slip, outputHDFflg: 0, ints: , floats: }
wall:   {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }
axis:   {physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }
"""


def write_ic(h5, axis_conv):
    with h5py.File(h5, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
        n = cc.shape[0]
        assert n == nc.shape[0], "node 変換でない"
        x, r = nc[:, 0].copy(), nc[:, 1].copy()   # 値の位置 = ノード
        if axis_conv == "centroid":
            r = cc[:, 1].copy()                    # 値の位置 = 双対重心 (r̄)
            x = cc[:, 0].copy()
        u, v = U0 - 2*A*x, A*r
        ro = np.full(n, RHO0)
        roe = P0/(GAMMA-1) + 0.5*ro*(u*u + v*v)
        for k, val in {"ro": ro, "roUx": ro*u, "roUy": ro*v, "roUz": 0*ro, "roe": roe}.items():
            f["/VALUE/"+k][:] = val.astype(np.float32)
        return nc, cc


def run_one(tag, h, conv, lsq, midfx, axis_conv, extra_mesh="", force=False):
    rd = HERE / f"h{h:g}_{tag}"
    if rd.exists() and not force and (rd / "res_1.h5").exists():
        return rd
    rd.mkdir(exist_ok=True)
    (rd / "solverConfig.yaml").write_text(solver_cfg(conv, lsq, midfx, extra_mesh))
    (rd / "bcondConfig.yaml").write_text(BCOND)
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    coords, quads, bedges = duct_mesh(h)
    write_msh41_2d(rd / "duct.msh", coords, quads, bedges)
    subprocess.run([str(BUILD / "convertGmshToForge"), "duct.msh", "duct.h5"], cwd=rd, env=ENV,
                   check=True, capture_output=True, text=True)
    write_ic(rd / "duct.h5", axis_conv)
    for f in rd.glob("res_*"): f.unlink()
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
    if not (rd / "res_1.h5").exists():
        print(r.stdout[-2000:], r.stderr[-2000:]); raise RuntimeError(f"{rd}: res_1.h5 無し")
    return rd


def evaluate(rd):
    with h5py.File(rd / "duct.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3); cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(rd / "res_0.h5") as f0, h5py.File(rd / "res_1.h5") as f1:
        ro0 = f0["/VALUE/ro"][:].astype(np.float64); ro1 = f1["/VALUE/ro"][:].astype(np.float64)
    e = (ro1 - ro0) / DT / (RHO0 * A)
    x, r = nc[:, 0], nc[:, 1]
    rpos = np.unique(np.round(r, 9))
    r1 = rpos[1]
    inner_x = (x > 0.2 + 1e-9) & (x < 0.8 - 1e-9)
    axis = inner_x & (r < 1e-9)
    row1 = inner_x & (np.abs(r - r1) < 1e-9)
    inter = inner_x & (r > r1 + 1e-9) & (r < 0.15)
    st = lambda m: (float(np.abs(e[m]).max()), float(np.sqrt(np.mean(e[m]**2))))
    return {"axis": st(axis), "row1": st(row1), "interior": st(inter), "n": int(nc.shape[0]),
            "ccy_axis": float(cc[axis, 1].mean()), "r1": float(r1),
            "e_axis_mean": float(e[axis].mean()), "e_row1_mean": float(e[row1].mean())}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--hs", default="0.02,0.01,0.005")
    ap.add_argument("--variants", default="all")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    hs = [float(s) for s in a.hs.split(",")]
    # (tag, conv, lsq, midfx, axis_conv)
    V = [("c1_gg_fx0_node", 1, 0, 0, "node"),
         ("c1_gg_fx1_node", 1, 0, 1, "node"),
         ("c1_lsq_fx0_node", 1, 2, 0, "node"),
         ("c1_lsq_fx1_node", 1, 2, 1, "node"),
         ("c1_gg_fx0_cent", 1, 0, 0, "centroid"),
         ("c1_gg_fx1_cent", 1, 0, 1, "centroid"),
         ("c1_lsq_fx1_cent", 1, 2, 1, "centroid"),
         ("c0_fx0_node", 0, 0, 0, "node"),
         ("c0_fx0_cent", 0, 0, 0, "centroid"),
         # 値位置=ノード座標 (nodeValueAtNode=1, 試作) — 値規約は node が整合
         ("van_c1_gg_fx0_node", 1, 0, 0, "node", ", nodeValueAtNode: 1"),
         ("van_c1_lsq_fx0_node", 1, 2, 0, "node", ", nodeValueAtNode: 1"),
         ("van_c0_fx0_node", 0, 0, 0, "node", ", nodeValueAtNode: 1")]
    if a.variants != "all":
        V = [v for v in V if v[0] in a.variants.split(",")]
    out = {}
    for v in V:
        tag, conv, lsq, midfx, ac = v[:5]
        extra = v[5] if len(v) > 5 else ""
        for h in hs:
            rd = run_one(tag, h, conv, lsq, midfx, ac, extra_mesh=extra, force=a.force)
            out[f"{tag}@h{h:g}"] = evaluate(rd)
            m = out[f"{tag}@h{h:g}"]
            print(f"{tag:18s} h={h:<6g} axis max {m['axis'][0]:.3e} (mean {m['e_axis_mean']:+.3e}) | row1 max {m['row1'][0]:.3e} | interior max {m['interior'][0]:.3e}  [ccy_axis={m['ccy_axis']:.4g}, r1={m['r1']:.3g}]")
    (HERE / "optest_results.json").write_text(json.dumps(out, indent=1))
