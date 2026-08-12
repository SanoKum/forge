#!/usr/bin/env python3
"""SU2 の `flow.vtu` から「実測/壁法則」比を出す — forge 側と**同一定義**の外部クロスチェック。

`rep_normal_interp.py` は forge の HDF5 (`res_*.h5` + `rep_*` 診断) 専用なので、
SU2 の結果を同じ土俵に載せるための独立ツール。

## 何をするか (forge 側と揃えてある点)

1. **$u_\\tau$ は独立に決める** — 壁第一節点の**分子粘性勾配**から
   $u_\\tau=\\sqrt{\\nu\\,\\partial U/\\partial y|_w}$。
   **壁法則を逆解きして $u_\\tau$ を作ることはしない** (それをやると比が恒等的に 1 になる)。
2. **同じ $y_{\\rm rep}$ へ要素ベース補間** (quad→2 三角形・重心座標、外挿しない)。
3. 速度**ベクトル**を補間してから接線成分 $U_t=|U-(U\\cdot n)n|$ を取る。
4. 壁法則予測 $U^+_{\\rm law}(y^+)\\cdot u_\\tau$ と比べる。

## 使い方

    python3 tools/su2_wall_law_ratio.py <run>/flow.vtu --yrep 1.7681e-4 \\
        --stations 0.451,0.602,0.694,0.800 [--law reichardt|spalding|log]

`--yrep` は比較したい forge 側代表点高さ (case/26 ny52 では 1.7681e-4 m)。
"""
import re, struct, argparse
import numpy as np

KAPPA = 0.41


def u_plus(yp, law="reichardt"):
    if law == "reichardt":
        return (np.log(1.0 + KAPPA * yp) / KAPPA
                + 7.8 * (1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0)))
    if law == "log":
        return np.log(yp) / KAPPA + 5.0
    if law == "spalding":          # y+ = f(u+) を Newton で反転
        up = np.log(max(yp, 1e-9)) / KAPPA + 5.0
        for _ in range(200):
            e = KAPPA * up
            f = up + np.exp(-KAPPA * 5.0) * (np.exp(e) - 1 - e - e**2 / 2 - e**3 / 6) - yp
            df = 1 + np.exp(-KAPPA * 5.0) * KAPPA * (np.exp(e) - 1 - e - e**2 / 2)
            s = f / df
            up -= s
            if abs(s) < 1e-13:
                break
        return up
    raise SystemExit(f"未知の壁法則: {law}")


def read_vtu(path):
    """appended raw (header_type=UInt64, 非圧縮) の VTU を読む。"""
    raw = open(path, "rb").read()
    hdr = raw[:raw.index(b"<AppendedData")].decode("utf8", "ignore")
    off0 = raw.index(b"_", raw.index(b"<AppendedData")) + 1
    tmap = {"Float32": "f4", "Float64": "f8", "Int32": "i4", "Int64": "i8", "UInt8": "u1"}
    out = {}
    for t, n, nc, off in re.findall(
            r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents=\s*"(\d+)" offset="(\d+)"', hdr):
        o = off0 + int(off)
        nb = struct.unpack("<Q", raw[o:o + 8])[0]
        a = np.frombuffer(raw[o + 8:o + 8 + nb], dtype=np.dtype("<" + tmap[t]))
        out[n if n else "pts"] = a.reshape(-1, int(nc)) if int(nc) > 1 else a
    return out


class TriMesh:
    """要素ベース補間 (quad→2 三角形、重心座標)。包含要素が無ければ None を返す。"""

    def __init__(self, P, conn, offs):
        tris, s = [], 0
        for o in offs:
            nd = conn[s:o]
            s = o
            if len(nd) == 4:
                tris += [[nd[0], nd[1], nd[2]], [nd[0], nd[2], nd[3]]]
            elif len(nd) == 3:
                tris.append(list(nd))
        self.P, self.tri = P, np.array(tris)
        from scipy.spatial import cKDTree
        self.tree = cKDTree(self.P[self.tri].mean(axis=1))

    def interp(self, p, f, ncand=64):
        _, cand = self.tree.query(p, k=min(ncand, len(self.tri)))
        for t in np.atleast_1d(cand):
            a, b, c = self.P[self.tri[t]]
            v0, v1, v2 = b - a, c - a, p - a
            den = v0[0] * v1[1] - v1[0] * v0[1]
            if abs(den) < 1e-30:
                continue
            u = (v2[0] * v1[1] - v1[0] * v2[1]) / den
            v = (v0[0] * v2[1] - v2[0] * v0[1]) / den
            if u >= -1e-9 and v >= -1e-9 and u + v <= 1.0 + 1e-9:
                w = np.array([1.0 - u - v, u, v])
                return np.array([float(np.dot(w, f[:, k][self.tri[t]])) for k in range(f.shape[1])])
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtu")
    ap.add_argument("--yrep", type=float, required=True, help="比較する代表点高さ [m]")
    ap.add_argument("--stations", default="0.451,0.602,0.694,0.800")
    ap.add_argument("--law", default="reichardt", choices=["reichardt", "spalding", "log"])
    a = ap.parse_args()

    d = read_vtu(a.vtu)
    P = d["pts"][:, :2]
    T = TriMesh(P, d["connectivity"], d["offsets"])
    U = np.column_stack([d["Velocity"][:, 0], d["Velocity"][:, 1]])
    ro, mu = d["Density"], d["Laminar_Viscosity"]
    x, y = P[:, 0], P[:, 1]
    xa = np.unique(np.round(x[x > 1e-6], 6))

    print(f"# {a.vtu}")
    print(f"# u_tau は壁第一節点の分子勾配から独立に決める (壁法則の逆解きはしない)")
    print(f"# y_rep={a.yrep:.4e} m へ要素ベース補間、壁法則={a.law}")
    print()
    print(f'{"x":>7s} {"Re_theta":>9s} {"u_tau":>8s} {"y+":>7s} '
          f'{"実測U_t":>9s} {"法則予測":>9s} {"実測/則":>8s}')
    for x0 in [float(v) for v in a.stations.split(",")]:
        xc = xa[np.argmin(abs(xa - x0))]
        i = np.sort(np.where(np.abs(x - xc) < 1e-4)[0])
        i = i[np.argsort(y[i])]
        j = i[1] if y[i[0]] < 1e-12 else i[0]          # 壁の 1 つ上
        nu = mu[j] / ro[j]
        ut = np.sqrt(nu * U[j, 0] / y[j])              # 分子勾配 (壁で u=0)
        # 運動量厚さ
        ue = U[i, 0].max()
        th = np.trapz(ro[i] * U[i, 0] / (ro[i].max() * ue) * (1 - U[i, 0] / ue), y[i])
        ret = ro[i].max() * ue * th / mu[i][-1]
        V = T.interp(np.array([xc, a.yrep]), U)
        if V is None:
            print(f"{xc:7.3f}  (包含要素なし — 外挿になるので除外)")
            continue
        Ut = abs(V[0])                                 # 平板なので n=(0,1) → 接線成分 = |U_x|
        yp = ut * a.yrep / nu
        pred = u_plus(yp, a.law) * ut
        print(f"{xc:7.3f} {ret:9.0f} {ut:8.4f} {yp:7.2f} {Ut:9.3f} {pred:9.3f} {Ut/pred:8.4f}")


if __name__ == "__main__":
    main()
