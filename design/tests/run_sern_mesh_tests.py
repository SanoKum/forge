#!/usr/bin/env python3
"""⑤ SERN 2 バンド構造メッシュ (スリットカウル) の整合テスト (S2)。"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec  # noqa: E402
from forge_design.meshing.mesh_sern import PHYS_SERN, SernMeshParams, generate_sern_mesh, write_msh41_named  # noqa: E402

FAIL = 0
def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1

k = PlanarMOC(SernKernelSpec(M_in=2.5, theta_r0=np.deg2rad(15.0), theta_c0=np.deg2rad(5.0), L_cowl=1.0,
                             x_max=6.0, nj=201, dx=4e-3, p_ext_over_p_in=0.05)).march()
d = k.design_ramp(M_c=3.3, f=0.45)
prm = SernMeshParams(ni_up=8, ni_noz=40, ni_plume=60, nj_top=31, nj_bot=21, scale=0.1, interface_angle=float(k.TH[-1, 0]),
                     top_ext_angle=d.info["theta_e"])
coords, quads, bedges, info, y_mid = generate_sern_mesh(d, prm)
ni, njt, njb, ite = info["ni"], info["nj_top"], info["nj_bot"], info["i_te"]
check("セル数 = (ni−1)(nj_top−1 + nj_bot−1)", quads.shape[0] == (ni - 1) * (njt - 1 + njb - 1))
check("ノード数 = ni·nj_bot + ni·(nj_top−1) + 重複 station 数", coords.shape[0] == ni * njb + ni * (njt - 1) + len(info["dup_stations"]))
# 反時計回り (符号付き面積 > 0)
xy = coords[:, :2]
a = np.zeros(len(quads))
for kk in range(4):
    p0, p1 = xy[quads[:, kk]], xy[quads[:, (kk + 1) % 4]]
    a += p0[:, 0] * p1[:, 1] - p1[:, 0] * p0[:, 1]
check("全 quad が反時計回り・非退化", np.all(a > 1e-14), f"min area {a.min():.2e}")
# 境界辺の数
check("cowl_in と cowl_out の辺数 = i_te", len(bedges["cowl_in"]) == ite and len(bedges["cowl_out"]) == ite)
check("ramp + top_out の辺数 = ni−1", len(bedges["ramp"]) + len(bedges["top_out"]) == ni - 1)
check("outlet 辺数 = nj_top−1 + nj_bot−1", len(bedges["outlet"]) == njt - 1 + njb - 1)
# スリット: cowl_in/cowl_out のノードは座標一致・ID 相異
ci = np.array(bedges["cowl_in"]); co = np.array(bedges["cowl_out"])
check("スリット: 上下ノードの座標一致", np.allclose(coords[ci], coords[co]))
check("スリット: 上下ノード ID は別 (TE 点のみ共有)", not np.any(ci[:, 0] == co[:, 0]) and ci[-1, 1] == co[-1, 1])
# 全境界辺の総和 = 閉境界 (各セル辺で共有されない辺の数)
from collections import Counter
cnt = Counter()
for q in quads:
    for kk in range(4):
        e = tuple(sorted((int(q[kk]), int(q[(kk + 1) % 4]))))
        cnt[e] += 1
free = {e for e, c in cnt.items() if c == 1}
tagged = {tuple(sorted(map(int, e))) for nm in bedges for e in bedges[nm]}
check("境界辺タグの集合 = quad の非共有辺の集合", free == tagged, f"free {len(free)} tagged {len(tagged)}")
# ランプ輪郭がメッシュ上境界に一致
rx = coords[[e[0] for e in bedges["ramp"]], 0] / prm.scale
ry = coords[[e[0] for e in bedges["ramp"]], 1] / prm.scale
m = rx >= 0
check("ランプ上境界 = 設計輪郭 (補間 1e-9)", np.allclose(ry[m], np.interp(rx[m], d.ramp_xy[:, 0], d.ramp_xy[:, 1]), atol=1e-9))
out = Path(__file__).with_name("_sern_mesh_test.msh"); write_msh41_named(out, coords, quads, bedges, PHYS_SERN)
txt = out.read_text(); check("msh4.1 書き出し (PhysicalNames 9 個)", "$PhysicalNames\n9" in txt); out.unlink()
print("ALL PASS" if FAIL == 0 else f"{FAIL} FAIL"); sys.exit(1 if FAIL else 0)
