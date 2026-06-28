#!/usr/bin/env python3
"""入口分布 (壁法則 BL/チャネル) プロファイル CSV を生成する。

forge の inlet 分布機能 (bcondConfig の `ints: {inletProfile: 1}`) 用テーブルを作る。
出力: `inlet_profile_<physID>.csv`。1 行目ヘッダ = 補間方向 + 量名。
本スクリプトは 1D-y 壁法則を生成 (チャネル: 上下壁、または片側 BL)。
Reichardt 合成則 u+ = (1/k)ln(1+k y+) + 7.8[1 - exp(-y+/11) - (y+/11)exp(-y+/3)]。

使い方:
  python3 gen_inlet_walllaw.py --physID 1 --ylo 1.0 --yhi 3.0 --Uc 49.0 --nu 8.5e-4 \
      [--mode channel|bl] [--ny 101] [--Uy 0 --Uz 0]
  nu = visc/ro (動粘性)。例: visc=0.001, ro=1.18 -> nu=8.5e-4。
  --mode channel: 上下壁 (ylo,yhi) の両方から壁法則、中央で Uc。
  --mode bl:      ylo 壁のみから BL、delta(=yhi-ylo) 厚で Uc。
"""
import argparse, numpy as np

KAPPA = 0.41
def reichardt(yp):
    yp = np.maximum(yp, 1e-12)
    return (1.0/KAPPA)*np.log(1.0+KAPPA*yp) + 7.8*(1.0 - np.exp(-yp/11.0) - (yp/11.0)*np.exp(-yp/3.0))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--physID', type=int, required=True)
    ap.add_argument('--ylo', type=float, required=True, help='下壁 y')
    ap.add_argument('--yhi', type=float, required=True, help='上壁 y (channel) / BL 外縁 (bl)')
    ap.add_argument('--Uc', type=float, required=True, help='中央/外縁 速度 [m/s]')
    ap.add_argument('--nu', type=float, required=True, help='動粘性 visc/ro [m2/s]')
    ap.add_argument('--mode', choices=['channel','bl'], default='channel')
    ap.add_argument('--ny', type=int, default=201)
    ap.add_argument('--Uy', type=float, default=0.0)
    ap.add_argument('--Uz', type=float, default=0.0)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    y = np.linspace(a.ylo, a.yhi, a.ny)
    if a.mode == 'channel':
        h = 0.5*(a.yhi - a.ylo)          # 半高
        d = np.minimum(y - a.ylo, a.yhi - y)  # 最近壁距離
        d = np.maximum(d, 0.0)
        dref = h
    else:  # bl
        d = np.maximum(y - a.ylo, 0.0)
        dref = a.yhi - a.ylo             # BL 厚

    # u_tau を「壁法則の最大 (中央/外縁) = Uc」になるよう求める
    # u(d)=u_tau*reichardt(d*u_tau/nu); d=dref で Uc。1 変数の単調方程式を二分法で解く。
    def umax(utau):
        return utau*reichardt(dref*utau/a.nu)
    lo, hi = 1e-6, a.Uc  # u_tau < Uc
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if umax(mid) < a.Uc: lo = mid
        else: hi = mid
    utau = 0.5*(lo+hi)

    yp = d*utau/a.nu
    Ux = utau*reichardt(yp)
    Ux = np.minimum(Ux, a.Uc)            # 中央でクランプ
    Ux[0] = 0.0; Ux[-1] = 0.0 if a.mode=='channel' else Ux[-1]  # 壁で 0

    out = a.out or f'inlet_profile_{a.physID}.csv'
    with open(out, 'w') as f:
        f.write('y Ux Uy Uz\n')
        for yi, ui in zip(y, Ux):
            f.write(f'{yi:.6f} {ui:.6f} {a.Uy:.6f} {a.Uz:.6f}\n')
    print(f'wrote {out}: mode={a.mode}, u_tau={utau:.4f}, Uc={a.Uc}, y+wall(1cell)~{(y[1]-y[0])*utau/a.nu:.2f}')
    print(f'  Ux range {Ux.min():.3f}..{Ux.max():.3f} ({a.ny} rows, header "y Ux Uy Uz")')

if __name__ == '__main__':
    main()
