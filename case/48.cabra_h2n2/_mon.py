import sys, h5py, numpy as np
r, s = sys.argv[1], int(sys.argv[2])
with h5py.File(f"{r}/cabra.h5") as m: xyz=m["MESH/COORD"][:].reshape(-1,3)
with h5py.File(f"{r}/res_{s}.h5") as h:
    V=h["VALUE"]; T=V["T"][:]; oh=V["Y4"][:]; TQ=V["cmc_TQmax"][:]; P=V["P"][:]; qd=V["chemQdot"][:]
N=len(T); x=xyz[:N,0]; y=xyz[:N,1]; D=0.00457
sel=oh>2e-4; ign=x[sel].min()/D if sel.any() else float('nan')
ax=np.where(np.abs(y)<1e-9)[0]; o=np.argsort(x[ax]); Tax=lambda zd: np.interp(zd*D,x[ax][o],T[ax][o])
print(f"{r.split('/')[-1]} step {s}: Tmax={T.max():.0f} K, ignx/d(OH)={ign:.1f}, TQmax={np.nanmax(TQ):.0f} K, P {P.min()/1e3:.1f}-{P.max()/1e3:.1f} kPa, Qdot {qd.min():.1e}..{qd.max():.1e}, axis T z/d14/20/26/30={Tax(14):.0f}/{Tax(20):.0f}/{Tax(26):.0f}/{Tax(30):.0f} (exp 620/1120/1410/1490), NaN={np.isnan(T).sum()}", flush=True)
