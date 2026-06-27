# bump 三者ライン比較 (y=0.25): SU2 vol_solution.vtu vs forge node/cell res_*.h5。
# forge 側 run は temp (破棄済)。再現するには run_dual_loM_{cell,node}_m2 を複製し
#   solver: "HLLE" にして cell=40000 / node=150000 step 回し、下の res パスを差し替える。
# SU2: gmsh -2 mesh/bump2d.geo -format su2 -o mesh/bump2d.su2 ; SU2_CFD bump_euler.cfg
import numpy as np, h5py, re
from scipy.interpolate import griddata
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

def load_su2(path):
    raw=open(path,'rb').read()
    hdr=raw[:raw.find(b'<AppendedData')].decode('latin1')
    arrs={}
    for m in re.finditer(r'Name="([^"]*)"\s*NumberOfComponents=\s*"(\d+)"\s*offset="(\d+)"',hdr):
        arrs[m.group(1)]=(int(m.group(2)),int(m.group(3)))
    # points array has empty Name
    mpts=re.search(r'<DataArray type="Float32" Name=""\s*NumberOfComponents=\s*"(\d+)"\s*offset="(\d+)"',hdr)
    npts=int(re.search(r'NumberOfPoints="(\d+)"',hdr).group(1))
    i=raw.find(b'<AppendedData'); i=raw.find(b'_',i)+1
    def rd(off,ncomp,n):
        ln=int.from_bytes(raw[i+off:i+off+8],'little')
        a=np.frombuffer(raw[i+off+8:i+off+8+ln],dtype=np.float32)
        return a.reshape(n,ncomp) if ncomp>1 else a
    pts=rd(int(mpts.group(2)),3,npts)
    P=rd(arrs['Pressure'][1],1,npts); T=rd(arrs['Temperature'][1],1,npts)
    V=rd(arrs['Velocity'][1],3,npts)
    return pts[:,0],pts[:,1],P,T,np.sqrt(V[:,0]**2+V[:,1]**2)

def load_forge(path):
    f=h5py.File(path,'r'); coord=f['MESH/COORD'][:].reshape(-1,3)
    P=f['VALUE/P'][:];T=f['VALUE/T'][:];Ux=f['VALUE/Ux'][:];Uy=f['VALUE/Uy'][:]
    nval=len(P)
    if nval==coord.shape[0]:
        x,y=coord[:,0],coord[:,1]                      # node-centered: 値=ノード
    else:
        conn=f['MESH/CONNE'][:]                         # [5,n0,n1,n2,n3]*ncell (5=XDMF quad)
        cells=conn.reshape(-1,5)[:,1:]                  # 4 node ids/cell
        cc=coord[cells].mean(axis=1)                    # セル重心
        x,y=cc[:,0],cc[:,1]
        assert len(x)==nval, (len(x),nval)
    return x,y,P,T,np.sqrt(Ux**2+Uy**2)

src={'SU2':load_su2('vol_solution.vtu'),
     'forge-node':load_forge('../run_tmp_hlle_node/res_150000.h5'),
     'forge-cell':load_forge('../run_tmp_hlle_cell/res_40000.h5')}

xline=np.linspace(0.02,2.98,200); yL=0.25
def sample(d):
    x,y,P,T,U=d
    pts=np.column_stack([x,y]); q=np.column_stack([xline,np.full_like(xline,yL)])
    return (griddata(pts,P,q,method='linear'),griddata(pts,T,q,method='linear'),griddata(pts,U,q,method='linear'))
S={k:sample(v) for k,v in src.items()}

fig,ax=plt.subplots(1,3,figsize=(15,4.2))
labels=['Pressure [Pa]','Temperature [K]','|U| [m/s]']
col={'SU2':'k','forge-node':'tab:blue','forge-cell':'tab:red'}
sty={'SU2':'-','forge-node':'--','forge-cell':':'}
for j,t in enumerate(labels):
    for k in src: ax[j].plot(xline,S[k][j],sty[k],color=col[k],lw=1.8,label=k)
    ax[j].set_xlabel('x'); ax[j].set_title(t+f'  (y={yL})'); ax[j].grid(alpha=.3); ax[j].legend()
plt.tight_layout(); plt.savefig('bump_line_compare.png',dpi=120); print('saved bump_line_compare.png')

print('\n=== y=0.25 ライン: SU2 基準の最大相対差 (%) ===')
for j,t in enumerate(['P','T','|U|']):
    su=S['SU2'][j]; 
    for k in ['forge-node','forge-cell']:
        d=S[k][j]; mask=np.isfinite(su)&np.isfinite(d)&(np.abs(su)>1e-9)
        rel=np.abs(d[mask]-su[mask])/np.abs(su[mask])
        print(f'  {t:3s} {k:11s}: max {100*rel.max():.2f}%  mean {100*rel.mean():.2f}%')
print('\n=== 代表値 (x=1.5, bump頂上直上 y=0.25) ===')
i15=np.argmin(np.abs(xline-1.5))
for k in src: print(f'  {k:11s}: P={S[k][0][i15]:.1f}  T={S[k][1][i15]:.2f}  |U|={S[k][2][i15]:.2f}')
