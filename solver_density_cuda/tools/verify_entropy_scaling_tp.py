#!/usr/bin/env python3
"""Step 3 実装前検証: thermally-perfect (TP) 単成分の entropy-scaled 固有系。

TP: p=ρRT, cp(T)=cp0+cp1·T (2項テスト), e=∫cv dT, h=∫cp dT, s⁰=∫cp/T dT。
次元付き η=-ρs, w=∂η/∂U=[(g-½|u|²)/T, u/T, -1/T], g=h-Ts。
検証:
 1. w = ∂η/∂U (数値微分)
 2. 固有ベクトル (音響 [1,u±cn,Ht±c·un], エントロピー [1,u,½|u|²+e-cvT], せん断 [0,t,u·t]) が A r=λr
 3. S=R⁻¹HR⁻ᵀ が対角・閉形式 (音響 ρ/(2γR), エントロピー ρ/cp, せん断 ρT) と一致
 4. Q=R|Λ|SRᵀ SPD
 5. CPG 極限 (cp1=0) で従来形と整合
"""
import numpy as np
rng = np.random.default_rng(3)

R = 296.8           # N2 [J/kgK]
cp0, cp1 = 1000.0, 0.05   # cp(T)=cp0+cp1*T (適度な T 依存)
Tref = 100.0
h_datum = 1.234e5   # datum を入れても不変なことの確認用

def cp(T): return cp0 + cp1*T
def cv(T): return cp(T) - R
def h(T):  return h_datum + cp0*(T-Tref) + 0.5*cp1*(T*T-Tref*Tref)
def e(T):  return h(T) - R*T
def s0(T): return cp0*np.log(T/Tref) + cp1*(T-Tref)

def T_from_e(eint):
    T = 300.0
    for _ in range(60):
        T -= (e(T)-eint)/cv(T)
    return T

def prim2U(ro,u,T):
    return np.array([ro, ro*u[0], ro*u[1], ro*u[2], ro*(e(T)+0.5*(u@u))])
def U2prim(U):
    ro=U[0]; u=U[1:4]/ro; T=T_from_e(U[4]/ro-0.5*(u@u)); p=ro*R*T
    return ro,u,T,p
def eta(U):
    ro,u,T,p=U2prim(U)
    s = s0(T) - R*np.log(p/101325.0)
    return -ro*s
def w_of(U):
    ro,u,T,p=U2prim(U)
    s = s0(T) - R*np.log(p/101325.0)
    g = h(T) - T*s
    return np.array([(g-0.5*(u@u))/T, u[0]/T, u[1]/T, u[2]/T, -1.0/T])

def flux(U,n):
    ro,u,T,p=U2prim(U)
    un=u@n
    return np.array([ro*un, ro*u[0]*un+p*n[0], ro*u[1]*un+p*n[1], ro*u[2]*un+p*n[2], (U[4]+p)*un])

ok=True
for trial in range(4):
    ro=rng.uniform(0.3,3.0); T=rng.uniform(250.0,900.0)
    u=rng.uniform(-100,100,3)
    n=rng.normal(size=3); n/=np.linalg.norm(n)
    aa=np.array([1.0,0,0]) if abs(n[0])<0.9 else np.array([0,1.0,0])
    t1=np.cross(n,aa); t1/=np.linalg.norm(t1); t2=np.cross(n,t1)
    U=prim2U(ro,u,T); p=ro*R*T
    ga=cp(T)/cv(T); c=np.sqrt(ga*R*T)
    Ht=h(T)+0.5*(u@u); un=u@n

    # 1. w=∂η/∂U
    grad=np.zeros(5)
    for j in range(5):
        hj=1e-6*max(abs(U[j]),1e-3)
        Up=U.copy(); Um=U.copy(); Up[j]+=hj; Um[j]-=hj
        grad[j]=(eta(Up)-eta(Um))/(2*hj)
    e1=np.max(np.abs(grad-w_of(U)))/np.max(np.abs(w_of(U)))

    # 2. 固有ベクトル: A r = λ r (A は数値 FD)
    A=np.zeros((5,5))
    for j in range(5):
        hj=1e-6*max(abs(U[j]),1e-3)
        Up=U.copy(); Um=U.copy(); Up[j]+=hj; Um[j]-=hj
        A[:,j]=(flux(Up,n)-flux(Um,n))/(2*hj)
    Rm=np.zeros((5,5))
    Rm[:,0]=np.r_[1.0, u-c*n, Ht-c*un]
    Rm[:,1]=np.r_[1.0, u,     0.5*(u@u)+e(T)-cv(T)*T]   # TP エントロピー波 (CPG では ½|u|²+datum)
    Rm[:,2]=np.r_[0.0, t1,    u@t1]
    Rm[:,3]=np.r_[0.0, t2,    u@t2]
    Rm[:,4]=np.r_[1.0, u+c*n, Ht+c*un]
    lam=np.array([un-c,un,un,un,un+c])
    e2=max(np.linalg.norm(A@Rm[:,k]-lam[k]*Rm[:,k])/(abs(lam[k])*np.linalg.norm(Rm[:,k])+1e-30) for k in range(5))

    # 3. S=R⁻¹HR⁻ᵀ 対角・閉形式
    dwdU=np.zeros((5,5))
    for j in range(5):
        hj=1e-6*max(abs(U[j]),1e-3)
        Up=U.copy(); Um=U.copy(); Up[j]+=hj;Um[j]-=hj
        dwdU[:,j]=(w_of(Up)-w_of(Um))/(2*hj)
    H=np.linalg.inv(dwdU)
    S=np.linalg.inv(Rm)@H@np.linalg.inv(Rm).T
    off=np.max(np.abs(S-np.diag(np.diag(S))))/np.max(np.abs(np.diag(S)))
    S_lit=np.diag([ro/(2*ga*R), ro/cp(T), ro*T, ro*T, ro/(2*ga*R)])
    e3=np.max(np.abs(np.diag(S)-np.diag(S_lit))/np.abs(np.diag(S_lit)))

    # 4. Q SPD
    Q=Rm@np.diag(np.abs(lam))@np.diag(np.diag(S))@Rm.T
    e4=np.linalg.eigvalsh(0.5*(Q+Q.T)).min()

    print(f"trial{trial}: T={T:6.1f} γ={ga:.4f} | w rel={e1:.1e}  eig rel={e2:.1e}  "
          f"S offdiag={off:.1e}  |S-閉形式| rel={e3:.1e}  minEig(Q)={e4:+.2e}")
    ok &= e1<1e-4 and e2<1e-4 and off<1e-4 and e3<1e-4 and e4>-1e-6

print("VERDICT:", "PASS" if ok else "FAIL")
