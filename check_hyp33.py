"""Explicit trajectory pairs for Hypothesis 3.3' (equal x_t, different memory m_t, different F outputs).
Admissible trajectories = arbitrary initial state AND arbitrary initial memory (pre-history)."""
import math
def A(a,b,Em,al=0.5):  # one step, memory m_a=m, m_b=-m ; returns (a,b)
    return (a+b+al*Em, a+b-al*Em)
ok=True
def report(name,x1,x1p,m,mp,F1,F1p):
    global ok
    same=all(abs(u-v)<1e-12 for u,v in zip(x1,x1p)); dm=abs(m-mp)>1e-9; dF=any(abs(u-v)>1e-9 for u,v in zip(F1,F1p))
    st="PASS" if (same and dm and dF) else "FAIL"; ok&=(st=="PASS")
    print(f"{name}: same x_t={same}, m differs={dm}, F-output differs={dF} -> {st}")
# Type A, alpha=.5. State (a,b), memory E_t=D_{t-1}. Vary D_0 keeping S_0,E_0.
def stepA(a,b,E,al=0.5): return (a+b+al*E, a+b-al*E)
def pair_AC(f,al,name):
    a0,b0,E0=2.0,3.0,0.7; S=a0+b0
    for D0 in (-1.0, 0.4):
        pass
    res=[]
    for D0 in (-1.0,0.4):
        a=(S+D0)/2; b=(S-D0)/2
        x1=(a+b+al*f(E0), a+b-al*f(E0))      # a_1,b_1 (m_a=E_0)
        # memory at t=1 is D_0; F output at t=1:
        F=(x1[0]+x1[1]+al*f(D0), x1[0]+x1[1]-al*f(D0))
        res.append((x1,D0,F))
    report(name,res[0][0],res[1][0],res[0][1],res[1][1],res[0][2],res[1][2])
pair_AC(lambda e:e,0.5,"Type A")
pair_AC(math.tanh,0.3,"Type C")
# Type B: a'=ab(1+al m), b'=ab(1+al(-m)) with m=E
def pair_B(al=0.1):
    out=[]
    for D0 in (-1.0,2.0):
        P0=6.0; a=(D0+math.sqrt(D0*D0+4*P0))/2; b=a-D0   # a*b=P0
        E0=0.5
        a1=a*b*(1+al*E0); b1=a*b*(1-al*E0)
        F=(a1*b1*(1+al*D0), a1*b1*(1-al*D0))
        out.append(((a1,b1),D0,F))
    report("Type B",out[0][0],out[1][0],out[0][1],out[1][1],out[0][2],out[1][2])
pair_B()
# Type C-asym: x=(a,b), m=a_{t-1}; a'=a+b^2+al*m, b'=a*b
def pair_Casym(al=0.3):
    a0,b0,c0=1.0,2.0,0.5
    a0p=2.0; b0p=a0*b0/a0p
    x1=(a0+b0**2+al*c0, a0*b0)
    c0p=(x1[0]-a0p-b0p**2)/al
    x1p=(a0p+b0p**2+al*c0p, a0p*b0p)
    F=(x1[0]+x1[1]**2+al*a0, x1[0]*x1[1]); Fp=(x1p[0]+x1p[1]**2+al*a0p, x1p[0]*x1p[1])
    report("Type C-asym",x1,x1p,a0,a0p,F,Fp)
pair_Casym()
# Type D: x=(S,D), memory M. Same x_2, different M_2.
def pair_D(lam=0.5,tau=0.05):
    g=math.exp(-lam); res=[]
    for (M0,D0) in ((0.3,-1.0),(0.9,-1.0+0.3-0.9)):  # chosen so M_1=g(M0+D0) equal
        D1=2*tau*M0; M1=g*(M0+D0); D2=2*tau*M1; M2=g*(M1+D1)
        res.append(((D2,),M2,(2*tau*M2,)))   # S doubles identically; compare D-coordinate
    report("Type D",res[0][0],res[1][0],res[0][1],res[1][1],res[0][2],res[1][2])
pair_D()
print("ALL PASS" if ok else "SOME FAIL")
