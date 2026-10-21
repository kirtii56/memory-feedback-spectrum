import math, sys, platform, datetime
import sympy as sp
print("Python", sys.version.split()[0], "| SymPy", sp.__version__, "|", platform.platform())
print("UTC", datetime.datetime.now(datetime.timezone.utc).isoformat())
rows=[]
def chk(case,q,exp,comp,tol):
    d=abs(comp-exp); rel=d/abs(exp) if exp else float('nan')
    st="PASS" if d<=tol else "FAIL"; rows.append((case,q,exp,comp,d,rel,tol,st))
# Type A
al=0.5; D,E=-1.0,0.0; trajA=[D]
for _ in range(4): D,E=2*al*E,D; trajA.append(D)
for t,(c,e) in enumerate(zip(trajA,[-1,0,-1,0,-1])): chk("A",f"D{t}",e,c,1e-12)
a=sp.Rational(1,2); mu=sp.symbols('mu')
JA=sp.Matrix([[0,2*a],[1,0]]); ev=sorted(JA.eigenvals().keys(),key=float)
chk("A","eig-",-1,float(ev[0]),1e-12); chk("A","eig+",1,float(ev[1]),1e-12)
# Type B
al=0.1; P,D,E=6.0,-1.0,0.0; trajB=[D]
for _ in range(4): P,D,E=P**2*(1-al**2*E**2),2*al*P*E,D; trajB.append(D)
for t,e in enumerate([-1,0,-7.2,0]): chk("B",f"D{t}",e,trajB[t],1e-9)
chk("B","D4 (order of magnitude, tol 1% rel)",-2370000.0,trajB[4],2370000*0.01)
x=sp.symbols('x'); Pp,Dd,Ee=sp.symbols('P D E'); a=sp.Rational(1,10)
F=sp.Matrix([Pp**2*(1-a**2*Ee**2),2*a*Pp*Ee,Dd]); J=F.jacobian([Pp,Dd,Ee]).subs({Pp:1,Dd:0,Ee:0})
evs=sorted([float(k) for k in J.eigenvals()])
for e,c in zip([-math.sqrt(0.2),math.sqrt(0.2),2.0],evs): chk("B","eig@(1,0,0)",e,c,1e-12)
# Type C
al=0.3; D,E=-1.0,0.0; tr=[D]
for _ in range(6): D,E=2*al*math.tanh(E),D; tr.append(D)
D,E=-1.0,0.0; trC=[]
# recompute properly (Run1 logic: D_{t+1}=2*alpha*tanh(E_t), E_{t+1}=D_t)
D,E=-1.0,0.0; seq=[D]
for _ in range(6): D,E=2*al*math.tanh(E),D; seq.append(D)
chk("C","D2",-0.4570,seq[2],5e-5); chk("C","D4",-0.2566,seq[4],5e-5); chk("C","D6",-0.1506,seq[6],5e-5)
# C-asym
al_=sp.Rational(3,10)
Ja=sp.Matrix([[1,0,al_],[0,0,0],[1,0,0]])
cp=sp.expand((Ja-mu*sp.eye(3)).det())
print("det(J-mu I) =",sp.factor(cp), "| expected -mu(mu^2-mu-alpha):", sp.simplify(cp-(-mu*(mu**2-mu-al_)))==0)
A=sp.Matrix([[1,0],[0,0]]); B=sp.Matrix([al_,0]); C=sp.Matrix([[1,0]])
Nmu=sp.simplify((C*(A-mu*sp.eye(2)).adjugate()*B)[0]); print("N(mu) =",Nmu)
chk("C-asym","N(mu)=-alpha*mu at mu=2",-0.6,float(Nmu.subs(mu,2)),1e-12)
lemma=sp.simplify(cp-((A-mu*sp.eye(2)).det()*(0-mu)-Nmu))==0; print("Lemma 3.1 identity:",lemma)
mp=max(float(r) for r in sp.solve(cp,mu)); chk("C-asym","mu+",1.2416,mp,5e-5)
# Type D
lam,tau=0.5,0.05; g=math.exp(-lam); D,M=-1.0,0.0; dd=[D]
for _ in range(5): D,M=2*tau*M,g*(M+D); dd.append(D)
mup=(g+math.sqrt(g*(g+8*tau)))/2; chk("D","mu+",0.694,mup,5e-4)
chk("D","tau_c",0.3244,(math.exp(lam)-1)/2,5e-5)
chk("D","D3/D2",g,dd[3]/dd[2],1e-12)
gs,ts=sp.symbols('g t',positive=True); JD=sp.Matrix([[0,2*ts],[gs,gs]])
print("Type D char poly:",sp.factor((JD-mu*sp.eye(2)).det()))
mn=min(float(r) for r in sp.solve(((JD-mu*sp.eye(2)).det()).subs({gs:g,ts:tau}),mu)); print("Type D mu-:",mn)
# Type D second-order recursion check D_{t+1}=g D_t + 2 tau g D_{t-1}
chk("D","2nd-order rec at t=4",dd[4],g*dd[3]+2*tau*g*dd[2],1e-14)
chk("D","2nd-order rec at t=5",dd[5],g*dd[4]+2*tau*g*dd[3],1e-14)
print("\n| Case | Quantity | Expected | Computed | AbsErr | RelErr | Tol | Status |")
for r in rows: print("| %s | %s | %.10g | %.10g | %.3e | %.3e | %.1e | %s |"%r)
nf=sum(r[-1]=="FAIL" for r in rows); print("\nTOTAL",len(rows),"checks;",nf,"FAIL")
print("Lemma 3.1 identity ok:",lemma)
sys.exit(1 if nf or not lemma else 0)
