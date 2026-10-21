import cypari2, math, sys
pari=cypari2.Pari(); pari.allocatemem(4*10**9); PREC=int(sys.argv[1]) if len(sys.argv)>1 else 60
pari.default('realprecision',PREC)
N=389; D=-20; l=17; Dl=l*l*D
E=pari.ellinit([0,1,1,-2,0])
h=int(pari.qfbclassno(Dl)); print('h(O_17)=',h)
beta=[b for b in range(2*N) if (b*b-Dl)%(4*N)==0][0]
# enumerate Heegner forms, one per class, also record best Im tau representative
cls={}; k=1
while len(cls)<h and k<3000:
    A=N*k
    for B in range(-A+1,A+1):
        if (B-beta)%(2*N) or (B*B-Dl)%(4*A): continue
        C=(B*B-Dl)//(4*A)
        if math.gcd(math.gcd(A,abs(B)),C)!=1: continue
        key=str(pari.qfbred(pari.Qfb(A,B,C)))
        if key not in cls or A<cls[key][0]: cls[key]=(A,B,C)
    k+=1
assert len(cls)==h, len(cls)
keys=sorted(cls); print('max A used:',max(v[0] for v in cls.values()))
# group structure via composition
def comp(a,b): return str(pari.qfbred(pari.qfbcomp(pari.Qfb(*cls[a]) if False else pari(a),pari(b))))
ident=[k_ for k_ in keys if str(pari.qfbred(pari(k_)))==str(pari.qfbred(pari.Qfb(1,beta%2,(beta%2-Dl)//4)))]
# identity: principal form reduced
idf=str(pari.qfbred(pari.Qfb(1, Dl%2, (Dl%2 - Dl)//4))); print('identity class',idf, idf in cls)
def order(g):
    x=g;n=1
    while x!=idf: x=comp(x,g); n+=1
    return n
# genus character (-4/A) with A odd: pick representative with odd A coprime to 2
def genus(k_):
    f=pari(k_); A,B,C=int(f[0]),int(f[1]),int(f[2])
    for v in (A,C,A+B+C):
        if v%2: return int(pari.kronecker(-4,v))
ker=[k_ for k_ in keys if genus(k_)==1]; print('kernel of Pic(O_17)->Pic(O_K): size',len(ker))
gens=[g for g in ker if order(g)==18]; g=gens[0]; print('generator of G=Z/18 found; #elements of order 18:',len(gens))
G=[idf]; 
for i in range(17): G.append(comp(G[-1],g))
assert len(set(G))==18
# Heegner points numerically
an=pari.ellan(E,2000000)
def z_of(k_):
    A,B,C=cls[k_]; tau=(-B+pari.sqrt(pari(Dl)))/(2*A); y=float(pari.imag(tau))
    nterms=int(PREC*math.log(10)/(2*math.pi*y))+20
    if nterms>len(an): raise Exception(f'need {nterms} terms (A={A})')
    q=pari.exp(2*pari.Pi()*pari('I')*tau)
    s=pari(0); qn=pari(1)
    for n in range(1,nterms+1):
        qn*=q; s+=an[n-1]*qn/n
    return s
Z={k_:z_of(k_) for k_ in keys}; print('36 points computed')
w1,w2=E.omega()[0],E.omega()[1]
def red(z):
    import numpy as np
    M=np.array([[float(pari.real(w1)),float(pari.real(w2))],[float(pari.imag(w1)),float(pari.imag(w2))]])
    c=np.linalg.solve(M,[float(pari.real(z)),float(pari.imag(z))]); return c
# sanity: trace to K of y_17 should equal a_17 * y_1 - ... (Heegner relation Tr y_l = a_l y_1 for l inert); y_1=P_K=0 => trace in torsion
tot=sum(Z.values()); print('trace of all 36 (should be ~ lattice point since a_17*P_K=0):',red(tot))
# Kolyvagin derivative for each coset h: w_h = sum_i i * z_{h g^i}
W={}
for hcls in keys:
    s=pari(0)
    for i in range(18): s+= i*Z[comp(hcls,G[i])]
    W[hcls]=s
w=W[idf]; print('D_17 y_17 lattice coords:',red(w))
# predicted: kappa_17 in line <P1-P2>: check whether w - P in 3E(K[17]) ... first step: is w itself close to a lattice point / E(Q) point?
P1=pari.ellztopoint(E,w); print('D_17 y_17 as complex point x=',pari.real(P1[0]),'+',pari.imag(P1[0]),'i')
import pickle; pickle.dump({'keys':keys,'cls':cls,'Z':{k:(str(pari.real(v)),str(pari.imag(v))) for k,v in Z.items()},'G':G,'W':{k:(str(pari.real(v)),str(pari.imag(v))) for k,v in W.items()}},open('koly17.pkl','wb'))
