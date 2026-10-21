import cypari2, sys, cmath, math, numpy as np
pari=cypari2.Pari(); pari.allocatemem(2*10**9); pari.default('realprecision',40)
def heegner_discs(N,k=8,start=-3):
    out=[]; D=start
    while len(out)<k:
        if pari.isfundamental(D) and all(int(pari.kronecker(D,int(p)))==1 for p in pari.factor(N)[0]):
            out.append(int(D))
        D-=1
    return out
def heegner_forms(N,D):
    h=int(pari.qfbclassno(D)); betas=[b for b in range(2*N) if (b*b-D)%(4*N)==0]
    beta=betas[0]
    seen={}; k=1
    while len(seen)<h and k<500:
        A=N*k
        B=beta
        while B<=A:
            for Bs in (B,B-2*N) if False else (B,):
                pass
            B+=2*N
        for B in range(-A+1,A+1):
            if (B-beta)%(2*N): continue
            if (B*B-D)%(4*A): continue
            C=(B*B-D)//(4*A)
            if math.gcd(math.gcd(A,abs(B)),C)!=1: continue
            cls=str(pari.qfbred(pari.Qfb(A,B,C)))
            if cls not in seen: seen[cls]=(A,B,C)
        k+=1
    assert len(seen)==h,(len(seen),h)
    return list(seen.values())
def heegner_point(E,N,D,an):
    pts=[]
    for (A,B,C) in heegner_forms(N,D):
        tau=(-B+cmath.sqrt(D))/(2*A); q=cmath.exp(2j*math.pi*tau)
        nterms=min(len(an),int(40/(2*math.pi*tau.imag))+10)
        z=sum(complex(an[n-1])/n*q**n for n in range(1,nterms+1))
        pts.append(z)
    return pts
label=sys.argv[1]; E=pari.ellinit(pari(sys.argv[2])); N=int(pari.ellglobalred(E)[0])
tors=pari.elltors(E); print(label,'N=',N,'torsion',tors)
w1,w2=complex(E.omega()[0]),complex(E.omega()[1]); M=np.array([[w1.real,w2.real],[w1.imag,w2.imag]])
an=pari.ellan(E,400000)
for D in heegner_discs(N,int(sys.argv[3]) if len(sys.argv)>3 else 6):
    zs=heegner_point(E,N,D,an); z=sum(zs)
    c=np.linalg.solve(M,[z.real,z.imag])
    print(f"D={D:5d} h={len(zs):2d}  z ≡ {c[0]%1:.5f}·w1 + {c[1]%1:.5f}·w2  (mod Λ)")
