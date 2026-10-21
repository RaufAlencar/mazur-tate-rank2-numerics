import cypari2, sys, cmath, math
exec(open('heegner.py').read().split("label=sys.argv")[0])
pari.default('realprecision',60)
label=sys.argv[1]; E=pari.ellinit(pari(sys.argv[2])); N=int(pari.ellglobalred(E)[0])
an=pari.ellan(E,600000)
for D in heegner_discs(N,int(sys.argv[3])):
    forms=heegner_forms(N,D); h=len(forms)
    if h<2: continue
    zs=heegner_point(E,N,D,an)
    print(f"--- D={D} h={h}  K=Q(sqrt{D})")
    for i,z in enumerate(zs[:3]):
        P=pari.ellztopoint(E,pari(z.real)+pari('I')*pari(z.imag))
        x=P[0]
        # recognise x as algebraic of degree <= 2h over Q
        pol=pari.algdep(x,2*h)
        ht=pari.ellheight(E,P) if False else None
        # numerical canonical height via complex: use ellheight on algebraic point if recognised
        rts=pari.polroots(pol)
        print(f"  z_{i}: x≈{complex(x).real:+.6f}{complex(x).imag:+.6f}i  minpoly deg {pari.poldegree(pol)}: {pol}")
    # heights: canonical height of individual point over H via sum of local heights ~ use archimedean formula check: height of z_0 relative to that of torsion (0)
    # cheap test: multiply z_0 by m for m=1..12 and see if it ever hits lattice
    z0=zs[0]; import numpy as np
    tor=[m for m in range(1,25) if np.allclose([((m*z0).real), ((m*z0).imag)], np.round(np.linalg.solve(M,[(m*z0).real,(m*z0).imag]))@M.T if False else [0,0]) ]
    c=np.linalg.solve(M,[z0.real,z0.imag]); 
    hits=[m for m in range(1,61) if np.allclose((m*c)%1,0,atol=1e-6) or np.allclose((m*c)%1,1,atol=1e-6)]
    print(f"  z_0 lattice coords ({c[0]:.6f},{c[1]:.6f}); m·z_0 ∈ Λ for m in {hits[:5] or 'none ≤ 60  → ordem infinita (ou torção > 60, impossível por Mazur sobre K-gráu)'}")
