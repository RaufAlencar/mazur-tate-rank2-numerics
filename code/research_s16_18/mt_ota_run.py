import cypari2,sys,itertools; pari=cypari2.Pari(); pari.allocatemem(6*10**9)
pari('read("mt_full.gp")'); pari('read("mt_ota.gp")')
ainvs=sys.argv[1]; mmax=int(sys.argv[2]); gmin=int(sys.argv[3]) if len(sys.argv)>3 else 2
pari(f'E=ellinit({ainvs}); Ms=msfromell(E,1); G=ellsaturation(E,ellrank(E)[4],1000); N0=ellglobalred(E)[1]; Den=if(sign(E.disc)>0,1,2);')
print(ainvs,'gens',pari('G'),'Den',pari('Den'))
ms=[int(q) for q in pari.primes(400) if 7<=int(q)<=mmax and int(pari('N0'))%int(q)]
from math import gcd
pairs=[(a,b) for a,b in itertools.combinations(ms,2) if gcd((a-1)//2,(b-1)//2)>=gmin]
print(len(pairs),'pairs')
for m1,m2 in pairs:
    r=pari(f'otacomp(E,Ms,{m1},{m2},G,Den)')
    if int(r[0])<0: print(m1,m2,'FAIL',r); continue
    n1,n2,g,a20,a11,a02,M,D1,D12,D2,ap1,ap2,p1,p2=[r[i] for i in range(14)]
    n1,n2,g,a20,a11,a02,M,D1,D12,D2,ap1,ap2=map(int,(n1,n2,g,a20,a11,a02,M,D1,D12,D2,ap1,ap2))
    c21=int(p1[1]); c22=int(p2[1])
    # sanity: x1^2 coefficient should be (a_{m2}-2)*c2(m1) mod n1 ; x2^2 similarly
    s1=(a20-(ap2-2)*c21)%n1; s2=(a02-(ap1-2)*c22)%n2
    # new: 4 M^3 a11 ≡ -u D12 (mod g) for the prime-level constant u (test u in a small set)
    lhs=(4*pow(M,3,g)*a11)%g
    fits=[u for u in (-16,-8,-4,-2,-1,1,2,4,8,16) if (lhs+u*D12)%g==0]
    print(f'm1={m1} m2={m2} g={g}  x1^2 norm-check {"OK" if s1==0 else "BAD "+str(s1)}  x2^2 {"OK" if s2==0 else "BAD "+str(s2)}  | 4M^3 a11={lhs} D12={D12%g} mod {g}  u fits: {fits}',flush=True)
