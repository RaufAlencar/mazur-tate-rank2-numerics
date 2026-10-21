import cypari2,sys; pari=cypari2.Pari(); pari.allocatemem(8*10**9); pari('read("mt_full.gp")'); pari('read("mt_rank.gp")')
from math import gcd
from collections import Counter
ainvs=sys.argv[1]; mmax=int(sys.argv[2]); r=int(sys.argv[3])
pari(f'E=ellinit({ainvs}); G=ellsaturation(E,ellrank(E)[4],1000); N0=ellglobalred(E)[1]; cT=prod(i=1,#ellglobalred(E)[5],ellglobalred(E)[5][i][4]);')
print(ainvs,'rank gens',pari('G'),'Tamagawa prod',pari('cT'),'disc sign',pari('sign(E.disc)'))
assert int(pari('#G'))==r
N0=int(pari('N0')); hits=Counter(); tested=0; notin=0
cands=[(s*num,dn) for s in (1,-1) for num in (1,2,4,8,16,32,64,128,256) for dn in (1,2,4,8,16)]
for m in [int(q) for q in pari.primes(200) if 23<=int(q)<=mmax and N0%int(q)]:
    tc=pari(f'thetaclass(E,{m},{r})')
    if int(tc[0])<0: notin+=1; print('  m',m,'theta NOT in I^',r,'?',tc); continue
    n=int(tc[0]); cr=int(tc[1]); Mm=int(pari(f'cT*ellcard(E,{m})'))
    d=int(pari(f'pairdet(E,{m},G,{Mm},{n},znprimroot({m}))'))
    lhs=(2**r*pow(Mm,2*r-1,n)*cr)%n     # 2^r M^(2r-1) c_r  vs u * 2^r det h = u*d
    if lhs==0 and d==0: continue
    tested+=1
    for num,dn in cands:
        if gcd(dn,n)==1 and (lhs*dn-num*d)%n==0: hits[f'{num}/{dn}' if dn!=1 else str(num)]+=1
print(f'rank {r}: tested {tested}, theta not in I^r: {notin}; best u (2^r M^(2r-1) c_r = u * 2^r det h):',hits.most_common(6))
