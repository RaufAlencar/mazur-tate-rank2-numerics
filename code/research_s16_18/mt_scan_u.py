import cypari2,sys; pari=cypari2.Pari(); pari.allocatemem(6*10**9); pari('read("mt_full.gp")')
from math import gcd
from collections import Counter
ainvs=sys.argv[1]; mmax=int(sys.argv[2])
pari(f'E=ellinit({ainvs}); G=ellsaturation(E,ellrank(E)[4],1000); P1=G[1]; P2=G[2]; N0=ellglobalred(E)[1];'); N0=int(pari('N0'))
print(ainvs,'N=',N0,'tamagawa',pari('vector(#ellglobalred(E)[5],i,ellglobalred(E)[5][i][4])'),'tors',pari('elltors(E)[1]'),'omega',pari('E.omega[1]'),'disc sign',pari('sign(E.disc)'),'gens',pari('[P1,P2]'))
hits=Counter(); tested=0
cands=[(s*num,dn) for s in (1,-1) for num in (1,2,4,8,16,32,3,6,12,24,9) for dn in (1,2,4,8,3)]
for m in [int(q) for q in pari.primes(300) if 23<=int(q)<=mmax and N0%int(q)]:
    r=pari(f'mtfullu(E,{m},P1,P2,[-1])')
    if int(r[0])<0: continue
    n,c2,N=int(r[0]),int(r[1]),int(r[2]); fourdet=int(r[3][0])
    lhs=(4*pow(N,3,n)*c2)%n
    if lhs==0 and fourdet==0: continue
    tested+=1
    for num,dn in cands:
        if gcd(dn,n)==1 and (lhs*dn-num*fourdet)%n==0: hits[f'{num}/{dn}' if dn!=1 else str(num)]+=1
print('tested',tested,' best constants u (N^3 c2 = u*det):',hits.most_common(6))
