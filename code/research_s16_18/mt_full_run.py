import cypari2,sys; pari=cypari2.Pari(); pari.allocatemem(8*10**9); pari('read("mt_full.gp")')
from math import gcd
ainvs=sys.argv[1]; mmax=int(sys.argv[2]); U=int(sys.argv[3]) if len(sys.argv)>3 else 8
pari(f'E=ellinit({ainvs}); G=ellsaturation(E,ellrank(E)[4],1000); P1=G[1]; P2=G[2]; N0=ellglobalred(E)[1];'); N0=int(pari('N0'))
ok=bad=triv=weak=0; fails=[]
for m in [int(q) for q in pari.primes(400) if 5<=int(q)<=mmax and N0%int(q)]:
    r=pari(f'mtfullu(E,{m},P1,P2,[-1])')
    if int(r[0])<0: print('m',m,'theta not in I^2 ?!',r); bad+=1; continue
    n,c2,N=int(r[0]),int(r[1]),int(r[2]); fourdet=int(r[3][0])
    lhs=(4*pow(N,3,n)*c2)%n; rhs=(U*fourdet)%n   # identity: 4 N^3 c2 = U * 4 det
    if lhs==0 and fourdet==0: triv+=1; continue
    if lhs==rhs:
        ok+=1
        if gcd(N,n)>1: weak+=1
    else: bad+=1; fails.append((m,n,lhs,rhs,N%n))
print(f'{ainvs}: identity N^3*c2 = {U}*det(<R_i,R_j>) mod (m-1)/2 :  OK {ok} (of which {weak} with gcd(N,n)>1), trivial(both 0) {triv}, FAIL {bad}')
if fails: print('  fails:',fails[:12])
