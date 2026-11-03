# Conferencia independente (registro): ker red_5 ∩ ker red_17 em E(Q)/3, 11642a1
import cypari2
pari = cypari2.Pari()
pari('E=ellinit([1,-1,0,-16,28]); P=[[2,-2],[3,-2],[4,2]]')
pari('red(Q,l)=if(Q==[0]||valuation(denominator(Q[1]),l)>0,[0],[Mod(numerator(Q[1]),l)/Mod(denominator(Q[1]),l),Mod(numerator(Q[2]),l)/Mod(denominator(Q[2]),l)])')
pari('inker(v,l)=my(El=ellinit(E[1..5],l),N=ellcard(El),c=ellgroup(El),Q=[0]);for(i=1,3,Q=elladd(E,Q,ellmul(E,P[i],v[i])));if(N%3,1,#c>1&&c[2]%3==0,-1,ellmul(El,red(Q,l),N/3)==[0])')
for l in (5, 17):
    print("l", l, "E(F_l) =", pari(f'ellgroup(ellinit(E[1..5],{l}))'))
    ker = [v for v in [(a,b,c) for a in range(3) for b in range(3) for c in range(3)] if int(pari(f'inker({list(v)},{l})')) == 1]
    print("  ker red mod 3:", ker)
k = [v for v in [(a,b,c) for a in range(3) for b in range(3) for c in range(3)] if int(pari(f'inker({list(v)},5)'))==1 and int(pari(f'inker({list(v)},17)'))==1]
print("intersecao:", k)
print("kronecker(-7,l):", int(pari('kronecker(-7,5)')), int(pari('kronecker(-7,17)')), " a_l:", pari('[ellap(E,5),ellap(E,17)]'))
