# Check at the Bley–Macías Castillo levels: det H mod 3 and the full identity mod n.
# usage: python3 bmc_levels.py            (needs mt_full.gp, mt_rank.gp in cwd)
import cypari2; pari=cypari2.Pari(); pari.allocatemem(3*10**9); pari('read("mt_full.gp")'); pari('read("mt_rank.gp")')
for a,ms,lab in [('[0,1,1,-2,0]',[7,13,43],'389a1'),('[1,-1,0,-4,4]',[19],'446d1'),('[1,0,0,0,1]',[7,13,31,37],'433a1')]:
    pari(f'E=ellinit({a}); G=ellsaturation(E,ellrank(E)[4],1000); cT=prod(i=1,#ellglobalred(E)[5],ellglobalred(E)[5][i][4]);')
    print(lab,'gens',pari('G'),'tam',pari('cT'))
    for m in ms:
        tc=pari(f'thetaclass(E,{m},2)'); n=(m-1)//2; Mm=int(pari(f'cT*ellcard(E,{m})'))
        d=int(pari(f'pairdet(E,{m},G,{Mm},{n},znprimroot({m}))'))
        lhs=(4*Mm**3*int(tc[1]))%n
        print(f'  m={m} n={n} c2={int(tc[1])} detH mod n={d} detH mod 3={d%3}  4M^3c2={lhs}  8detH={(8*d)%n}  4detH={(4*d)%n}')
