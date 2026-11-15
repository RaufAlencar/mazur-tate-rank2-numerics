# Tests for mazur_tate_pairing.gp.
#   T1  the worked example of the note's Appendix: 389a1, m = 23, g = 5  ->  H = (9, 7; 7, 9) mod 11
#   T2  mtmatrix == the half-pairing recipe of the independent implementation (independent/mt_indep.gp),
#       entry by entry, on the ten curves of the note's table, 7 <= m <= 100
#   T3  mtpair == Mazur-Tate 1987 / Portillo-Bobadilla 2019 eq. (18) up to sign, same curves and levels
#   T4  error handling: non-minimal model, point outside the formal group, m dividing N
# Run (WSL): ~/bsdvenv/bin/python test_mazur_tate_pairing.py
import cypari2, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
pari('read("mazur_tate_pairing.gp")'); pari('read("independent/mt_indep.gp")')
pari('dd(T)=sqrtint(denominator(T[1]))')
pari('gPB(E,R,S,Pp,m)={my(a=dd(elladd(E,Pp,R)),b=dd(elladd(E,Pp,S)),c=dd(Pp),e=dd(elladd(E,Pp,elladd(E,R,S))));Mod(a*b,m)/Mod(c*e,m)}')
CURVES = [("389a1", "[0,1,1,-2,0]", "[[-1,1],[0,0]]"), ("709a1", "[0,-1,1,-2,0]", None), ("446d1", "[1,-1,0,-4,4]", None),
          ("571b1", "[0,1,1,-4,2]", None), ("5077a1", "[0,0,1,-7,6]", None), ("655a1", "[0,0,1,-13,18]", None),
          ("1132a1", "[0,1,0,-5,4]", None), ("916c1", "[0,0,0,-4,1]", None), ("1088j1", "[0,1,0,-25,39]", None),
          ("433a1", "[1,0,0,0,1]", None)]
fails = []
# T1
H = pari('E=ellinit([0,1,1,-2,0]); mtmatrix(E, [[-1,1],[0,0]], 23, Mod(5,23))[1]')
t1 = str(H) == "[9, 7; 7, 9]"
print("T1 appendix 389a1 m=23:", H, "OK" if t1 else "FAIL"); fails += [] if t1 else ["T1"]
# T2, T3
n2 = n3 = 0
for lab, a, base in CURVES:
    pari(f"E=ellinit({a}); N=ellglobalred(E)[1]; ct=ellglobalred(E)[3]")
    pari(f"P={base}" if base else "P=ellsaturation(E, ellrank(E)[4], 1000)")
    pari("Pp=P[1]; for(k=2,#P,Pp=elladd(E,Pp,P[k]))")
    for m in [int(q) for q in pari("primes([7,100])")]:
        if int(pari(f"N % {m}")) == 0: continue
        pari(f"g=znprimroot({m}); X=mtmatrix(E,P,{m},g); M=X[2]; Href=pair_det(E,P,M,{m},g)[2]")
        ok2 = int(pari(f"my(n=({m}-1)/2,o=1);for(i=1,#P,for(j=1,#P,if((X[1][i,j]-Href[i,j])%n,o=0)));o"))
        n2 += 1; fails += [] if ok2 else [f"T2 {lab} m={m}"]
        ok3 = int(pari(f"my(R=vector(#P,i,ellmul(E,P[i],M)),o=1);for(i=1,#P,for(j=i,#P,my(u=mtpair(E,R[i],R[j],{m}),v=gPB(E,R[i],R[j],Pp,{m}));if(u!=v&&u!=-v,o=0)));o"))
        n3 += 1; fails += [] if ok3 else [f"T3 {lab} m={m}"]
    print(lab, "done", flush=True)
print(f"T2 mtmatrix == independent recipe: {n2 - sum(f.startswith('T2') for f in fails)}/{n2}")
print(f"T3 mtpair == Portillo-Bobadilla eq. (18) up to sign: {n3 - sum(f.startswith('T3') for f in fails)}/{n3}")
# T4
def raises(expr):
    try: pari(expr); return False
    except cypari2.PariError: return True
t4 = [raises('mtmatrix(ellinit([0,0,0,-4*16,16^3]), [[0,64]], 13)'),          # non-minimal model of 916c1-like curve
      raises('E=ellinit([0,1,1,-2,0]); mtpair(E, [-1,1], [0,0], 23)'),          # points not in formal group
      raises('E=ellinit([0,1,1,-2,0]); mtmatrix(E, [[-1,1],[0,0]], 389)')]     # m = conductor
print("T4 errors raised:", t4); fails += [] if all(t4) else ["T4"]
print("ALL PASS" if not fails else f"FAILURES: {fails[:10]}"); sys.exit(1 if fails else 0)
