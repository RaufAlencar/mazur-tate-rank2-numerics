# Registro: reproduz a comparação Portillo-Bobadilla (2019) eq. (18) x nossa fórmula, 389a1, m = 23..41, três P'.
import cypari2, pathlib
pari = cypari2.Pari(); pari.allocatemem(10**9)
src = (pathlib.Path(__file__).parent / "portillo.gp")
pari('read("%s")' % src.as_posix())
pari('E=ellinit([0,1,1,-2,0]); P=[[-1,1],[0,0]]')
ok = tot = 0
for Pp in ("[1,0]", "[3,5]", "[-2,0]"):
    if not int(pari(f'ellisoncurve(E,{Pp})')): print("P' fora da curva", Pp); continue
    for m in (23, 29, 31, 37, 41):
        r = pari(f'my(N=ellcard(E,{m}),R=ellmul(E,P[1],N),S=ellmul(E,P[2],N)); [gPB(E,R,R,{Pp},{m})/ours(E,R,R,{m}), gPB(E,S,S,{Pp},{m})/ours(E,S,S,{m}), gPB(E,R,S,{Pp},{m})/ours(E,R,S,{m})]')
        vals = [int(pari.lift(x)) for x in r]
        good = [v in (1, m - 1) for v in vals]
        ok += sum(good); tot += 3
        print(f"P'={Pp} m={m} razao PB/nossa = {vals} (+-1? {good})")
print(f"razao = +-1 em {ok}/{tot}")
