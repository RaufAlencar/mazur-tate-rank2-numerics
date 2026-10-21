# Compara o símbolo [0]^+ = mseval(..., [0, oo]) do PARI com L(E,1)/omega_1 (sinal da convenção).
import cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9); pari.default("realprecision", 30)
for lab, a in [("67a1","[0,1,1,-12,-21]"),("75a1","[0,-1,1,-8,-7]"),("571a1","[0,-1,1,-929,-10595]"),("11a1","[0,-1,1,-10,-20]")]:
    pari(f"E = ellinit({a}); MS = msfromell(E, 1)")
    s0 = pari("mseval(MS[1], MS[2], [0, oo])"); s0r = pari("mseval(MS[1], MS[2], [oo, 0])")
    L = pari("lfun(E, 1)"); w1 = pari("E.omega[1]")
    print(lab, "[0,oo] ->", s0, " [oo,0] ->", s0r, " L(E,1)/omega1 =", pari(f"{L}/{w1}"), " MS[3] (normalizacao) =", pari("MS[3]") if int(pari("#MS")) >= 3 else "-")
