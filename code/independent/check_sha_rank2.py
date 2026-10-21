# Conferência independente das curvas de posto 2 com #Ш > 1 da sessão sha-rank2.
import cypari2
pari = cypari2.Pari(); pari.allocatemem(3 * 10**9); pari.default("realprecision", 30)
C = {"53295337.a1": "[1,-1,0,94,9]", "136116413.a1": "[0,0,1,-314,-2067]",
     "499290.e2": "[1,1,0,-11458,-452978]", "226005.b5": "[1,1,1,-268695,67104702]"}
for lab, a in C.items():
    pari(f"E = ellinit({a})")
    N = pari("ellglobalred(E)[1]"); rk = pari("ellrank(E)")
    G = pari("ellsaturation(E, ellrank(E)[4], 1000)")
    an = pari("ellanalyticrank(E)")
    r = int(an[0]); Lr = an[1] / pari(f"factorial({r})")
    reg = pari.matdet(pari.ellheightmatrix(pari("E"), G)) if len(G) == 2 else None
    sha = Lr / (pari("ellbsd(E)") * reg) if reg else None
    print(lab, "N =", N, "| ellrank [r,R,s,...] =", [rk[0], rk[1], rk[2]], "| #gens =", len(G),
          "| r_an =", r, "| Sha_an =", pari.round(sha) if sha is not None else None, "(", sha, ")",
          "| tors =", pari("elltors(E)[1]"), "| c_Tam =", pari("ellglobalred(E)[3]"), flush=True)
