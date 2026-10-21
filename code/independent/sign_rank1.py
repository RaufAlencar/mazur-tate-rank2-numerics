# Decide o sinal de u_E em posto 1. Testa u = -2c e u = +2c e conta acertos de cada um (nada é ajustado).
import cypari2, json, os
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
pari('read("mt_indep.gp")')
C = {"37a1":"[0,0,1,-1,0]","43a1":"[0,1,1,0,0]","53a1":"[1,-1,1,0,0]","57a1":"[0,-1,1,-2,2]",
     "58a1":"[1,-1,0,-1,1]","61a1":"[1,0,0,-2,1]","77a1":"[0,0,1,2,0]","79a1":"[1,1,1,-2,0]"}
out = {}
for lab, a in C.items():
    pari(f'E = ellinit({a}); MS = msinit_curve(E); N = ellglobalred(E)[1]; P = ellsaturation(E, ellrank(E)[4], 1000)')
    assert int(pari('#P')) == 1 and int(pari('elltors(E)[1]')) == 1
    cinf = 2 if pari('E.disc') > 0 else 1; ctam = int(pari('ellglobalred(E)[3]'))
    hits = {"-2c": 0, "+2c": 0}; nontriv = 0; notin = 0
    for m in [int(q) for q in pari('primes([7, 400])')]:
        if int(pari(f'N % {m}')) == 0: continue
        res = pari(f'mt_check(E, MS, P, {ctam}, {m})')
        if 'NOT' in str(res[2]): notin += 1; continue
        n, lhs, dH = int(res[1]), int(res[4]), int(res[3])
        if lhs == 0 and dH == 0: continue
        nontriv += 1
        if lhs == (-2 * cinf * dH) % n: hits["-2c"] += 1
        if lhs == (2 * cinf * dH) % n: hits["+2c"] += 1
    out[lab] = dict(c_inf=cinf, c_Tam=ctam, nao_triviais=nontriv, acertos=hits, theta_fora_de_I=notin)
    print(lab, out[lab], flush=True)
json.dump(out, open("sign_results_rank1.json", "w"), indent=1)
