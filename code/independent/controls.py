# Controles negativos: o teste precisa FALHAR com constante errada, sem c_Tam, ou com theta de outra curva.
import cypari2, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "..", "results", "independent")
os.chdir(HERE)
pari = cypari2.Pari(); pari.allocatemem(3 * 10**9)
pari('read("mt_indep.gp")')
pred = json.load(open(os.path.join(RES, "blind_predictions.json")))["curvas"]

def run(label, u, ctam, theta_from=None, mmax=200):
    c = pred[label]
    pari(f'E = ellinit({c["ainvs"]}); P = {c["base"]}; N = ellglobalred(E)[1]')
    src = pred[theta_from]["ainvs"] if theta_from else c["ainvs"]
    pari(f'E2 = ellinit({src}); MS = msinit_curve(E2); N2 = ellglobalred(E2)[1]')
    ok = bad = 0
    for m in pari(f'primes([7, {mmax}])'):
        m = int(m)
        if int(pari(f'(N*N2) % {m}')) == 0: continue
        res = pari(f'mt_check(E, MS, P, {ctam}, {m})')
        if 'NOT' in str(res[2]): bad += 1; continue
        n, lhs, dH = int(res[1]), int(res[4]), int(res[3])
        rhs = (u * dH) % n
        if lhs == 0 and rhs == 0: continue
        if lhs == rhs: ok += 1
        else: bad += 1
    return ok, bad

tests = [
  ("563a1 com u certo = 4", "563a1", 4, 1, None),
  ("563a1 com u = -4", "563a1", -4, 1, None),
  ("563a1 com u = 8", "563a1", 8, 1, None),
  ("563a1 com u = 2", "563a1", 2, 1, None),
  ("664a1 certo (c_Tam = 4)", "664a1", 4, 4, None),
  ("664a1 SEM c_Tam (=1)", "664a1", 4, 1, None),
  ("664a1 com c_Tam = 2", "664a1", 4, 2, None),
  ("563a1 pontos + theta de 643a1", "563a1", 4, 1, "643a1"),
  ("11197a1 com u = +16", "11197a1", 16, 1, None),
]
for name, lab, u, ct, tf in tests:
    ok, bad = run(lab, u, ct, tf, 150 if lab == "11197a1" else 200)
    print(f"{name:34s} acertos {ok:3d}  falhas {bad:3d}", flush=True)
