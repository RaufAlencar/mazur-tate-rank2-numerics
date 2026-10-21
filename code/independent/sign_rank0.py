# Posto 0: testa M^(-1) * a_0 == u (mod n) para as duas hipóteses registradas. a_0 = soma dos símbolos [a/m]+.
import cypari2, json, os
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
C = {"11a1":"[0,-1,1,-10,-20]","14a1":"[1,0,1,4,-6]","15a1":"[1,1,1,-10,-10]","32a1":"[0,0,0,4,0]",
     "67a1":"[0,1,1,-12,-21]","75a1":"[0,-1,1,-8,-7]","571a1":"[0,-1,1,-929,-10595]","681b1":"[1,1,0,-1154,-15345]"}
SHA = {"571a1": 4, "681b1": 9}
out = {}
for lab, a in C.items():
    pari(f'E = ellinit({a}); MS = msfromell(E, 1); N = ellglobalred(E)[1]')
    cinf = 2 if pari('E.disc') > 0 else 1; ctam = int(pari('ellglobalred(E)[3]')); tors = int(pari('elltors(E)[1]'))
    sha = SHA.get(lab, 1)
    pA, pB = cinf * sha, -cinf * sha          # previsões para u * tors^2
    hit = {"H_A(+)": 0, "H_B(-)": 0}; tested = 0
    for m in [int(q) for q in pari('primes([7, 400])')]:
        if int(pari(f'N % {m}')) == 0: continue
        n = (m - 1) // 2
        a0 = pari(f'sum(a = 1, {m}-1, mseval(MS[1], MS[2], [a/{m}, oo]))')
        num, den = int(pari.numerator(a0)), int(pari.denominator(a0))  # com torção, a_0 é racional
        M = ctam * (m + 1 - int(pari(f'ellap(E, {m})')))
        if gcd(M, n) != 1 or gcd(den, n) != 1: continue
        lhs = (pow(M, -1, n) * num * pow(den, -1, n) * tors * tors) % n
        if lhs == 0 and pA % n == 0: continue
        tested += 1
        if lhs == pA % n: hit["H_A(+)"] += 1
        if lhs == pB % n: hit["H_B(-)"] += 1
    out[lab] = dict(c_inf=cinf, c_Tam=ctam, tors=tors, Sha=sha, niveis=tested, acertos=hit)
    print(lab, out[lab], flush=True)
json.dump(out, open("sign_results_rank0.json", "w"), indent=1)
