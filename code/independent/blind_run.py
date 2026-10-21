# PASSO 2 do teste cego: roda a identidade com as previsões já gravadas em blind_predictions.json.
import cypari2, json, sys, time, os
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "..", "results", "independent")
os.chdir(HERE)
pari = cypari2.Pari()
pari.allocatemem(3 * 10**9)
pari('read("mt_indep.gp")')
pred = json.load(open(os.path.join(RES, "blind_predictions.json")))["curvas"]
MMAX = {2: int(sys.argv[1]) if len(sys.argv) > 1 else 200, 3: int(sys.argv[2]) if len(sys.argv) > 2 else 150}
summary = {}
for label, c in pred.items():
    t0 = time.time()
    pari(f'E = ellinit({c["ainvs"]}); MS = msinit_curve(E); P = {c["base"]}; N = ellglobalred(E)[1]')
    u, r, ctam = c["u_E_previsto"], c["r"], c["c_Tam"]
    ok = bad = triv = notin = 0; fails = []
    for m in pari(f'primes([7, {MMAX[r]}])'):
        m = int(m)
        if int(pari(f'N % {m}')) == 0: continue
        res = pari(f'mt_check(E, MS, P, {ctam}, {m})')
        if str(res[2]) == '"NOT_IN_I^r"' or str(res[2]) == 'NOT_IN_I^r':
            notin += 1; fails.append((m, "theta fora de I^r")); continue
        n, lhs, dH = int(res[1]), int(res[4]), int(res[3])
        rhs = (u * dH) % n
        if lhs == 0 and rhs == 0: triv += 1
        elif lhs == rhs: ok += 1
        else: bad += 1; fails.append((m, lhs, rhs, n))
    summary[label] = dict(r=r, u_previsto=u, nao_triviais_ok=ok, falhas=bad, triviais=triv, fora_de_Ir=notin,
                          detalhes_falhas=fails[:10], segundos=round(time.time() - t0))
    print(label, summary[label], flush=True)
json.dump(summary, open(os.path.join(RES, "blind_results_rerun.json"), "w"), indent=1)
