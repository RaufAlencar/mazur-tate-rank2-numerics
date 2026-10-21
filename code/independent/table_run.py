# Roda a implementação independente (mt_indep.gp) nas 10 curvas da tabela da nota (§4),
# com u_E = (-2)^r * c_inf / #tors^2 FIXO pela fórmula (sem ajustar constante aos dados).
import cypari2, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
pari = cypari2.Pari(); pari.allocatemem(4 * 10**9)
pari('read("mt_indep.gp")')
TABLE = [  # (rótulo, a-invariantes, base usada na nota ou None = saturada do ellrank, m máximo)
    ("389a1",  "[0,1,1,-2,0]",    "[[-1,1],[0,0]]",            600),
    ("709a1",  "[0,-1,1,-2,0]",   None,                        400),
    ("446d1",  "[1,-1,0,-4,4]",   None,                        400),
    ("571b1",  "[0,1,1,-4,2]",    None,                        330),
    ("5077a1", "[0,0,1,-7,6]",    None,                        330),
    ("655a1",  "[0,0,1,-13,18]",  None,                        330),
    ("1132a1", "[0,1,0,-5,4]",    None,                        330),
    ("916c1",  "[0,0,0,-4,1]",    None,                        330),
    ("1088j1", "[0,1,0,-25,39]",  None,                        300),
    ("433a1",  "[1,0,0,0,1]",     None,                        200),
]
only = sys.argv[1:]  # opcional: rótulos
out = {}
for label, a, base, mmax in TABLE:
    if only and label not in only: continue
    t0 = time.time()
    pari(f'E = ellinit({a}); MS = msinit_curve(E); N = ellglobalred(E)[1]')
    r = int(pari.ellanalyticrank(pari('E'))[0])
    cinf = 2 if pari('E.disc') > 0 else 1
    tors = int(pari('elltors(E)[1]'))
    ctam = int(pari('ellglobalred(E)[3]'))
    pari(f'P = {base}' if base else 'P = ellsaturation(E, ellrank(E)[4], 1000)')
    assert int(pari('#P')) == r
    num = (-2) ** r * cinf
    ok = bad = triv = notin = weak = 0; fails = []
    for m in [int(q) for q in pari(f'primes([7, {mmax}])')]:
        if int(pari(f'N % {m}')) == 0: continue
        res = pari(f'mt_check(E, MS, P, {ctam}, {m})')
        if 'NOT' in str(res[2]): notin += 1; fails.append((m, 'theta fora de I^r')); continue
        n, lhs, dH = int(res[1]), int(res[4]), int(res[3])
        # u_E = num / tors^2: compara tors^2 * lhs com num * det (evita dividir mod n)
        L, R = (tors * tors * lhs) % n, (num * dH) % n
        if L == 0 and R == 0: triv += 1
        elif L == R: ok += 1
        else: bad += 1; fails.append((m, lhs, dH, n))
    out[label] = dict(r=r, c_inf=cinf, tors=tors, c_Tam=ctam, u_E=f"{num}/{tors*tors}", m_max=mmax,
                      ok=ok, falhas=bad, triviais=triv, fora_de_Ir=notin, detalhes=fails[:10], seg=round(time.time()-t0))
    print(label, out[label], flush=True)
json.dump(out, open("table_results.json", "w"), indent=1)
