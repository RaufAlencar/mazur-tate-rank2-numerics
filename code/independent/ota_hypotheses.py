# Para os níveis m em que p | J_m (p = 5, 7), confere as hipóteses do Thm. 1.2 de Ota (2018) que dá para checar:
#  (a) E sem CM; (b) p não divide 6·N·#E(F_p)·prod c_l; (c) E(F_m)[p] cíclico;
#  (d) rho_p sobrejetiva: as curvas aqui têm condutor primo (semiestáveis), e para E semiestável rho_p é
#      sobrejetiva ou redutível (Serre 1972, Prop. 21); redutível <=> existe isogenia racional de grau p.
import cypari2, json, os
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
pari = cypari2.Pari(); pari.allocatemem(3 * 10**9)
pari('read("mt_indep.gp")')
pari('''subgrp_order(E, P, m) = {
  my(Em = ellinit(E, m), S = Map(), Q = List(), z = [0]);
  mapput(S, z, 1); listput(Q, z);
  my(G = vector(#P, i, [Mod(P[i][1], m), Mod(P[i][2], m)]), k = 1);
  while (k <= #Q, my(x = Q[k]); k++;
    for (i = 1, #G, my(y = elladd(Em, x, G[i])); if (!mapisdefined(S, y), mapput(S, y, 1); listput(Q, y))));
  #Q;
}''')
CURVES = [("389a1", "[0,1,1,-2,0]", "[[-1,1],[0,0]]", 10000), ("709a1", "[0,-1,1,-2,0]", None, 10000),
          ("433a1", "[1,0,0,0,1]", None, 3000), ("563a1", "[1,1,1,-15,16]", None, 3000),
          ("643a1", "[1,0,0,-4,3]", None, 3000), ("997b1", "[0,-1,1,-5,-3]", None, 3000),
          ("571b1", "[0,1,1,-4,2]", None, 3000), ("5077a1", "[0,0,1,-7,6]", None, 1500)]
rows = []
for lab, a, base, mmax in CURVES:
    pari(f'E = ellinit({a}); N = ellglobalred(E)[1]')
    pari(f'P = {base}' if base else 'P = ellsaturation(E, ellrank(E)[4], 1000)')
    N = int(pari('N'))
    CMJ = {0, 1728, -3375, 8000, -32768, 54000, 287496, -884736, -12288000, 16581375, -884736000, -147197952000, -262537412640768000}
    j = pari('E.j'); cm = 1 if (pari.denominator(j) == 1 and int(j) in CMJ) else 0
    prime_cond = bool(pari(f'isprime({N})'))
    for p in (5, 7):
        ap = int(pari(f'ellap(E, {p})')); Np = p + 1 - ap
        isog = int(pari(f'#ellisomat(E, {p}, 1)[1]')) > 1   # existe isogenia racional de grau p?
        for m in [int(q) for q in pari(f'primes([7, {mmax}])')]:
            n = (m - 1) // 2
            if N % m == 0 or n % p: continue
            Nm = int(pari(f'ellcard(E, {m})'))
            if Nm % p: continue
            J = Nm // int(pari(f'subgrp_order(E, P, {m})'))
            if J % p: continue
            cyc = pari(f'ellgroup(ellinit(E, {m}))')
            p_rank = sum(1 for d in cyc if int(d) % p == 0)
            rows.append(dict(curva=lab, p=p, m=m, sem_CM=(cm == 0), p_nao_divide_6N=(6 * N) % p != 0,
                             p_nao_divide_EFp=Np % p != 0, condutor_primo=prime_cond, isogenia_grau_p=isog,
                             EFm_p_ciclico=(p_rank <= 1)))
for r in rows: print(r, flush=True)
ok = all(r["sem_CM"] and r["p_nao_divide_6N"] and r["p_nao_divide_EFp"] and r["condutor_primo"]
         and not r["isogenia_grau_p"] and r["EFm_p_ciclico"] for r in rows)
print("níveis:", len(rows), "| todas as hipóteses conferidas valem:", ok)
json.dump(rows, open("ota_hypotheses_results.json", "w"), indent=1)
