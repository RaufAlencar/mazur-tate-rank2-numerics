# Testa a implicação do teorema de Ota (provada para p >= 5, sob hipóteses):
#   termo líder c_r(m) != 0 mod p  ==>  p não divide J_m,
# ou, na forma contrapositiva testável:  p | J_m  ==>  c_r(m) == 0 mod p.
# Nível primo m com p | n = (m-1)/2. Curvas com c_Tam = 1 e torção trivial, para que J_m = #E(F_m) / #<P_i mod m>.
import cypari2, json, os, sys
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
CURVES = [("389a1", "[0,1,1,-2,0]", "[[-1,1],[0,0]]", 3000), ("709a1", "[0,-1,1,-2,0]", None, 3000),
          ("433a1", "[1,0,0,0,1]", None, 3000), ("563a1", "[1,1,1,-15,16]", None, 3000),
          ("643a1", "[1,0,0,-4,3]", None, 3000), ("997b1", "[0,-1,1,-5,-3]", None, 3000),
          ("571b1", "[0,1,1,-4,2]", None, 3000), ("5077a1", "[0,0,1,-7,6]", None, 1500)]
PS = [3, 5, 7, 11, 13]
out = {}
for lab, a, base, mmax in CURVES:
    pari(f'E = ellinit({a}); MS = msinit_curve(E); N = ellglobalred(E)[1]')
    pari(f'P = {base}' if base else 'P = ellsaturation(E, ellrank(E)[4], 1000)')
    assert int(pari('ellglobalred(E)[3]')) == 1 and int(pari('elltors(E)[1]')) == 1
    r = int(pari('#P')); stats = {p: dict(p_divide_J=0, c_zero_quando_p_divide_J=0, violacoes=[]) for p in PS}
    for m in [int(q) for q in pari(f'primes([7, {mmax}])')]:
        n = (m - 1) // 2
        if int(pari(f'N % {m}')) == 0 or not any(n % p == 0 for p in PS): continue
        Np = int(pari(f'ellcard(E, {m})'))
        if not any(n % p == 0 and Np % p == 0 for p in PS): continue   # p | J_m exige p | #E(F_m)
        J = Np // int(pari(f'subgrp_order(E, P, {m})'))
        if not any(n % p == 0 and J % p == 0 for p in PS): continue
        cl = pari(f'theta_class(theta_coefs(MS, {m}, znprimroot({m})), {n}, {r})')
        if not int(cl[0]): continue
        c = int(cl[1])
        for p in PS:
            if n % p: continue
            s = stats[p]
            if J % p == 0:
                s["p_divide_J"] += 1
                if c % p == 0: s["c_zero_quando_p_divide_J"] += 1
                else: s["violacoes"].append(m)
    out[lab] = {"r": r, **{str(p): v for p, v in stats.items()}}
    print(lab, out[lab], flush=True)
json.dump(out, open("ota_implication_results.json", "w"), indent=1)
