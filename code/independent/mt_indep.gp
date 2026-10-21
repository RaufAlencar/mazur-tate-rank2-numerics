\\ Implementação independente da identidade de Mazur–Tate da nota (rev 80),
\\ escrita só a partir do texto da nota (§1–§3), sem usar os scripts da sessão de pesquisa.
\\ Teste: 2^r * M^(2r-1) * c_r(m) == u_E * det(H)  (mod n),  n = (m-1)/2.

\\ Símbolo modular plus [a/m]^+ normalizado como msfromell (sinal 1).
msinit_curve(E) = msfromell(E, 1);

\\ Coeficientes c_t de theta_m = sum_a [a/m]^+ sigma_a em Z[G_m^+], t = log_g(a) mod n.
theta_coefs(MS, m, g) =
{
  my(M = MS[1], xp = MS[2], n = (m-1)/2, c = vector(n));
  for (a = 1, m-1,
    my(t = znlog(Mod(a, m), g) % n, v = mseval(M, xp, [a/m, oo]));
    c[t+1] += v);
  c;
}

\\ Classe c_r(m) em I^r/I^(r+1) ≅ Z/n pelas reduções sucessivas de §1.
\\ Devolve [ok, c_r]: ok = 0 se theta_m não está em I^r.
theta_class(c, n, r) =
{
  my(K = r + n + 2, a = vector(K + 1));
  for (k = 0, K, a[k+1] = sum(t = 0, n-1, c[t+1] * binomial(t, k)));
  if (a[1] != 0, return([0, 0]));
  for (j = 1, r - 1,
    if (denominator(a[j+1]) != 1 || a[j+1] % n != 0, return([0, j]));
    my(kj = a[j+1] / n);
    for (i = 2, K - j + 1, a[j+i] -= kj * binomial(n, i)));
  if (denominator(a[r+1]) != 1, return([0, -1]));
  [1, a[r+1] % n];
}

\\ Meia-pareação h(R) = log_g[t/d] = log_g[-A/B] mod n, com x = A/d^2, y = B/d^3.
halfpair(E, R, m, g) =
{
  my(x = R[1], y = R[2], d = sqrtint(denominator(x)), A, B, n = (m-1)/2);
  if (d^2 != denominator(x), error("denominador de x não é quadrado"));
  A = numerator(x); B = y * d^3;
  if (denominator(B) != 1, error("B não inteiro"));
  if (A % m == 0 || B % m == 0, error("A ou B divisível por m"));
  znlog(Mod(-A, m) / Mod(B, m), g) % n;
}

\\ Matriz H da polarização para R_i = M P_i, e o det mod n.
pair_det(E, P, M, m, g) =
{
  my(r = #P, n = (m-1)/2, R = vector(r, i, ellmul(E, P[i], M)), h = vector(r, i, halfpair(E, R[i], m, g)), H = matrix(r, r));
  for (i = 1, r, H[i,i] = 2*h[i]);
  for (i = 1, r, for (j = i+1, r,
    H[i,j] = halfpair(E, elladd(E, R[i], R[j]), m, g) - h[i] - h[j]; H[j,i] = H[i,j]));
  [matdet(H) % n, H];
}

\\ Teste num nível primo m. Devolve [m, n, c_r, detH, lhs, rhs_sem_u, M].
mt_check(E, MS, P, ctam, m) =
{
  my(r = #P, n = (m-1)/2, g = znprimroot(m), Np = m + 1 - ellap(E, m), M = ctam * Np, cl, dH, lhs);
  cl = theta_class(theta_coefs(MS, m, g), n, r);
  if (!cl[1], return([m, n, "NOT_IN_I^r", cl[2]]));
  dH = pair_det(E, P, M, m, g)[1];
  lhs = (2^r * M^(2*r-1) * cl[2]) % n;
  [m, n, cl[2], dH, lhs, M];
}
