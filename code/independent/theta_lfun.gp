\\ θ_m pelas L-funções torcidas do PARI (lfuntwist), independente da série escrita pela sessão sha-rank2.
\\ S_k := Σ_a [a/m]^+ χ_k(a) = κ · τ(χ̄_k) · L(E,χ_k,1) / ω1, com κ a calibrar contra msfromell.
\\ Inversão: [a/m]^+ = (1/(m-1)) Σ_k S_k χ̄_k(a)  (k ímpar dá 0; k = 0 dá Σ[a/m] = (a_m-2)[0] = 0 em posto ≥ 1).

\\ Caractere χ_k relativo ao gerador g de (Z/m)^×: χ_k(g^j) = e(kj/(m-1)).
chik(k, a, m, g) = exp(2*Pi*I*k*znlog(Mod(a,m), g)/(m-1));
gauss_bar(k, m, g) = sum(a = 1, m-1, conj(chik(k, a, m, g)) * exp(2*Pi*I*a/m));

\\ Vetor de S_k (k par, 2 ≤ k ≤ m-3) a partir das L-funções torcidas, sem o fator κ.
Svec_L(E, m) =
{
  my(G = znstar(m, 1), g = Mod(G.gen[1], m), w1 = real(E.omega[1]), S = vector(m-1));
  for (k = 1, m-2, if (k % 2, next);
    my(L = lfun(lfuntwist(E, [G, [k]]), 1));
    S[k+1] = gauss_bar(k, m, g) * L / w1);
  [S, g];
}

\\ Vetor de S_k a partir dos símbolos de msfromell (para calibrar κ).
Svec_ms(E, MS, m, g) = vector(m-1, kk, my(k = kk-1); sum(a = 1, m-1, mseval(MS[1], MS[2], [a/m, oo]) * chik(k, a, m, g)));

\\ Símbolos [a/m]^+ (a = 1..m-1) a partir de S e κ, reconhecidos como racionais de denominador ≤ 1000.
symbols_from_S(S, kappa, m, g) =
{
  my(v = vector(m-1, a, real(sum(kk = 1, m-1, kappa*S[kk]*conj(chik(kk-1, a, m, g)))/(m-1))));
  [vector(m-1, a, bestappr(v[a], 1000)), vecmax(vector(m-1, a, abs(v[a] - bestappr(v[a], 1000))))];
}

\\ Coeficientes c_t de θ_m em Z[G_m^+] com t = log_g(a) mod n, dados os símbolos.
theta_coefs_from_symbols(sym, m, g) =
{
  my(n = (m-1)/2, c = vector(n));
  for (a = 1, m-1, my(t = znlog(Mod(a, m), g) % n); c[t+1] += sym[a]);
  c;
}
