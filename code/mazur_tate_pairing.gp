\\ mazur_tate_pairing.gp: the Mazur-Tate pairing at a prime level m, by the closed formula of
\\ R. Alencar, "Numerical verification of the refined Mazur-Tate conjecture in ranks 0-3", v0.7, §2
\\ (doi:10.5281/zenodo.23964110):
\\
\\     <R,R>_m = [A(R)]^(-1),     <R,S>_m = [gamma_{R,S} / A(R+S)]     in F_m^x / {+-1},
\\
\\ where x(T) = A(T)/d(T)^2 in lowest terms and gamma_{R,S} is the constant term of the primitive integral
\\ line alpha*y + beta*x + gamma through R and S (the tangent line if R = S). It agrees with the
\\ auxiliary-point formula of Mazur-Tate (1987, p. 747) used by Portillo-Bobadilla (2019, eq. (18)),
\\ needs no auxiliary point, and in our measurements is 2.5 to 4.5 times faster (benchmark_formula).
\\ Status: derived from the biextension definition granting the Neron-model description of the
\\ canonical local splittings (note §2 and §7.1); not a theorem until that identification is checked.
\\ Orientation: this is the orientation of the note; Tan's (1995) is its inverse.
\\
\\ Hypotheses checked at run time: E is a global minimal model; R, S (and R+S) lie in E_S(Q), i.e. in
\\ the formal group at m (m divides d) and in E_0(Q_l) at every bad l (the caller guarantees the latter,
\\ e.g. by multiplying by c_Tam, as mtmatrix does).
\\
\\ PARI/GP 2.17. License MIT, as the rest of the repository.

mt_isminimal(E) = ellglobalred(E)[2] == [1, 0, 0, 0];

\\ A(T): numerator of x(T); d(T) = sqrt of the denominator of x(T).
mt_A(T) = numerator(T[1]);
mt_d(T) = my(q = denominator(T[1]), d = sqrtint(q)); if (d^2 != q, error("mt: denominator of x is not a square")); d;

\\ <R,S>_m as an element of F_m^x, well defined up to sign. R, S in E_S(Q), R, S, R+S != O.
mtpair(E, R, S, m) =
{
  my(lam, c, den, L, RS);
  if (R == [0] || S == [0], error("mt: point at infinity"));
  if (mt_d(R) % m || mt_d(S) % m, error("mt: point not in the formal group at m"));
  if (R == S, return(Mod(mt_A(R), m)^-1));
  lam = (S[2] - R[2]) / (S[1] - R[1]); c = R[2] - lam * R[1];
  den = lcm(denominator(lam), denominator(c));
  L = [den, -den * lam, -den * c]; L = L / content(L);    \\ den*y - den*lam*x - den*c = 0, primitive
  RS = elladd(E, R, S);
  if (RS == [0], error("mt: R + S = O"));
  Mod(L[3], m) / Mod(mt_A(RS), m);
}

\\ Discrete logarithm of a class of F_m^x/{+-1} in Z/n, n = (m-1)/2 (the sign is killed mod n).
mtlog(x, m, g) = znlog(x, g) % ((m - 1) / 2);

\\ For a basis P of E(Q)/tors and a prime m of good reduction: R_i = M P_i with M = c_Tam * #E(F_m),
\\ returns [H, M, n] with H_ij = log_g <R_i, R_j>_m in Z/n. Requires m odd, m not dividing N.
mtmatrix(E, P, m, g = 0) =
{
  my(r = #P, n = (m - 1) / 2, M, R, H = matrix(r, r));
  if (!mt_isminimal(E), error("mt: E must be a global minimal model (ellminimalmodel)"));
  if (ellglobalred(E)[1] % m == 0, error("mt: m divides the conductor"));
  if (!g, g = znprimroot(m));
  M = ellglobalred(E)[3] * ellcard(E, m);
  R = vector(r, i, ellmul(E, P[i], M));
  for (i = 1, r, H[i, i] = mtlog(mtpair(E, R[i], R[i], m), m, g));
  for (i = 1, r, for (j = i + 1, r, H[i, j] = mtlog(mtpair(E, R[i], R[j], m), m, g); H[j, i] = H[i, j]));
  [H, M, n];
}
