# Resultados numéricos desta entrega (reproduzíveis com os scripts ao lado)

## Constante u_E (mt_rank_run.py; u em 2^r M^{2r-1} c_r ≡ u·det H mod n)
| curva | ainvs | r | c_Tam | c_∞ | tors | u | níveis | não triviais | falhas |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 655a1 | [0,0,1,-13,18] | 2 | 2 | 1 | 1 | 4 | 23–330 | 51 | 0 |
| 1132a1 | [0,1,0,-5,4] | 2 | 3 | 1 | 1 | 4 | 23–330 | 45 | 0 |
| 916c1 | [0,0,0,-4,1] | 2 | 3 | 2 | 1 | 8 | 23–330 | 46 | 0 |
| 571b1 | [0,1,1,-4,2] | 2 | 1 | 1 | 1 | 4 | 23–330 | 54 | 0 |
| 1088j1 | [0,1,0,-25,39] | 2 | 4 | 2 | 2 | 2 | 23–300 | 44 | 0 |
| 433a1 | [1,0,0,0,1] | 2 | 1 | 1 | 1 | 4 | 7–113 | 21 | 0 |
Lei: u_E = (−2)^r · c_∞ / #tors². (389a1, 709a1, 446d1: 8; 5077a1 posto 3: −16, logs anteriores.)

## Níveis de Bley–Macías Castillo (bmc_levels.py)
389a1: m=7 detH≡2 (mod 3), m=13 ≡1, m=43 ≡2 — identidade mod n OK nos três (4M³c₂ = 8detH: 1=1, 2=2, 19=19).
446d1: m=19 detH≡1 (mod 3), identidade OK (2=2).
433a1: m=7 ≡2, m=13 ≡2, m=31 ≡1, m=37 ≡1 — identidade OK com u=4 (2=2, 2=2, 10=10, 10=10).
Todos não nulos mod 3, como o código deles exige (det Ψ ≠ 0 mod 3).

## Nível composto (mt_ota_run.py)
389a1, 7 ≤ m₁ < m₂ ≤ 61, gcd(n₁,n₂) ≥ 3: 28 pares; θ_S ∈ I²; checks de norma x₁², x₂² OK em todos; coeficiente misto:
4M³a₁₁ ≡ −8·D₁₂ (mod g) em todos, não trivial em 10 pares (g ∈ {3,5,6,7}); (11,31) mod 5 e (29,43) mod 7 excluem u = −8 e ±4, ±16.
433a1 (Den=2): 26/28 pares no momento da entrega, 10 não triviais, todos consistentes com 4M³a₁₁ ≡ −8·D₁₂ (= −Den·u_E), 0 falhas.

## Biextensão (§2 da nota)
- biext.gp: E·c_R/A(2R) ≡ ±(t/d)² (mod m) em 389a1, 571b1, 446d1, m ∈ {23,29,31,37}, ambos geradores (24/24).
- chkg2.gp: −A³/(A(2R)·c_L) ≡ ±(t/d)² e δ = c_L = 1 em 349 pontos R_i = M·P_i, 7 ≤ m ≤ 113, curvas 389a1, 446d1, 1088j1, 916c1, 655a1, 5077a1.
- chord.gp: [γ_{R,S}/A(R+S)] = polarização (t/d)(R+S)/((t/d)(R)(t/d)(S)) em 129 pares, 5 curvas, 11 ≤ m ≤ 113.
- chkg.gp (fora de E₀): 446d1 ímpar·P₁ dá δ = c_L = 2; 1088j1 ímpar·P_i dá δ = c_L = 8 — a hipótese R ∈ E₀ é necessária.
- gsweep.py 400: 453 curvas de posto ≥ 1, 3631 pontos c_Tam·k·P (k ≤ 8), δ = c_L = 1 em todos.

## Posto 4
234446a1: msfromell estourou 7 GB (OOM) antes do primeiro nível (rank4.log).
