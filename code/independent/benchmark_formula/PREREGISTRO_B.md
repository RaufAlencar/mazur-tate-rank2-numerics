# Pré-registro B: velocidade da forma fechada contra a fórmula com ponto auxiliar

Registrado em 15/11/2026 pela sessão de registro do BSD, antes de qualquer medição de tempo. Nenhum
tempo das duas fórmulas foi medido antes deste arquivo (a conferência de valores 45/45 de 04/11 não
mediu tempo).

## Objeto
Calcular a matriz do emparelhamento de Mazur–Tate H = (⟨R_i, R_j⟩) em F_m^×/±1, com R_i = M·P_i,
M = c_Tam·#E(F_m), por dois métodos, a partir dos MESMOS pontos R_i já calculados:
- **A (forma fechada, nota v0.7 §2):** diagonal ⟨R,R⟩ = [A(R)]⁻¹; fora da diagonal ⟨R,S⟩ = [γ/A(R+S)],
  γ = termo constante da reta primitiva por R e S. Custo: 1 soma de pontos por par fora da diagonal.
- **B (Mazur–Tate 1987, p. 747; Portillo-Bobadilla 2019, eq. (18)):** g(R,S,P′) = d(P′+R)·d(P′+S) /
  (d(P′)·d(P′+R+S)), com P′ = P₁ + P₂ (ponto auxiliar fixo, fora do grupo formal). Custo: somas
  P′+R_i, R_i+R_j (inclui duplicação na diagonal) e P′+(R_i+R_j).

## Dados
Curvas 389a1, 433a1, 709a1 (posto 2) e 5077a1 (posto 3), base saturada do ellrank (389a1: [[-1,1],[0,0]]).
Níveis: primos 7 ≤ m ≤ 200 com m ∤ N. Implementação: PARI/GP 2.17 via cypari2, WSL, 1 núcleo.

## Medidas
- t_pts: tempo de calcular os R_i (comum aos dois métodos; descritivo).
- t_A, t_B: tempo da matriz inteira por cada método, mínimo de 3 repetições, mesmo processo.
- razão ρ = t_B / t_A por nível.

## Hipóteses e critério
- **C1 (consistência, deve valer sempre):** os dois métodos dão a mesma matriz em F_m^×/±1, entrada a
  entrada, em todos os níveis. Qualquer desacordo é investigado antes de falar de tempo.
- **CN (controle negativo):** a variante errada da forma fechada com [A(R)] no lugar de [A(R)]⁻¹ deve
  discordar de B em pelo menos 80% das entradas diagonais (se não, o teste de C1 não discrimina).
- **B1 (velocidade):** a mediana de ρ sobre todos os níveis é ≥ 2. CONFIRMADO se ≥ 2; NÃO CONFIRMADO
  se < 2. O valor exato é relatado de qualquer jeito, por curva e por faixa de m.
- **B2 (descritivo, sem veredito):** peso da economia no cálculo total, (t_pts + t_B)/(t_pts + t_A).
  Esperamos que t_pts domine e que a vantagem total seja pequena; isso será dito como estiver.

## Execução
Rodada única com o script congelado (sha256 em HASHES.txt). Orçamento: 30 min de parede, timeout
externo 45 min. Saída: resultado_B.json + resultado_B.md.

## Honestidade
É medida de desempenho de duas fórmulas para o MESMO objeto, não teste de conjectura. A vantagem, se
houver, é de implementação; não muda nada da matemática.
