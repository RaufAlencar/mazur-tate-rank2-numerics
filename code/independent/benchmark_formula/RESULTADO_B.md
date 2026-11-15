# Resultado B: forma fechada x fórmula com ponto auxiliar (15/11/2026)

Rodada única, script conferido por hash (HASHES.txt: PREREGISTRO_B.md 417293fb…, bench.py 4ce101f9…).
Saídas: resultado_B_m200.json (sha256 e04fae38…), saida_B.txt (4e861551…).

## Veredito
| Item | Resultado |
|---|---|
| C1: mesma matriz em F_m^×/±1, entrada a entrada | 172 de 172 níveis, 0 desacordos |
| CN: variante errada [A(R)] discorda na diagonal | 329/387 (≥ 80%: passa) |
| **B1: mediana de ρ = t_B / t_A ≥ 2** | **3.45: CONFIRMADO** (mín. 1.85, máx. 5.37) |
| B2 (descritivo): razão do tempo total, incluindo o cálculo dos R_i | 2.55 |

## Por curva e faixa de m (mediana de ρ)
| Curva | m | níveis | ρ |
|---|---|---|---|
| 389a1 | 7 a 50 | 12 | 4.64 |
| 389a1 | 51 a 100 | 10 | 4.39 |
| 389a1 | 101 a 200 | 21 | 4.51 |
| 433a1 | 7 a 50 | 12 | 3.83 |
| 433a1 | 51 a 100 | 10 | 3.57 |
| 433a1 | 101 a 200 | 21 | 3.50 |
| 709a1 | 7 a 50 | 12 | 3.57 |
| 709a1 | 51 a 100 | 10 | 3.44 |
| 709a1 | 101 a 200 | 21 | 3.30 |
| 5077a1 | 7 a 50 | 12 | 2.43 |
| 5077a1 | 51 a 100 | 10 | 2.50 |
| 5077a1 | 101 a 200 | 21 | 2.45 |

## Leitura
- A forma fechada calcula a matriz do emparelhamento cerca de 3.4 vezes mais rápido que a fórmula com ponto
  auxiliar, nos mesmos pontos, com resultado idêntico. A razão vem de evitar somas de pontos enormes
  (a diagonal não precisa de soma nenhuma: é o inverso do numerador de x).
- Contra o esperado no pré-registro, o cálculo dos R_i NÃO domina: o ganho no tempo total é de 2.6 vezes.
  As somas extras da fórmula B são feitas em pontos do tamanho final (até 13790 dígitos no denominador aqui),
  e custam tanto quanto multiplicar o ponto.
- A vantagem diminui no posto 3 (5077a1), onde a fração da matriz fora da diagonal é maior.
- Limites: 4 curvas, m ≤ 200, 1 núcleo, PARI 2.17; ponto auxiliar P′ = P₁ + P₂ fixo (outra escolha de P′
  muda o custo de B). É desempenho de implementação, não matemática nova.
