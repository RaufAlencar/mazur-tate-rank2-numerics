# Pré-registro: classes de Kolyvagin de duas pernas em posto 3 (p = 3)

Data/hora do registro: 2026-11-03T12:31:23-03:00 (antes de qualquer cálculo de κ_n ou κ_ℓ em posto 3).

O que já tinha rodado antes deste registro:
- o controle de posto 2 (389a1, D = −7, ℓ = 5, kgen3 inalterado; resultado: reta prevista = calculada);
- um teste de depuração do código novo em 11642a1, em precisão 120, que não chegou ao nível n (falhou antes, nas pernas).

Nenhum valor de κ foi visto.

## Arquivos congelados

| arquivo | sha256 |
|---|---|
| prereg/previsoes.json | 1f0920daca7e4f89530732da1f84503c2ec882f2518d846f96875e1b5e496b34 |
| prereg/previsao.py | 984fcbdc91a7863e9ea3dbfd66ccfb5372d40fba87ed4b8a9daf4d84eae98285 |

Código de cálculo no momento do registro (ainda em depuração; toda mudança posterior fica no log da execução e na seção "Adendos" do relatório):

| arquivo | sha256 |
|---|---|
| kgen2p.py | 012c627b812ac862116f4653591e92925ba286c15d5ef47adde36e418a7d79b5 |
| two.gp | 4c02c61fe3d66651d97e0ae832905a81fa2af4d5d8e92b3091f11d693d24be5d |
| hreps.py | c3a544a35aaeea432eec6dd408ec30165aeb9fdb825bee10f9ce8cd5b9642aab |
| divp.gp | 36b62d043fffd7a98b67c5a518a0f3e52d23caf24694c72a662d3f6c60b3aaaf |

A previsão não depende do código de cálculo: sai só da redução dos geradores saturados de `ellrank` em E(F_ℓ).

## Hipóteses verificadas (work/hyp.py)

Nos três casos:
- E[3] é irredutível (ψ₃ sem raiz em K);
- 3 ∤ c_Tam(E);
- Ш(E/Q) = 1 (tabela de Cremona);
- E^D tem posto analítico 0, torção trivial e Ш_an(E^D) ∈ {1, 121, 100}, nenhum divisível por 3;
- K = Q(√D) tem h = 1 e D < −4.

Logo Sel₃(E/K)⁺ = E(Q)/3 ≅ F₃³ e Sel₃(E/K)⁻ = 0.

## Previsões (base = geradores saturados de ellrank, coordenadas (a,b,c) ↔ aP₁+bP₂+cP₃)

| # | curva, K, n | previsão para κ_n | por quê |
|---|---|---|---|
| P1 | 11642a1, D=−7, 5·17 | κ_n ∈ ⟨(1,2,2)⟩: ou é 0, ou gera essa reta | ker red₅ ∩ ker red₁₇ tem dimensão 1 |
| P2 | 5077a1, D=−163, 11·23 | κ_n ∈ ⟨(1,2,1)⟩: ou é 0, ou gera essa reta | idem |
| P3 | 18097a1, D=−163, 11·29 | κ_n = 0 (forçado) | a interseção tem dimensão 2, então λ(n) ≥ 2 |
| C1 | os três casos, cada perna ℓ sozinha | κ_ℓ = 0 (forçado) | κ_ℓ ∈ Sel₃(E^D/Q) = 0, porque y_K é de torção e ∂_ℓκ_ℓ = φ(loc κ₁) = 0 |

Por que P1 e P2 admitem κ_n = 0: o lado E de Sel_{F(n)} tem dimensão 1, mas o lado E^D tem dimensão 0 ou 2. Isso depende de a classe de Sel_{F(ℓ₁)}(E^D) ser localmente não nula em ℓ₂, e não sei calcular isso aqui. Pela heurística, κ_n ≠ 0 em cerca de 2/3 dos pares.

## O que conta como refutação (critérios fixados agora)

| # | critério | limite |
|---|---|---|
| R1 | κ_n ≠ 0 fora da reta prevista (P1, P2) | o conjunto de candidatos, depois de ≥ 15 primos q válidos, não contém a reta prevista e contém outra |
| R2 | κ_n ≠ 0 em P3 | W_n fora de 3E(F_q) em algum q válido |
| R3 | κ_ℓ ≠ 0 em C1 | W_ℓ fora de 3E(F_q) em algum q válido |
| R4 | classe fora de E(Q)/3 | conjunto de candidatos vazio com W ≠ 0 |

Um q só é "válido" se, nele, valem todos os testes internos:
- todos os pontos reduzidos estão em E(F_q);
- a ação de Galois forma ciclos únicos;
- vale a relação de traço Σ_t y_{n,(s,t)} = a_{ℓ₂}·y_{ℓ₁,s} (e a simétrica).

Se algum teste interno falhar de forma sistemática, a execução é **inválida** (defeito de cálculo): não aprova nem refuta.
