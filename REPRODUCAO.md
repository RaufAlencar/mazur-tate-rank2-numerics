# Reprodução independente — 21/10/2026

Ambiente: WSL Ubuntu nesta máquina, Python 3.12.15 (via uv), cypari2 2.2.4 / PARI/GP 2.17.2 (mesma versão da sessão original).
Preparação: `setup_wsl.sh`. Scripts rodados a partir de `projeto/scripts` (cópia em `~/bsdwork` para não alterar o pacote).

| Teste | Seção | Script | Resultado original | Reproduzido aqui | Log |
| --- | --- | --- | --- | --- | --- |
| BSD forte, 47 curvas de posto 2, N ≤ 1500 | 8 | `run_bsd2.sh` → `bsd2.py 1500` | Ш = 1 em 47/47 | **Ш = 1 em 47/47**, 5 s; Reg(389a1) = 0,15246018 | `repro_bsd2_1500.log` |
| Classe de Kolyvagin, 389a1, K = Q(√−7), ℓ = 5, p = 3 | 12 | `run_kgen3.sh` → `kgen3.py` | reta prevista ⟨P₁+P₂⟩, κ = −(P₁+P₂) | **prevista = calculada**, reta ⟨P₁+P₂⟩, 20 primos, traço ok em todos, < 1 s | `repro_kgen3_389a1_D-7_l5.log` |

## Ponto de atenção encontrado

`kgen3.py` não usa os geradores do documento, P₁ = (−1, 1) e P₂ = (0, 0). Ele usa os que `ellrank` + `ellsaturation` devolvem: G₁ = (0, −1) = −P₂ e G₂ = (1, 0) = P₁ + P₂ (verificado por `check_basis.sh`).
Por isso a saída diz "(0, 1, 0)", que significa G₂ = P₁ + P₂ na base do documento. O resultado bate, mas a notação engana.

Antes de publicar: conferir em qual base estão escritas as linhas da seção 16 do `resumo_resultados.csv` (casos rodados com `kgen3.py`, p = 5 e 571b1) e fixar uma base única em todo o texto.

---

# Reprodução independente da identidade de Mazur–Tate (nota rev 80) — 21/10/2026

Implementação escrita **do zero, só a partir do texto da nota** (§1–§3), sem usar os scripts da sessão de pesquisa: `repro_indep/mt_indep.gp`.

## 1. Apêndice da nota (389a1, m = 23)
Reproduzido exatamente: símbolos [a/23]⁺, coeficientes c_t, c₂(23) = 6, H = (20, −15; −15, 20), det H ≡ 10, lado esquerdo 2²·28³·6 ≡ 3 ≡ 8·10 (mod 11). Script: `check_appendix.py`.

## 2. Teste cego (previsão registrada antes de rodar)
Previsões u_E = (−2)ʳ·c∞/#tors² gravadas em `blind_predictions.json` às 11:46, para 7 curvas **fora da tabela da nota**; depois rodado `blind_run.py` (7 ≤ m ≤ 200 em posto 2; m ≤ 150 em posto 3).

| Curva | r | c∞ | c_Tam | u_E previsto | não triviais OK | falhas | θ_m ∉ Iʳ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 563a1 | 2 | 1 | 1 | 4 | 34 | 0 | 0 |
| 643a1 | 2 | 1 | 1 | 4 | 36 | 0 | 0 |
| 664a1 | 2 | 1 | 4 | 4 | 36 | 0 | 0 |
| 681c1 | 2 | 1 | 2 | 4 | 37 | 0 | 0 |
| 707a1 | 2 | 2 | 2 | 8 | 33 | 0 | 0 |
| 997b1 | 2 | 2 | 1 | 8 | 35 | 0 | 0 |
| 11197a1 | 3 | 2 | 1 | −16 | 26 | 0 | 0 |

**Total: 237 verificações não triviais, 0 falhas**, incluindo uma segunda curva de posto 3. Logs: `blind_run.log`, `blind_results.json`.

## 3. Controles negativos (`controls.py`, `controls.log`)
O teste falha quando deve: constante errada (563a1 com u = −4, 8, 2: 0–2 acertos, 32–34 falhas); sem c_Tam (664a1: 2 acertos, 37 falhas); c_Tam = 2 em vez de 4 (1/36); pontos de 563a1 com θ de 643a1 (1/38); posto 3 com sinal trocado u = +16 (0/26).

## Conclusão
A identidade da nota, incluindo a regra para u_E e o papel de c_Tam, foi confirmada por uma implementação independente e por uma previsão cega em curvas novas. Ressalvas que continuam valendo (nota §7): sinal (−1)ʳ vs (−1)^{ν(m)} em aberto; lema δ = 1 em ℓ = 2 não provado; fórmula do emparelhamento derivada assumindo a descrição de Néron dos splittings. Reprodução no mesmo PARI/GP 2.17.2 — independente no código, não na biblioteca.

## 4. As 10 curvas da tabela da nota, com a implementação independente (`table_run.py`, `table_run.log`)
Constante FIXA pela fórmula u_E = (−2)ʳ·c∞/#tors² (para 1088j1, testado tors²·lhs ≡ (−2)ʳ·c∞·det, sem dividir mod n). Base: a da nota para 389a1; saturada do `ellrank` nas demais.

| Curva | r | c∞ | tors | c_Tam | u_E | m ≤ | não triviais OK | falhas | θ_m ∉ Iʳ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 389a1 | 2 | 2 | 1 | 1 | 8 | 600 | 90 | 0 | 0 |
| 709a1 | 2 | 2 | 1 | 1 | 8 | 400 | 65 | 0 | 0 |
| 446d1 | 2 | 2 | 1 | 2 | 8 | 400 | 62 | 0 | 0 |
| 571b1 | 2 | 1 | 1 | 1 | 4 | 330 | 55 | 0 | 0 |
| 5077a1 | 3 | 2 | 1 | 1 | −16 | 330 | 47 | 0 | 0 |
| 655a1 | 2 | 1 | 1 | 2 | 4 | 330 | 55 | 0 | 0 |
| 1132a1 | 2 | 1 | 1 | 3 | 4 | 330 | 43 | 0 | 0 |
| 916c1 | 2 | 2 | 1 | 3 | 8 | 330 | 45 | 0 | 0 |
| 1088j1 | 2 | 2 | 2 | 4 | 8/4 = 2 | 300 | 47 | 0 | 0 |
| 433a1 | 2 | 1 | 1 | 1 | 4 | 200 | 37 | 0 | 0 |

**Total da implementação independente: 17 curvas (10 da nota + 7 cegas), 783 níveis não triviais, 0 falhas.**

## 5. Sinal em posto 1 (depois do envio da v0.1) — `sign_rank1.py`, `sign_rank1.log`
Motivação: um revisor externo notou que a frase da nota (§4(b), §7.2) "um fator (−1)^{r(r−1)/2} resolveria o sinal" está errada (em r = 3 ela trocaria o sinal certo) e sugeriu (−1)^{r(r+1)/2}. Com só r = 2, 3, ambas as hipóteses (−1)^{r+1} e (−1)^{r(r+1)/2} cabem; elas divergem em r = 1.
Hipóteses registradas antes de rodar em `sign_predictions_rank1.json` (13:29): H_A (lei da nota) u_E = −2c∞; H_B (revisor) u_E = +2c∞.

| Curva | c∞ | c_Tam | níveis não triviais | acertos u = −2c∞ | acertos u = +2c∞ |
| --- | --- | --- | --- | --- | --- |
| 37a1 | 2 | 1 | 67 | 67 | 2 |
| 43a1 | 1 | 1 | 66 | 66 | 2 |
| 53a1 | 1 | 1 | 68 | 68 | 1 |
| 57a1 | 1 | 2 | 69 | 69 | 1 |
| 58a1 | 1 | 2 | 67 | 67 | 1 |
| 61a1 | 1 | 1 | 69 | 69 | 1 |
| 77a1 | 1 | 2 | 62 | 62 | 1 |
| 79a1 | 2 | 1 | 65 | 65 | 3 |

**H_A confirmada (533/533), H_B refutada.** A lei u_E = (−2)ʳc∞/#tors² vale em r = 1, 2, 3; o fator entre a convenção da nota e o −1 fixo de Ota é (−1)^{r+1}, não (−1)^{r(r−1)/2} (nota) nem (−1)^{r(r+1)/2} (revisor).

## 6. Posto 0 e a origem do sinal (`sign_rank0.py`, `ms_sign.py`)
Sugerido pelo revisor externo. Hipóteses registradas às 13:54 em `sign_predictions_rank0.json`: H_A (lei da nota) u·tors² = +c∞·Ш; H_B (relação de Hecke Σ[a/m]⁺ = (a_m − 2)[0]⁺ + BSD, supondo [0]⁺ = +L/ω₁) u·tors² = −c∞·Ш. Em posto 0 a identidade é M⁻¹·a₀ ≡ u (mod n), a₀ = Σ_a [a/m]⁺ (racional quando há torção; tratado com aritmética exata — uma primeira rodada com um atalho errado para a₀ racional foi descartada).

| Curva | c∞ | c_Tam | tors | Ш | níveis | H_A (+) | H_B (−) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 11a1 | 1 | 5 | 5 | 1 | 29 | 29 | 0 |
| 14a1 | 1 | 6 | 6 | 1 | 18 | 18 | 0 |
| 15a1 | 2 | 8 | 8 | 1 | 30 | 30 | 0 |
| 32a1 | 1 | 4 | 4 | 1 | 39 | 39 | 0 |
| 67a1 | 1 | 1 | 1 | 1 | 39 | 39 | 0 |
| 75a1 | 1 | 1 | 1 | 1 | 27 | 27 | 0 |
| 571a1 | 1 | 1 | 1 | **4** | 37 | 37 | 0 |
| 681b1 | 2 | 4 | 4 | **9** | 30 | 30 | 0 |

**H_A: 249/249.** A lei u_E = (−2)ʳ·c∞·#Ш/#tors² vale agora em r = 0, 1, 2, 3, e em posto 0 o fator #Ш (4 e 9) aparece exatamente — algo que os postos ≥ 2 não podiam testar, pois lá #Ш = 1 em todas as curvas.

**Por que H_B errou — e o que isso resolve.** `ms_sign.py`: o símbolo do PARI dá `mseval(…, [0, oo])` = **−L(E,1)/ω₁** (67a1: −1 contra L/ω₁ = 1; 571a1: −4 contra 4; 11a1: −1/5 contra 1/5). Logo o θ_m calculado com `mseval([a/m, oo])` (nos scripts da pesquisa e nos meus) é **−θ_m na convenção de Mazur–Tate** que a nota declara em §1 ([0]⁺ = +L/ω₁). Consequência: na convenção de Mazur–Tate a constante é −(−2)ʳ·c∞Ш/tors², isto é, sinal −(−1)ʳ; invertendo a orientação do mapa de Artin nas unidades locais (que, pela §4(b) da nota, multiplica o det por (−1)ʳ), o sinal fica **−1 para todo r**, que é exatamente o (−1)^{ν(m)} de Ota em nível primo. O "sinal em aberto" reduz-se então a uma convenção verificada (o sinal do símbolo do PARI) e a uma escolha binária de orientação do Artin, consistente em r = 0, 1, 2, 3. A segunda parte ainda precisa ser conferida no texto de Ota antes de virar afirmação.

## 7. Implicação de Ota em p = 3, 5, 7, 11, 13 (`ota_implication.py`, `ota_implication.log`)
Sugerido pelo revisor externo depois de constatarmos que Ota (2018) prova seus resultados só com 2 e 3 invertidos em R (Thm. 1.2, condição (i): p | 6N·|E(F_p)|·∏ c_ℓ ⇒ p invertível; também exige curva sem CM, ρ̄ sobrejetiva e E(F_ℓ)[p] cíclico; §3 e §5 fixam p ≥ 5). A nota (§4, "The vanishing levels") usava a implicação "termo líder ≢ 0 mod p ⇒ p ∤ J_m" em p = 3, fora dessas hipóteses.
Teste: para primos m com p | n = (m−1)/2 e p | J_m (J_m = #E(F_m)/#⟨P_i mod m⟩, curvas com c_Tam = 1 e torção trivial), verificar c_r(m) ≡ 0 mod p. Curvas 389a1, 709a1, 433a1, 563a1, 643a1, 997b1, 571b1 (m ≤ 3000) e 5077a1 (m ≤ 1500).

| p | níveis com p ∣ J_m | c_r ≡ 0 mod p | violações |
| --- | --- | --- | --- |
| 3 (fora do teorema) | 62 | 62 | 0 |
| 5 | 5 | 5 | 0 |
| 7 | 4 | 4 | 0 |
| 11, 13 | 0 | — | — |

Nenhuma violação. Ressalvas: em p ≥ 5 os casos que testam a implicação são poucos (9) por natureza; as demais hipóteses de Ota (ρ̄ sobrejetiva, E(F_m)[p] cíclico) não foram conferidas caso a caso. Observação: a verificação de θ_m ∈ Iʳ feita nas seções anteriores é sobre Z, inclusive nas partes 2 e 3, que o Thm. 1.2 de Ota não cobre.
