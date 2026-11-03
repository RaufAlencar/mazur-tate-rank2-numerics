# Classes de Kolyvagin de duas pernas em posto 3: viabilidade e teste numérico

Sessão exploratória, 2026-11-03. Pasta: `BSD_Kolyvagin_2pernas`. Tudo aqui é para ser conferido pela sessão de registro antes de qualquer uso.

## Resumo

- **Viável para um caso**: 11642a1, K = Q(√−7), p = 3, n = 5·17 (grau [K[n]:K] = 108). A classe κ₈₅ ∈ E(Q)/3 ≅ F₃³ é **não nula** e gera exatamente a reta prevista, ⟨P₁+2P₂+2P₃⟩ = ker red₅ ∩ ker red₁₇, em 25/25 primos q, com todas as checagens internas passando. É o primeiro teste numérico da versão de duas pernas de que temos notícia.
- **Controles**: κ₅ = κ₁₇ = 0, como a teoria força, em 25/25 primos. O caso antigo de posto 2 (389a1, D=−7, ℓ=5) foi reproduzido com o kgen3 inalterado.
- **Replicação**: recálculo independente em 2600 dígitos, com outro método de reconhecimento, deu saída idêntica primo a primo.
- **Inviável nesta máquina**: 5077a1 (D=−163, n=11·23, grau 288) e o controle de zero forçado de duas pernas, 18097a1 (D=−163, n=11·29, grau 360). A altura dos pontos de Heegner exigiria ~9000 dígitos e ~10⁹–10¹⁰ termos da série modular (estimativa: ~100 h ou mais).
- Ressalva principal: é **um** par só. Ele confirma a direção, mas não a estatística (por exemplo, a fração de pares com κ_n = 0).

## 1. Viabilidade

### 1.1 Escolha dos casos (work/scan3b.py)

Varri as curvas de posto 3 com N < 20000 (tabelas de Cremona), exigindo:
- K com h = 1 e D < −4 (sem unidades extras);
- hipótese de Heegner;
- E[3] irredutível;
- E^D de posto analítico 0.

Para cada curva, procurei os menores primos de Kolyvagin ℓ (inerte em K, 3 | ℓ+1, 3 | a_ℓ; ℓ = 2 excluído). Só três casos têm grau (ℓ₁+1)(ℓ₂+1) ≤ 400:

| curva | D | (ℓ₁, ℓ₂) | grau | ker red_ℓ₁ ∩ ker red_ℓ₂ |
|---|---|---|---|---|
| 11642a1 | −7 | (5, 17) | 108 | dimensão 1, ⟨(1,2,2)⟩ |
| 5077a1 | −163 | (11, 23) | 288 | dimensão 1, ⟨(1,2,1)⟩ |
| 18097a1 | −163 | (11, 29) | 360 | dimensão 2: κ = 0 forçado |

Para 5077a1, os outros K com h = 1 dão grau ≥ 2016 ou têm E^D de posto 2. Uma primeira versão da varredura marcava erroneamente 11642a1 como tendo 3-isogenia; a checagem pelo fator linear de ψ₃ corrigiu isso.

Hipóteses verificadas (work/hyp.py):
- 3 ∤ c_Tam(E) e Ш(E/Q) = 1;
- E(K)[3] = 0 (ψ₃ sem raiz em K);
- E^D com posto analítico 0, torção trivial e Ш_an(E^D) = 1, 121, 100, nenhum divisível por 3.

Portanto Sel₃(E/K)⁺ = E(Q)/3 ≅ F₃³ e Sel₃(E/K)⁻ = 0.

### 1.2 Por que o método antigo não serve, e o que substituiu

O kgen3 reconhece o polinômio de x(y_n) sobre Q e roda polredbest no grau total; isso estoura perto do grau 200. Encontrei quatro obstáculos.

**(a) Série modular.** Com N = 11642 par e 2 decomposto em K, as classes do lado de 𝔭̄₂ só têm representante de Heegner com A ≈ 6322·N: seriam 24 milhões de termos por 100 dígitos. **Solução**: representantes via involuções de Atkin–Lehner W_Q, Q ∈ {2, 5821, 11642}, usando φ∘W_Q = w_Q·φ, sem constante, porque E(Q)_tors = 1. Os sinais w_Q foram conferidos numericamente (erro ~10⁻³⁰). O custo cai ~60×: pior classe com 414 mil termos por 100 dígitos (hreps.py).

**(b) Gerador pequeno das pernas.** O polredbest do polinômio de x(y₁₇), com denominadores de ~10^59, não termina nem com fatoração parcial. **Solução**:
- para Dℓ² ≡ 1 mod 8, polinômio de classe de Weber `polclass(Dℓ², 1)` (coeficientes de até 3 dígitos);
- para D = −163, γ₂ reduzido por polredbest parcial.

Cada raiz é casada à sua classe pelo valor de j: (f²⁴−16)³/f²⁴ = j(τ_c), ou γ₂³ = j, com erro ~10⁻¹⁵³. A ação de Galois mod q sai de resolventes dirigidos inteiros sobre O_K, como no kgen3.

**(c) Torre em vez do grau total.** x(y_n) e y(y_n) são escritos na base inteira {w_i(θ₁)·w_j(θ₂)} de K[5]·K[17] por um sistema linear 108×108 sobre os 108 conjugados complexos. O rótulo de cada conjugado é o par de projeções Pic(O₈₅) → Pic(O₅) × Pic(O₁₇), feitas sobre representantes coprimos a n.

**(d) Altura.** O denominador comum dos coeficientes de x(y_n) tem **717 dígitos**, e o de y(y_n) tem **1074**. O bestappr coordenada a coordenada precisaria de ~2×. **Solução**:
- denominador comum por LLL sobre 8 coordenadas (recc.gp), que exige ~1,15× os dígitos;
- para y, usar o d de x como fator conhecido.

Validei isso em dados sintéticos (d com 1100 dígitos e 1300 de precisão: 216/216 valores exatos). Com isso, 1200 dígitos bastam.

Outros dois defeitos meus, corrigidos durante a depuração:
- `pari.sqrt` e `pari.polroots` chamados como funções Python no cypari2 rodam em 64 bits, não na realprecision do GP;
- θ = B(x) pelo mapa inverso do polredbest perdia ~100 dígitos.

### 1.3 Por que 5077a1 e 18097a1 ficaram de fora

Na perna ℓ = 11 de 5077a1 (grau 12), os coeficientes de x já têm denominador 10^152, ou ~12,7 dígitos por conjugado. Em 11642a1 eram 3–4 dígitos por conjugado nas pernas e 6,6 no nível n. Extrapolação para o grau 288:
- d_x ~ 10^5800 e d_y ~ 10^8700;
- precisão necessária ≈ 9000 dígitos;
- 21,8 milhões de termos por 100 dígitos → ~2·10⁹ termos;
- tempo da ordem de 100 h.

Para 18097a1 (grau 360, 95 milhões de termos por 100 dígitos) é bem pior.

O que tornaria esses casos viáveis:
- uma máquina maior e dias de CPU; ou
- um método que evite reconhecer y_n: por exemplo, reduzir os pontos de Heegner diretamente mod q via curvas CM sobre F_q, o que exigiria a parametrização modular mod q, e não sei fazer isso para N = 5077; ou
- curvas de posto 3 com N maior (fora das tabelas que tenho, N < 20000) que tenham pares (ℓ₁, ℓ₂) pequenos com D pequeno.

## 2. Pré-registro

Arquivo `prereg/PREREGISTRO.md`, registrado em 2026-11-03T12:31:23-03:00. Seu sha256 está em `prereg/PREREGISTRO.sha256`: fbd139d0…0c3a3b. As previsões estão em `prereg/previsoes.json` (sha256 1f0920da…5e496b34) e foram geradas por `prereg/previsao.py`, independente do código de cálculo, antes de qualquer κ de posto 3 ser calculado.

| # | caso | previsão |
|---|---|---|
| P1 | 11642a1, n = 5·17 | κ_n ∈ ⟨(1,2,2)⟩: gera a reta ou é 0 |
| P2 | 5077a1, n = 11·23 | κ_n ∈ ⟨(1,2,1)⟩ ou 0 |
| P3 | 18097a1, n = 11·29 | κ_n = 0 (forçado: interseção de dimensão 2) |
| C1 | toda perna única | κ_ℓ = 0 (forçado: κ_ℓ ∈ Sel₃(E^D/Q) = 0) |

Por que P1 e P2 admitem 0: o lado E^D de Sel_{F(n)} tem dimensão 0 ou 2, e isso depende de uma classe ramificada que não sei calcular.

Critérios de refutação, fixados antes:
- **R1**: candidatos que, depois de ≥ 15 primos válidos, excluem a reta prevista e contêm outra;
- **R2**: P3 não nulo;
- **R3**: C1 não nulo;
- **R4**: conjunto de candidatos vazio com W ≠ 0.

Um primo q só é válido se todas as checagens internas passarem (pontos na curva, ciclos de Galois, relação de traço).

## 3. Resultados

### 3.1 Controle de posto 2 (validação do pipeline antigo)

Caso 389a1, D = −7, ℓ = 5, p = 3, kgen3 sem alterações (`ctrl/`, sha256 em `ctrl/sha_ctrl.txt`). Reta prevista = calculada (⟨P₂⟩ na base saturada atual do ellrank), em 20/20 primos. Isso reproduz uma linha da tabela da nota (lá rotulada ⟨P₁+P₂⟩ em outra base).

### 3.2 11642a1, n = 85 (P1 e C1)

Execução `runs/run_11642_p1200c.log`. Os pontos de nível n vêm do dump da execução `p1200b`, em 1200 dígitos; as pernas foram recalculadas. Hashes do código em `runs/run_11642_p1200c.sha256`.

- Reconhecimento:
  - pernas: d(x) = 10^19 e 10^75, d(y) = 10^27 e 10^112;
  - nível n: d_x = 10^717 (LLL direto), d_y = 10^1074 (LLL com o fator d_x).
- 25 primos q ≡ 1 mod 3, totalmente decompostos em K[85], entre 52879 e 453301. Em **todos**:
  - os 6 + 18 + 108 pontos reduzidos estão em E(F_q);
  - os ciclos de Galois são únicos;
  - valem as relações de traço Σ_t y_{n,(s,t)} = a₁₇·y_{5,s} e Σ_s y_{n,(s,t)} = a₅·y_{17,t};
  - Σ y_ℓ = 0 nas duas pernas.

  Nenhum primo foi descartado.
- **κ₈₅**: a classe é não nula em 22 dos 25 primos (W ∈ 3E(F_q) nos outros 3, compatível com localização nula). O conjunto de candidatos (a,b,c), módulo escalar, converge para {(1,2,2),(2,1,1)} no 7º primo e se mantém nos 18 seguintes.

  **Resultado: κ₈₅ gera ⟨P₁+2P₂+2P₃⟩, exatamente a reta prevista** (P₁ = (2,−2), P₂ = (3,−2), P₃ = (4,2), geradores saturados do ellrank). Uma classe errada ou aleatória esvaziaria a interseção em poucos primos.
- **κ₅ e κ₁₇**: W ∈ 3E(F_q) em 25/25 primos, nas duas pernas. **κ_ℓ = 0, como forçado.**

Nenhum critério de refutação foi atingido.

### 3.3 Replicação em 2600 dígitos

Execução `runs/run_11642_p2600.log` (sha256 ac05be60…c7a402; código em `runs/run_11642_p2600.sha256`):
- pontos de Heegner recalculados do zero em 2600 dígitos (6032 s para os 108 pontos);
- reconhecimento por `bestappr` direto, sem LLL;
- d_x = 10^717 e d_y = 10^1074, os mesmos valores obtidos via LLL.

Resultado: as 25 linhas por primo (classe W lida pelo emparelhamento de Tate, pareamentos dos geradores, todas as checagens) são **idênticas byte a byte** às da execução de 1200 dígitos. Mesmo resultado final: κ₈₅ gera ⟨(1,2,2)⟩ e κ₅ = κ₁₇ = 0.

### 3.4 Execuções inválidas (contam, mas não aprovam nem refutam)

| execução | o que aconteceu |
|---|---|
| precisão 120 | teste de depuração; nível n não reconhecido |
| precisão 400 | nível n saturou (d ≥ 10^199); 0 primos válidos de 181, porque os pontos não estavam na curva, como esperado. Gravada em kgen2p_results.jsonl com res "0" e contagens zero: **não é κ = 0, é execução inválida** |
| precisão 1200, sem LLL | saturou (d ≥ 10^599), abortada pela guarda |

Essas falhas mostram que a checagem "ponto na curva" discrimina reconhecimento errado.

## 4. Ressalvas

1. **Um único par.** A confirmação da direção é forte para esse par: 25 primos, todas as checagens, convergência no 7º primo. Mas não diz nada sobre a frequência de κ_n = 0 nem sobre outros pares, curvas ou primos p.
2. **Não testado**: o lado "não nulo ⇔ condição em E^D", que eu não sei calcular aqui, e o zero forçado de duas pernas (P3), que é inviável. O controle de zero forçado que existe é o de uma perna (C1).
3. **Reconhecimento numérico.** LLL e bestappr não são prova. A garantia vem de duas fontes:
   - as checagens mod q: 132 pontos por primo na curva e as relações de traço, em 25 primos;
   - a replicação em 2600 dígitos por outro caminho (bestappr direto), com saída idêntica.
4. **Hipóteses analíticas.** Ш(E) = 1 vem da tabela e Ш_an(E^D) foi calculado; ambos dependem de BSD/Kolyvagin para as partes 3-primárias. Kolyvagin já dá finitude para E^D de posto 0; não verifiquei Ш(E^D)[3] = 0 por descida.
5. **Teoria.** Que κ_n esteja em ker red_ℓ₁ ∩ ker red_ℓ₂ é quase elementar: segue dos axiomas de sistema de Kolyvagin (transversal em ℓ | n) mais a anulação de κ_ℓ. O conteúdo não trivial testado é κ_n ≠ 0 nesse par, isto é, o sistema de Kolyvagin de Heegner é não nulo mod 3 em ν = 2 para esta curva. É compatível com a conjectura de Kolyvagin e o resultado estrutural (Howard 2004, Mazur–Rubin 2004, W. Zhang). Não verifiquei literatura que já tenha calculado κ_{ℓ₁ℓ₂} explicitamente em posto 3.
6. Código novo, escrito nesta sessão, sem revisão independente. O pipeline antigo só foi revalidado no caso de posto 2.

## 5. Arquivos

| arquivo | o que é |
|---|---|
| `kgen2p.py`, `two.gp`, `hreps.py`, `recc.gp`, `divp.gp` | cálculo (uso: `python kgen2p.py ainvs D l1 l2 p PREC qmax maxq [KMAX]`; `FROMDUMP=1` retoma do dump) |
| `prereg/` | pré-registro, previsões e script de previsão |
| `runs/` | logs, hashes do código por execução, dumps numéricos das soluções |
| `ctrl/` | controle de posto 2 |
| `work/` | varreduras, testes e depuração |
| `kgen2p_results.jsonl` | uma linha por execução |
