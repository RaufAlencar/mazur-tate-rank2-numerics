# BSD — Mapa da Obstrução (posto ≥ 2): código, dados e resultados

Pacote de tudo que foi produzido na sessão de 20–21 out 2026 para o documento
"BSD — Mapa da Obstrução (posto ≥ 2)". Curva principal: **389a1**, `[0,1,1,-2,0]`,
posto 2, geradores P1 = (−1, 1), P2 = (0, 0). Tudo em PARI/GP via Python.

## Ambiente

- Python 3.11+
- **cypari2 2.2.4** (`pip install cypari2`), que embute **PARI/GP 2.17.2**
- numpy (só para `heegner.py`/`heegner_indiv.py`)
- Dados de Cremona: `allbsd.*` e `allcurves.*` de https://github.com/JohnCremona/ecdata
  (pasta `allbsd/` e `allcurves/`). Baixar para `data/` antes de rodar bsd2.py, p6.py, padic2.py, kgen*.py.
  Formato das linhas de `allbsd`: `N iso num ainvs rank tors cp omega L reg sha`.
- Memória: os scripts de Kolyvagin pedem até 3–8 GB de pilha PARI (`pari.allocatemem`). Em máquina com 7 GB, graus
  de corpo acima de ~120 estouram (`polredbest`); ajuste o `allocatemem` se tiver mais RAM.

### Armadilha de precisão (importante)

No cypari2, **chamadas de método** (`pari.sqrt(x)`, `pari.ellztopoint(...)`) ignoram `pari.default('realprecision', n)`
e rodam em ~19 dígitos. Só **strings GP** (`pari('sqrt(-5)')`, `pari('read("heeg.gp")')`) respeitam a precisão.
Por isso todo cálculo de alta precisão está em arquivos `.gp` chamados por `pari('read(...)')`.
Os scripts antigos `heegner.py`, `heegner_indiv.py`, `kolyvagin.py`, `koly2.py`, `koly3.py`, `kolyL.py` têm esse bug
(resultados só confiáveis a ~15 dígitos); ficam no pacote por histórico. Os resultados válidos vêm de `heeg.gp` + scripts posteriores.

## Scripts → seções do documento

| Seção do documento | Script(s) | O que faz | Saída |
| --- | --- | --- | --- |
| 8 (Exp. 1) BSD forte, 47 curvas posto 2 | `bsd2.py N` | r_an, geradores saturados, regulador, Ш analítico vs Cremona | stdout |
| 8 (Exp. 1b) varredura Ш em 279 056 curvas | one-liners awk sobre `data/allbsd.*` (ver seção 8 e 9 do doc) | contagem de Ш≠1 por posto | — |
| 8 (Exp. 2) BSD p-ádica | `padic.py N`, `padic2.py N` | `ellpadicbsd`, `ellpadicregulator`, valuações | `padic2.log`, `padic2_1500.json` |
| 8 (Exp. 3) Pontos de Heegner (baixa precisão) | `heegner.py label ainvs k`, `heegner_indiv.py` | P_K e pontos individuais (≈15 dígitos) | stdout |
| 9 Goldfeld, fator de Euler p-ádico | one-liners na sessão (ver seção 9) | — | — |
| 10 (P1 parte 1) primos de Kolyvagin e retas previstas | one-liner na sessão; reimplementado em `kgen3.py` (bloco PREDICTION) | ker(E(Q)/p→E(F_ℓ)/p), divisibilidade do gerador do twist | stdout |
| 10 (P2) reta de Darmon–Rotger | `p2.py N` | N(E,p)=h_p(P0,P0)/R_p com `ellpadicheightmatrix` (convenção H = A − s2·B) | stdout |
| 10 (P6) 2-Selmer | `p6.py Nmax [Nmin]` | `ellrank` limite superior = 2 em todas as curvas de posto 2 | `p6.log`, `p6b.log`, `p6full.log` |
| 11 κ₁₇ para 389a1, K=Q(√−5) (h=2) | `heeg.gp` (run2) → `koly_full.gp` + `twist.gp` + `divtest.gp`; dados em `data_derived/heeg_17*.pkl`, `resolvents.pkl`, `directed.pkl`, `zeta.pkl`, `ycoord.pkl`; `koly_consist.gp`, `koly_solve.gp`, `wcheck.gp` são diagnósticos | soma sobre os dois cosets com orientação por resolventes sobre K | stdout (ver doc) |
| 11 controle 43a1 (κ=0) | `pipeline.py` | pipeline h(K)=1 | stdout |
| 12 bateria 12 casos (389a1, p=3) | `kgen2.py D l PREC qmax maxq` + `battery.sh`, `battery2.sh`; `kgen.py` é a versão h=2 (estourou memória em ℓ=53) | classe κ_ℓ identificada em (a,b,c), κ = aP1+bP2+cQ | `results/result_D*_l*.json`, `logs/battery*.log`, `data_derived/kgen_*.pkl` |
| 12/16 bateria 2 (p=5, curvas 571b1 e 709a1) | `kgen3.py ainvs D l p PREC qmax maxq` + `battery3.py` | pipeline geral (qualquer curva, qualquer p); imprime PREDICTION e FINAL | `results/kgen3_results.jsonl`, `logs/battery3.log` |
| 13 (P3) ordem de θ_m de Mazur–Tate | `theta.gp` (`thetaord`) | ord_x de θ_m em F_p[x]/(x^{p^e}) | stdout |
| 15 termo líder de Mazur–Tate | `mt_lead.gp` (`thetacoefs`, `mtreg`) + `mt_batch.py ainvs mmax` | c₂(m) vs −#E(F_m)·det⟨P_i,P_j⟩_m com ⟨R,R⟩_m = [t(R)/d(R)] | `logs/mt_*.log` |

## Como rodar (exemplos)

```bash
pip install cypari2 numpy
mkdir -p data && cd data && curl -O https://raw.githubusercontent.com/JohnCremona/ecdata/master/allbsd/allbsd.00000-09999 && cd ..

# BSD forte nas curvas de posto 2 com N ≤ 1500
python3 bsd2.py 1500

# 2-Selmer em todas as curvas de posto 2 com N ≤ 6000
python3 p6.py 6000

# Classe de Kolyvagin: 389a1, K=Q(sqrt -7), l=5, p=3, 120 dígitos, primos q ≤ 300000, 20 primos
python3 kgen3.py "[0,1,1,-2,0]" -7 5 3 120 300000 20
# saída: PREDICTION (reta prevista) e FINAL (candidatos (a,b,c) a menos de escalar)

# Termo líder de Mazur–Tate em 389a1 para m ≡ 1 mod 9, m ≤ 1200
python3 mt_batch.py "[0,1,1,-2,0]" 1200
```

Convenções que custaram horas (ver seções 11 e 15 do documento):
- `elltatepairing(E, P, Q, m)`: **P** é o ponto de m-torção; o valor vem em F*/(F*)^m e precisa ser elevado a (q−1)/m.
- `ellpadicheightmatrix` devolve `[A, B]`; a matriz da altura ciclotômica é `A − ellpadics2(E,p,n)·B`.
- `ellglobalred(E)[3]` é a fatoração do condutor, **não** os números de Tamagawa (esses estão em `[4][,4]`).
- Para h(K) > 1 a classe de Kolyvagin é a soma das derivadas sobre representantes de Gal(K[n]/K[1]) em Gal(K[n]/K),
  com orientação relativa fixada por resolventes direcionadas sobre O_K (`koly_full.gp`).

## Resultados principais (todos no documento)

- BSD forte: 47/47 curvas de posto 2 (N ≤ 1500), Ш = 1; nenhuma curva de posto 2 com Ш ≠ 1 em 279 056 (N ≤ 500 000).
- 2-Selmer: limite superior = 2 em 348 672 curvas de posto 2 (N ≤ 500 000); Ш[2] = 0 em todas.
- Classes de Kolyvagin: direção = ker(E(Q)/p → E(F_ℓ)/p) e critério de anulamento confirmados em 17/17 casos
  (4 corpos, p = 3 e 5, curvas 389a1 e 571b1) — enunciado é Prop. 0.2 do rascunho de Stein; verificação numérica da direção parece inédita.
- Mazur–Tate: ord_I θ_m ≥ posto em 75/75 (postos 0–3); termo líder c₂(m) ≡ −#E(F_m)·det⟨P_i,P_j⟩_m (mod 3) em 20/20 primos de 389a1.

Formatos: `.pkl` são pickles Python com strings GP (polinômios g, A, ζ, resolventes) para reproduzir sem recomputar os
pontos de Heegner; `.json/.jsonl` são os vereditos finais por caso; `.log` são as saídas brutas.
