# Entrega — seções 16–18 (identidade de Mazur–Tate, níveis compostos, bateria 2 de Kolyvagin) e §2 da nota (biextensão)

Data: 2026-10-21. Ambiente: Python 3.11, cypari2 2.2.4 (PARI/GP 2.17.2), 7 GB RAM. Todos os scripts rodam a partir da pasta `scripts/`
(`cd scripts && python3 ...`); os `.gp` são lidos pelos `.py` via `pari('read(...)')`. Nenhum script precisa das tabelas de Cremona,
exceto `gsweep.py` (`data/allcurves.00000-09999` de https://github.com/JohnCremona/ecdata) e `kgen3.py` (ver README do zip anterior).

## Scripts → o que fazem → log

| Script | Função | Log / resultado |
| --- | --- | --- |
| `mt_full.gp` + `mt_full_run.py ainvs mmax U` | identidade prime-level posto 2 (`4·M³·c₂ ≡ U·det H mod n`), conta OK/trivial/FAIL | `logs/mt_full_389_1500.log` (389a1, 221 OK / 15 triviais / 0 falhas), `logs/mt_full_others.log` |
| `mt_scan_u.py ainvs mmax` | varre a constante U | — |
| `mt_rank.gp` + `mt_rank_run.py ainvs mmax r` | identidade para posto r qualquer (`2^r M^{2r−1} c_r ≡ u·det H`), imprime o melhor u | `logs/stress1_*.log`, `logs/stress2.log` (655a1, 1132a1, 916c1, 571b1, 1088j1, 5077a1) |
| `stress_runs.sh` | as linhas de comando exatas de todas as curvas da tabela da nota | — |
| `bmc_levels.py` | det H mod 3 e identidade completa nos níveis de Bley–Macías Castillo (389a1: 7,13,43; 433a1: 7,13,31,37; 446d1: 19) e a linha 433a1 | saída reproduzida em `RESULTS.md` |
| `mt_ota.gp` + `mt_ota_run.py ainvs mmax gmin` | fórmula de inclusão–exclusão de Ota em nível composto S=m₁m₂ (coeficientes x₁², x₁x₂, x₂²) | `logs/mt_ota_389.log` (28 pares), `logs/mt_ota_433.log` (em andamento no momento da entrega: 26/28 pares) |
| `mt_comp.gp` | versão exploratória anterior (sem inclusão–exclusão) — histórico | — |
| `biext.gp`, `chkg.gp`, `chkg2.gp`, `chord.gp` | computação na biextensão (§2 da nota): E·c/A(2R) vs (t/d)²; identidades inteiras δ, c_L; fórmula da secante | saída em `RESULTS.md` |
| `gsweep.gp` + `gsweep.py Nmax` | varredura δ = c_L = 1 em pontos de E₀, curvas de posto ≥1, condutor ≤ Nmax | `logs/gsweep.log` (453 curvas, 3631 pontos, 0 exceções) |
| `theta.gp`, `mt_lead.gp`, `mt_batch.py` | ordem de θ_m e termo líder mod 3 (seções 13/15, já entregues antes) | — |
| `kgen3.py`, `heeg.gp`, `battery3.py` | bateria 2 de Kolyvagin (p=5; 571b1, 709a1) | `logs/battery3.log`, `results/kgen3_results.jsonl` |
| `rank4.log` | tentativa 234446a1: OOM em `msfromell` (7 GB) | `logs/rank4.log` |

## Convenções que importam

- `pairdet`/`Bmat`: H_ii = 2·h(R_i), H_ij = h(R_i+R_j) − h(R_i) − h(R_j), com h(R) = log_g[(t/d) unit part mod m], t = −x/y, d² = denominador de x. H é a matriz do emparelhamento (polarização da forma quadrática h). det H é o que a nota chama det H.
- u impresso por `mt_rank_run.py` é o u_E da nota: `2^r M^{2r−1} c_r ≡ u·det H (mod n)`, M = c_Tam·#E(F_m).
- `mt_ota_run.py` escala θ por Den = 2 quando c_∞ = 1 (para ter inteiros); por isso em 433a1 a constante composta impressa é 8 = 2·u_E, e em 389a1 é 8 = u_E. "u fits" lista os u com `4M³a₁₁ + u·D₁₂ ≡ 0 mod g`; a previsão é u = Den·u_E (sinal trocado em relação ao nível primo).
- Base dos geradores: sempre `ellsaturation(ellrank(E)[4],1000)`; em 389a1 isso dá (0,−1), (1,0) (= −P₂, P₁+P₂ na base (−1,1),(0,0)). Os rótulos de reta da §5 da nota para 571b1 e 709a1 estão nessa base saturada.
