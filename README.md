# Numerical verification of the refined Mazur–Tate conjecture in ranks 2 and 3

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23964110.svg)](https://doi.org/10.5281/zenodo.23964110)

Code, data and logs accompanying the note *Numerical verification of the refined Mazur–Tate conjecture in ranks 2 and 3, with an elementary formula for the Mazur–Tate pairing* (Rauf Alencar, October 2026). The current version of the note is **v0.7**: [`note/Alencar_MazurTate_rank0-3_v0.7.pdf`](note/Alencar_MazurTate_rank0-3_v0.7.pdf) (source text [`note/note_v0.7_rev161.md`](note/note_v0.7_rev161.md)); earlier versions are kept unchanged (tags `v0.2`, `v0.3`, `v0.4`, `v0.4.1`, `v0.5`, `v0.6`, `v0.6.1`). Version 0.1, as sent to two experts on 21 October 2026, is kept unchanged in [`note/Alencar_MazurTate_rank2_3_v0.1.pdf`](note/Alencar_MazurTate_rank2_3_v0.1.pdf) (git tag `v0.1`). v0.2 corrects a wrong statement about the sign in v0.1, explains the sign (PARI's modular-symbol sign and the Artin-map orientation), restricts the use of Ota's theorem to p ≥ 5, and adds ranks 0 and 1; v0.3 adds the rank-0 table (#Ш = 4 and 9), states the identity exactly in Ota's form under the stated Artin-map assumption, and explains why the odd-ℓ argument fails at ℓ = 2; v0.4 fixes a cut sentence in §4, adds rank 0 to the abstract (about 1700 non-trivial levels on 33 curves of ranks 0–3), and corrects the count of Ota-implication levels at p ≥ 5 to 11 (9 of them within the checkable hypotheses of Ota's theorem); v0.4.1 corrects one sentence of §7.4 (no rank-2 curve with #Ш ≠ 1 among the *optimal* curves up to 500,000, but non-optimal ones with #Ш = 4 exist); v0.5 adds an exploratory test of the #Ш factor in rank 2 (two curves with #Ш = 9 and two with #Ш = 4, θ_m from twisted L-values; one curve reproduced independently); v0.6 proves the integrality step at ℓ = 2 (Lemma in §2: an integral identity shows ℓ | δ ⇔ R ∉ E₀(Q_ℓ) for every prime ℓ), so the only remaining hypothesis of the pairing formula is the identification of the local splittings with those of Mazur–Tate (1983), and adds two exploratory paragraphs: a two-legged Kolyvagin class in rank 3 (§5) and statistics of det H mod p on 15,034 curves with a registered 3-isogeny effect (§6); see [`REPRODUCAO.md`](REPRODUCAO.md) §5–§9. **v0.7 corrects a statement of novelty in v0.6:** v0.6 said we had not found the pairing formula in print, but a formula for the same pairing, a denominator cross-ratio with an auxiliary point, goes back to Mazur–Tate (1987, p. 747) and was used by Portillo-Bobadilla (2019, eq. (18)) in his rank-1 verification. Our closed form agrees with it up to sign (checked on 389a1, 45 of 45 cases); what remains ours is the closed form without an auxiliary point, the lemma at ℓ = 2, and its use in ranks 2 and 3. The "elementary formula" in the title refers to that closed form.

**Website (PT/EN, interactive, for a general audience):** https://raufalencar.github.io/mazur-tate-rank2-numerics/ (source in [`docs/`](docs/); the lab values come from [`code/independent/site_data/dump.py`](code/independent/site_data/dump.py), output in `results/independent/site_levels_389a1_5077a1.json`; the cover art is AI-generated and only illustrative).

**Archived on Zenodo:** DOI [10.5281/zenodo.23964110](https://doi.org/10.5281/zenodo.23964110) (all versions; v0.6.1 is [10.5281/zenodo.23964111](https://doi.org/10.5281/zenodo.23964111), and each later version has its own DOI on the Zenodo page). Please cite the DOI.

**Status: unreviewed numerical note.** Nothing here proves or claims progress on the Birch and Swinnerton-Dyer conjecture. The results are numerical checks of existing conjectures; the open points are listed in §7 of the note.

## What is checked

| Result | Scope | Where |
| --- | --- | --- |
| Leading-term identity M^{2r−1}·c_r(m) ≡ u_E·det⟨R_i,R_j⟩_m (mod (m−1)/2) for the Mazur–Tate theta element, with an elementary coordinate formula for the pairing | 10 curves in the note (ranks 2 and 3, about 700 non-trivial prime levels), no failures | note §2–§4, `code/research/` |
| Independent re-implementation from the text of the note, with u_E = (−2)ʳ·c∞/#tors² **predicted before running** | 7 further curves (563a1, 643a1, 664a1, 681c1, 707a1, 997b1, 11197a1 of rank 3), 237 non-trivial levels, no failures; the same code with u_E fixed by the formula also re-checked the 10 curves of the note (546 non-trivial levels, no failures); negative controls fail as they should | [`REPRODUCAO.md`](REPRODUCAO.md), `code/independent/`, `results/independent/` |
| Direction of Kolyvagin derived classes in rank 2: F_p·κ_ℓ = ker(E(Q)/p → E(F_ℓ)/p) | 20 cases, 3 curves, p ∈ {3, 5} | note §5, `kgen3.py` |
| Two-legged Kolyvagin class in rank 3 (exploratory, added 3 Nov 2026): κ₈₅ for 11642a1, K = Q(√−7), p = 3, n = 5·17, predicted before computing to span ker red₅ ∩ ker red₁₇ = ⟨P₁+2P₂+2P₃⟩ | 1 case, 25/25 primes, replicated at 2600 digits; predicted line re-derived independently; κ₈₅ itself not recomputed by independent code | [`code/kolyvagin_two_legs/`](code/kolyvagin_two_legs/), [`REPRODUCAO.md`](REPRODUCAO.md) §9 |
| Strong BSD, p-adic BSD, 2-Selmer bounds on rank-2 curves (consistency checks) | see `code/research/README_research_session.md` | `bsd2.py`, `padic2.py`, `p6.py` |

## How the work was done

The author is not a professional mathematician. The research computations and the drafting of the note were carried out with an AI system (Claude, Anthropic) in one session; a second, separate session re-implemented the main identity from the text of the note alone and ran the blind test and controls. All numerical claims were checked by running code. Errors found along the way are recorded in the note and in the logs rather than removed.

## Reproducing

Requirements: Python 3.12, `cypari2` 2.2.4 (bundles PARI/GP 2.17.2), `numpy`. On Windows, use WSL (cypari2 has no native Windows wheel). `setup/setup_wsl.sh` builds the environment used here.

```bash
bash setup/setup_wsl.sh                      # creates ~/bsdvenv
cd code/independent
~/bsdvenv/bin/python check_appendix.py      # worked example of the note, 389a1 at m = 23
~/bsdvenv/bin/python blind_run.py 200 150   # blind test on the seven extra curves (~1 min)
~/bsdvenv/bin/python controls.py            # negative controls
```

The research scripts and their conventions are documented in [`code/research/README_research_session.md`](code/research/README_research_session.md) (sections 8–15) and [`code/research_s16_18/README_delivery.md`](code/research_s16_18/README_delivery.md) (sections 16–18 and the biextension computation; exact command lines in `stress_runs.sh`). Cremona's tables (`allbsd.*`) are downloaded from [github.com/JohnCremona/ecdata](https://github.com/JohnCremona/ecdata) and not redistributed here.

Caveats on the code: some early research scripts (`heegner.py`, `heegner_indiv.py`, `kolyvagin.py`, `koly2.py`, `koly3.py`, `kolyL.py`) have a known precision bug and are kept only as history. `data_derived/*.pkl` are Python pickles; load them only if you trust this repository.

## License

Code: MIT (see `LICENSE`). Note, data and logs: CC BY 4.0.

---

**Resumo em português.** Código, dados e registros de uma nota numérica sobre a conjectura refinada de Mazur–Tate em postos 2 e 3. O autor não é matemático profissional; o trabalho foi feito com auxílio de IA (Claude) e verificado rodando o código, incluindo uma reimplementação independente com previsão cega. Não é uma prova nem um avanço na conjectura de Birch e Swinnerton-Dyer.
