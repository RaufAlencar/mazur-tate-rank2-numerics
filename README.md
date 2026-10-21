# Numerical verification of the refined Mazur–Tate conjecture in ranks 2 and 3

Code, data and logs accompanying the note *Numerical verification of the refined Mazur–Tate conjecture in ranks 2 and 3, with an elementary formula for the Mazur–Tate pairing* (Rauf Alencar, October 2026). The current version of the note is **v0.2**: [`note/Alencar_MazurTate_rank0-3_v0.2.pdf`](note/Alencar_MazurTate_rank0-3_v0.2.pdf) (source text [`note/note_v0.2_rev131.md`](note/note_v0.2_rev131.md)). Version 0.1, as sent to two experts on 21 October 2026, is kept unchanged in [`note/Alencar_MazurTate_rank2_3_v0.1.pdf`](note/Alencar_MazurTate_rank2_3_v0.1.pdf) (git tag `v0.1`). v0.2 corrects a wrong statement about the sign in v0.1, explains the sign (PARI's modular-symbol sign and the Artin-map orientation), restricts the use of Ota's theorem to p ≥ 5, and adds ranks 0 and 1; see [`REPRODUCAO.md`](REPRODUCAO.md) §5–§7.

**Status: unreviewed numerical note.** Nothing here proves or claims progress on the Birch and Swinnerton-Dyer conjecture. The results are numerical checks of existing conjectures; the open points are listed in §7 of the note.

## What is checked

| Result | Scope | Where |
| --- | --- | --- |
| Leading-term identity M^{2r−1}·c_r(m) ≡ u_E·det⟨R_i,R_j⟩_m (mod (m−1)/2) for the Mazur–Tate theta element, with an elementary coordinate formula for the pairing | 10 curves in the note (ranks 2 and 3, about 700 non-trivial prime levels), no failures | note §2–§4, `code/research/` |
| Independent re-implementation from the text of the note, with u_E = (−2)ʳ·c∞/#tors² **predicted before running** | 7 further curves (563a1, 643a1, 664a1, 681c1, 707a1, 997b1, 11197a1 of rank 3), 237 non-trivial levels, no failures; the same code with u_E fixed by the formula also re-checked the 10 curves of the note (546 non-trivial levels, no failures); negative controls fail as they should | [`REPRODUCAO.md`](REPRODUCAO.md), `code/independent/`, `results/independent/` |
| Direction of Kolyvagin derived classes in rank 2: F_p·κ_ℓ = ker(E(Q)/p → E(F_ℓ)/p) | 20 cases, 3 curves, p ∈ {3, 5} | note §5, `kgen3.py` |
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
