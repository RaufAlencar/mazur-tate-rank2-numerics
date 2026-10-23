# Reprodução independente do teste "#Ш > 1 em posto 2" (sessão sha-rank2).
# Passo 1: calibra κ comparando S_k das L-funções torcidas do PARI com S_k dos símbolos de msfromell (389a1, 433a1).
# Passo 2: em 53295337.a1 (#Ш = 9), alguns níveis primos: θ_m pelas L-funções torcidas, classe c_2 e det H com
#          mt_indep.gp, e teste da lei u = 4·c∞·#Ш/#tors² = 36 contra a alternativa sem #Ш (u = 4).
import cypari2, sys
pari = cypari2.Pari(); pari.allocatemem(3 * 10**9); pari.default("realprecision", 40)
pari('read("mt_indep.gp")'); pari('read("theta_lfun.gp")')

print("== calibração de κ (S_ms / S_L), deve ser constante = ±c∞ ==", flush=True)
for lab, a in [("389a1", "[0,1,1,-2,0]"), ("433a1", "[1,0,0,0,1]")]:
    pari(f"E = ellinit({a}); MS = msfromell(E, 1)")
    cinf = 2 if pari("E.disc") > 0 else 1
    for m in (11, 13, 17):
        pari(f"tmp = Svec_L(E, {m}); SL = tmp[1]; g = tmp[2]; SM = Svec_ms(E, MS, {m}, g)")
        ratios = [pari(f"SM[{k+1}]/SL[{k+1}]") for k in range(2, m - 2, 2)]
        print(lab, "c∞ =", cinf, "m =", m, "razões:", sorted(set(str(pari.round(pari.real(r) * 10**6) / 10**6) for r in ratios)), flush=True)

KAPPA_SIGN = int(sys.argv[1]) if len(sys.argv) > 1 else -1   # κ = KAPPA_SIGN · c∞ (decidido pela calibração)
print("\n== 53295337.a1 (#Ш = 9, tors 1, c_Tam 1) ==", flush=True)
pari("E = ellinit([1,-1,0,94,9]); P = ellsaturation(E, ellrank(E)[4], 1000)")
cinf = 2 if pari("E.disc") > 0 else 1; kappa = KAPPA_SIGN * cinf
ok = bad = triv = alt_fail = 0
for m in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    n = (m - 1) // 2
    pari(f"tmp = Svec_L(E, {m}); SL = tmp[1]; g = tmp[2]; sy = symbols_from_S(SL, {kappa}, {m}, g)")
    err = float(pari("sy[2]"))
    pari(f"c = theta_coefs_from_symbols(sy[1], {m}, g); cl = theta_class(c, {n}, 2)")
    if not int(pari("cl[1]")):
        print("m", m, "θ fora de I² ou não inteiro:", pari("cl"), "erro de arredondamento", err); continue
    cr = int(pari("cl[2]")); M = int(pari(f"ellcard(E, {m})"))   # c_Tam = 1
    dH = int(pari(f"pair_det(E, P, {M}, {m}, g)[1]"))
    lhs = (4 * pow(M, 3, n) * cr) % n
    law, alt = (36 * dH) % n, (4 * dH) % n
    tag = "trivial" if lhs == 0 and law == 0 else ("OK" if lhs == law else "FALHA")
    ok += tag == "OK"; bad += tag == "FALHA"; triv += tag == "trivial"
    if law != alt and lhs != alt: alt_fail += 1
    print(f"m {m:3d} n {n:3d} | erro {err:.1e} | c2 {cr:3d} detH {dH:3d} | lhs {lhs:3d} lei(36) {law:3d} sem Ш(4) {alt:3d} | {tag}", flush=True)
print(f"\nlei com #Ш: {ok} OK, {triv} triviais, {bad} falhas | alternativa sem #Ш falha em {alt_fail} níveis discriminantes")
