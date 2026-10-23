# Calibração direta: símbolos [a/m]^+ reconstruídos das L-funções torcidas (κ = -c∞) contra msfromell, curvas pequenas.
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9); pari.default("realprecision", 40)
pari('read("theta_lfun.gp")')
for lab, a in [("389a1", "[0,1,1,-2,0]"), ("433a1", "[1,0,0,0,1]"), ("5077a1", "[0,0,1,-7,6]")]:
    pari(f"E = ellinit({a}); MS = msfromell(E, 1)")
    cinf = 2 if pari("E.disc") > 0 else 1
    for m in (11, 13, 17, 19, 23):
        pari(f"tmp = Svec_L(E, {m}); sy = symbols_from_S(tmp[1], {-cinf}, {m}, tmp[2])")
        ms = pari(f"vector({m}-1, a, mseval(MS[1], MS[2], [a/{m}, oo]))")
        same = pari("sy[1]") == ms
        inv = all(pari(f"sy[1][{a}] == {ms}[lift(Mod({a},{m})^-1)]") for a in range(1, m))
        print(lab, "c∞", cinf, "m", m, "| idênticos:", same, "| idênticos após a -> a^-1:", inv, "| erro", float(pari("sy[2]")), flush=True)
