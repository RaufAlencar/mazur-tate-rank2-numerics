import cypari2
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9); pari.default("realprecision", 40)
pari('read("theta_lfun.gp")')
for lab, a in [("389a1", "[0,1,1,-2,0]"), ("5077a1", "[0,0,1,-7,6]")]:
    pari(f"E = ellinit({a}); MS = msfromell(E, 1)")
    for m in (13, 17):
        pari(f"tmp = Svec_L(E, {m}); sy = symbols_from_S(tmp[1], -1, {m}, tmp[2])")
        ok = all(pari(f"sy[1][{x}] == mseval(MS[1], MS[2], [lift(Mod({x},{m})^-1)/{m}, oo])") for x in range(1, m))
        print(lab, "m", m, "com κ = -1 (não -c∞) e a -> a^-1: idênticos =", ok, flush=True)
