# Quais tipos de Kodaira em p = 2 aparecem entre as curvas de posto positivo com N <= 400 (a amostra do gsweep da nota)?
import cypari2, collections, os
pari = cypari2.Pari(); pari.allocatemem(10**9)
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "projeto", "scripts", "data", "allbsd.00000-09999")
types = collections.Counter(); ncurves = 0
for line in open(DATA):
    f = line.split()
    N, ainvs, rank = int(f[0]), f[3], int(f[4])
    if N > 400 or rank == 0: continue
    ncurves += 1
    if N % 2: types["boa redução em 2"] += 1; continue
    k = int(pari(f"elllocalred(ellinit({ainvs}), 2)[2]"))
    name = {1: "I0", 2: "II", 3: "III", 4: "IV", -1: "I0*", -2: "II*", -3: "III*", -4: "IV*"}.get(k)
    if name is None: name = f"I{k-4}" if k > 4 else f"I{-k-4}*"
    types[name] += 1
print("curvas de posto > 0 com N <= 400 (todas as da tabela, não só ótimas):", ncurves)
for t, c in sorted(types.items(), key=lambda x: -x[1]): print(f"  {t:18s} {c}")
