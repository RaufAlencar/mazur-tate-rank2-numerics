#!/bin/bash
# Mostra os geradores que kgen3.py usa para 389a1 e escreve-os na base do documento P1=(-1,1), P2=(0,0).
~/bsdvenv/bin/python - <<'EOF'
import cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9)
E = pari('ellinit([0,1,1,-2,0])')
G = pari.ellsaturation(E, pari.ellrank(E)[3], 1000)
print("geradores do kgen3:", G[0], G[1])
P1, P2 = pari('[-1,1]'), pari('[0,0]')
for name, g in (("G1", G[0]), ("G2", G[1])):
    sols = [(a, b) for a in range(-4, 5) for b in range(-4, 5)
            if pari.elladd(E, pari.ellmul(E, P1, a), pari.ellmul(E, P2, b)) == g]
    print(name, "= a*P1 + b*P2 com (a,b) =", sols)
EOF
