# Confere a implementação independente contra o exemplo do apêndice da nota (389a1, m = 23).
import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari('read("mt_indep.gp")')
pari('E = ellinit([0,1,1,-2,0]); MS = msinit_curve(E); m = 23; g = znprimroot(m)')
print("g =", pari('g'))
print("[a/23]+ :", pari('vector(22, a, mseval(MS[1], MS[2], [a/23, oo]))'))
print("c_t     :", pari('theta_coefs(MS, 23, g)'), " (nota: 0,-2,-2,2,0,0,2,-2,-2,0,4)")
print("c_2(23) :", pari('theta_class(theta_coefs(MS, 23, g), 11, 2)'), " (nota: 6)")
print("det, H  :", pari('pair_det(E, [[-1,1],[0,0]], 28, 23, g)'), " (nota: det 10, H = (9,7;7,9) mod 11)")
print("check   :", pari('mt_check(E, MS, [[-1,1],[0,0]], 1, 23)'), " (nota: lhs 3, u_E*det = 8*10 = 3 mod 11)")
