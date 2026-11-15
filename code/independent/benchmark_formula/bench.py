# Pré-registro B (PREREGISTRO_B.md): forma fechada (A) x fórmula com ponto auxiliar (B), mesmos pontos.
import cypari2, json, os, sys, time, statistics
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
for _d in [
 "dd(T)=sqrtint(denominator(T[1]))",
 "A_pair(E,R,S,m)={my(lam,c,den,L,RS);if(R==S,return(Mod(numerator(R[1]),m)^-1));lam=(S[2]-R[2])/(S[1]-R[1]);c=R[2]-lam*R[1];den=lcm(denominator(lam),denominator(c));L=[den,-den*lam,-den*c];L=L/content(L);RS=elladd(E,R,S);Mod(L[3],m)/Mod(numerator(RS[1]),m)}",
 "A_wrong(E,R,m)=Mod(numerator(R[1]),m)",
 "B_pair(E,R,S,Pp,m)={my(a=dd(elladd(E,Pp,R)),b=dd(elladd(E,Pp,S)),c=dd(Pp),e=dd(elladd(E,Pp,elladd(E,R,S))));Mod(a*b,m)/Mod(c*e,m)}",
 "matA(E,RR,m)={my(r=#RR);matrix(r,r,i,j,A_pair(E,RR[i],RR[j],m))}",
 "matB(E,RR,Pp,m)={my(r=#RR);matrix(r,r,i,j,B_pair(E,RR[i],RR[j],Pp,m))}",
 "eqpm(x,y)=(x==y)||(x==-y)"]:
    pari(_d)
CURVAS = [("389a1", "[0,1,1,-2,0]", "[[-1,1],[0,0]]"), ("433a1", "[1,0,0,0,1]", None),
          ("709a1", "[0,-1,1,-2,0]", None), ("5077a1", "[0,0,1,-7,6]", None)]
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
def tmin(expr, k=3):
    best = 1e9
    for _ in range(k):
        t = time.perf_counter(); pari(expr); best = min(best, time.perf_counter() - t)
    return best
res = []; c1_bad = 0; cn_tot = cn_disc = 0
for lab, a, base in CURVAS:
    pari(f"E=ellinit({a}); N=ellglobalred(E)[1]; ct=ellglobalred(E)[3]")
    pari(f"P={base}" if base else "P=ellsaturation(E, ellrank(E)[4], 1000)")
    pari("Pp=P[1]; for(k=2,#P,Pp=elladd(E,Pp,P[k]))")
    for m in [int(q) for q in pari(f"primes([7,{MMAX}])")]:
        if int(pari(f"N % {m}")) == 0: continue
        t = time.perf_counter(); pari(f"M=ct*ellcard(E,{m}); RR=vector(#P,i,ellmul(E,P[i],M))"); tp = time.perf_counter() - t
        if int(pari("my(z=0);for(i=1,#RR,if(RR[i]==[0],z=1));z")): continue
        tA = tmin(f"HA=matA(E,RR,{m})"); tB = tmin(f"HB=matB(E,RR,Pp,{m})")
        ok = int(pari("my(r=#RR,o=1);for(i=1,r,for(j=1,r,if(!eqpm(HA[i,j],HB[i,j]),o=0)));o"))
        c1_bad += (ok == 0)
        for i in range(1, int(pari("#RR")) + 1):
            cn_tot += 1; cn_disc += int(pari(f"!eqpm(A_wrong(E,RR[{i}],{m}),HB[{i},{i}])"))
        dig = int(pari("#Str(denominator(RR[1][1]))"))
        res.append(dict(curva=lab, m=m, t_pts=tp, t_A=tA, t_B=tB, rho=tB / tA, C1=ok, digitos=dig))
    print(lab, "ok", flush=True)
rhos = [x["rho"] for x in res]
med = statistics.median(rhos)
tot = sum(x["t_pts"] + x["t_B"] for x in res) / sum(x["t_pts"] + x["t_A"] for x in res)
porcurva = {c: statistics.median([x["rho"] for x in res if x["curva"] == c]) for c, _, _ in CURVAS}
resumo = dict(niveis=len(res), C1_desacordos=c1_bad, CN_discrimina=f"{cn_disc}/{cn_tot}",
              CN_passa=cn_disc >= 0.8 * cn_tot, rho_mediana=med, rho_min=min(rhos), rho_max=max(rhos),
              rho_mediana_por_curva=porcurva, B2_razao_total=tot,
              B1="CONFIRMADO" if (c1_bad == 0 and med >= 2) else "NAO CONFIRMADO")
if os.environ.get("PILOTO"):  # piloto cego: só mecânica e tempo de parede, sem nenhuma razão de tempo
    print("PILOTO niveis", len(res), "C1_desacordos", c1_bad, "CN_passa", resumo["CN_passa"])
else:
    json.dump(dict(resumo=resumo, niveis=res), open(f"resultado_B_m{MMAX}.json", "w"), indent=1)
    print(json.dumps(resumo, indent=1))
