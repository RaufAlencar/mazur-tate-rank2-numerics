import cypari2, math, sys, time
pari=cypari2.Pari(); pari.allocatemem(4*10**9); pari.default('realprecision',30)
rows=[l.split() for l in open('data/allbsd.00000-09999')]
sel=[r for r in rows if r[4]=='2' and r[2]=='1' and int(r[0])<=int(sys.argv[1])]
print(len(sel),'optimal rank-2 curves, N<=',sys.argv[1])
bad=[]; shas={}
t0=time.time()
for r in sel:
    N,iso,num,a=r[0],r[1],r[2],r[3]
    E=pari.ellinit(pari(a))
    ran,Lr=pari.ellanalyticrank(E); ran=int(ran); Lr=Lr/math.factorial(ran)
    rk=pari.ellrank(E); gens=rk[3]
    if len(gens)<2: bad.append((N+iso,'gens',rk)); continue
    sat=pari.ellsaturation(E,gens,10**3)
    reg=pari.matdet(pari.ellheightmatrix(E,sat))
    sha=float(Lr/(pari.ellbsd(E)*reg))
    cre=float(r[10])
    ok = ran==2 and abs(sha-round(sha))<1e-6 and abs(sha-cre)<1e-6
    shas[round(sha)]=shas.get(round(sha),0)+1
    if not ok: bad.append((N+iso,ran,sha,cre))
    print(f"{N+iso:7s} r_an={ran} Reg={float(reg):.8f} Sha={sha:.8f} cremona={cre:.0f} {'OK' if ok else 'MISMATCH'}",flush=True)
print('time',round(time.time()-t0),'s'); print('Sha distribution',shas); print('problems',bad)
