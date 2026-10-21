import cypari2, sys
pari=cypari2.Pari(); pari.allocatemem(2*10**9); pari.default('realprecision',40)
def check(ainvs,label=''):
    E=pari.ellinit(ainvs)
    N=pari.ellglobalred(E)[0]
    r,Lr=pari.ellanalyticrank(E)
    rk=pari.ellrank(E); gens=rk[3]
    rk_lo,rk_hi=int(rk[0]),int(rk[1])
    sat=pari.ellsaturation(E,gens,10**4) if len(gens)>0 else gens
    reg=pari.matdet(pari.ellheightmatrix(E,sat)) if len(sat)>0 else pari(1)
    bsd=pari.ellbsd(E)   # Omega*prod c_p / |tors|^2
    tors=pari.elltors(E)[0]
    tam=1
    for c in pari.ellglobalred(E)[3] if False else []: pass
    sha=Lr/(bsd*reg)
    print(f"{label:8s} N={N} r_an={int(r)} rk_alg=[{rk_lo},{rk_hi}] L^(r)/r!={float(Lr):.10f} Reg={float(reg):.10f} tors={tors} Sha_est={float(sha):.10f}")
    return float(sha)
if __name__=='__main__':
    for lab,a in [('389a1',[0,1,1,-2,0]),('433a1',[1,0,0,-6,6]),('446d1',[1,-1,0,-4,4]),('563a1',[0,1,1,-2,1]),('571b1',[0,-1,1,-929,-10595]),('643a1',[1,0,1,-2,1]),('655a1',[0,-1,1,-5,2]),('664a1',[0,0,0,-8,4]),('681c1',[1,1,0,-13,-16]),('707a1',[0,1,1,-6,1])]:
        try: check(a,lab)
        except Exception as e: print(lab,'ERR',e)
