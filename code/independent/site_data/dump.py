import cypari2, json, pathlib
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
here = pathlib.Path(__file__).resolve().parent
pari('read("%s")' % (here.parent / "mt_indep.gp").as_posix())
pari('dump(ai,P,ctam,uE,mmax)=my(E=ellinit(ai),MS=msinit_curve(E),N=ellglobalred(E)[1],out=List(),v);forprime(m=5,mmax,if(N%m==0,next);v=mt_check(E,MS,P,ctam,m);if(type(v[3])=="t_STR",listput(out,[m,v[2],-1,-1,v[2]]);next);listput(out,[m,v[2],v[5],lift(Mod(uE*v[4],v[2])),ctam*(m+1-ellap(E,m))]));Vec(out)')
res = {}
for lab, ai, P, ctam, uE, mmax in [("389a1","[0,1,1,-2,0]","[[-1,1],[0,0]]",1,8,400),
                                    ("5077a1","[0,0,1,-7,6]","[[-2,3],[-1,3],[0,2]]",1,-16,250)]:
    v = pari(f"dump({ai},{P},{ctam},{uE},{mmax})")
    rows = [[int(x) for x in r] for r in v]
    res[lab] = rows
    ok = sum(1 for r in rows if r[2] >= 0 and r[2] == r[3]); bad = sum(1 for r in rows if r[2] >= 0 and r[2] != r[3])
    print(lab, "niveis", len(rows), "iguais", ok, "diferentes", bad, "fora", sum(1 for r in rows if r[2] < 0))
(here.parent.parent.parent / "results/independent/site_levels_389a1_5077a1.json").write_text(json.dumps(res))
print([r for r in res["389a1"] if r[0] == 23])
