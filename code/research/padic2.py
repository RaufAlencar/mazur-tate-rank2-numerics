import cypari2, math, sys, time, json
pari=cypari2.Pari(); pari.allocatemem(4*10**9); pari.default('realprecision',30)
rows=[l.split() for l in open('data/allbsd.00000-09999')]
NMAX=int(sys.argv[1]); PR=[3,5,7,11,13]
sel=[r for r in rows if r[4]=='2' and r[2]=='1' and int(r[0])<=NMAX]
out=[]; t0=time.time(); anomalies=[]
for r in sel:
    N,iso,a=r[0],r[1],r[3]; E=pari.ellinit(pari(a))
    rk=pari.ellrank(E); gens=pari.ellsaturation(E,rk[3],10**3)
    if len(gens)!=2: print(N+iso,'no gens',rk); continue
    rec={'label':N+iso,'N':int(N),'sha_arch':float(r[10]),'p':{}}
    for p in PR:
        if int(N)%p==0: continue
        ap=int(pari.ellap(E,p))
        if ap%p==0: rec['p'][p]={'ss':True,'ap':ap}; continue
        try:
            rp,Lp=pari.ellpadicbsd(E,p,7)
            R=pari.ellpadicregulator(E,p,7,gens)
            vR=int(pari.valuation(R,p)); vL=int(pari.valuation(Lp,p))
            sha=Lp/R
            v=int(pari.valuation(sha,p)); 
            sha_int = pari.lift(pari.truncate(sha)) if v>=0 else None
            # check sha ≡ integer: compare mod p^k
            sha_mod = int(pari.lift(sha+pari(f'O({p}^6)'))) if v>=0 else None
            rec['p'][p]={'ap':ap,'r_p':int(rp),'vR':vR,'vL':vL,'sha_mod_p6':sha_mod,'anom':(ap==1)}
            if int(rp)!=2 or sha_mod!=1 or vR>=2: anomalies.append((N+iso,p,int(rp),vR,vL,sha_mod,ap))
        except Exception as e:
            rec['p'][p]={'err':str(e)[:60]}
    out.append(rec)
    s=' '.join(f"p{p}:{'ss' if v.get('ss') else ('r'+str(v.get('r_p'))+' vR'+str(v.get('vR'))+' Sha='+str(v.get('sha_mod_p6'))+(' ap=1!' if v.get('anom') else '')) if 'err' not in v else 'ERR'}" for p,v in rec['p'].items())
    print(N+iso, s, flush=True)
json.dump(out,open(f'padic2_{NMAX}.json','w'))
print('done',len(out),'curves',round(time.time()-t0),'s'); print('ANOMALIES (r_p!=2 or Sha!=1 or vR>=2):'); [print(x) for x in anomalies]
