import cypari2, math, sys, time, json
pari=cypari2.Pari(); pari.allocatemem(4*10**9); pari.default('realprecision',30)
rows=[l.split() for l in open('data/allbsd.00000-09999')]
sel=[r for r in rows if r[4]=='2' and r[2]=='1']
NMAX=int(sys.argv[1]); sel=[r for r in sel if int(r[0])<=NMAX]
out=[]; t0=time.time()
for r in sel:
    N,iso,a=r[0],r[1],r[3]; cre_sha=float(r[10])
    E=pari.ellinit(pari(a)); rec={'label':N+iso,'N':int(N),'sha_arch':cre_sha,'padic':{}}
    for p in [3,5,7,11,13]:
        if int(N)%p==0: rec['padic'][p]='bad'; continue
        ap=int(pari.ellap(E,p))
        typ='ord' if ap%p!=0 else 'ss'
        try:
            res=pari.ellpadicbsd(E,p,12)
            rp=int(res[0]); shp=res[1]
            rec['padic'][p]={'type':typ,'ap':ap,'r_p':rp,'sha_p':str(shp)}
        except Exception as e:
            rec['padic'][p]={'type':typ,'ap':ap,'err':str(e)[:80]}
    out.append(rec)
    flag=[p for p,v in rec['padic'].items() if isinstance(v,dict) and v.get('r_p',2)!=2]
    print(rec['label'],{p:(v if v=='bad' else (v.get('r_p'),v.get('sha_p',v.get('err')))) for p,v in rec['padic'].items()},'<<< ANOMALY' if flag else '',flush=True)
json.dump(out,open(f'padic_{NMAX}.json','w'))
print('done',len(out),'curves in',round(time.time()-t0),'s')
