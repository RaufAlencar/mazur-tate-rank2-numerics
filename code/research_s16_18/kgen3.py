import cypari2, math, sys, pickle, time, json
pari=cypari2.Pari(); pari.allocatemem(3*10**9)
def log(*a): print(*a, flush=True)
def run(ainvs, D, l, p, PREC, qmax, maxq):
    E=pari(f'ellinit({ainvs})'); N=int(pari.ellglobalred(E)[0]); tag=f'{N}_{ainvs}_D{D}_l{l}_p{p}'.replace(' ','')
    pari.default('realprecision',PREC); pari('read("heeg.gp")'); pari('read("divp.gp")')
    gens=pari.ellsaturation(E,pari.ellrank(E)[3],1000); assert len(gens)==2; P1,P2=gens[0],gens[1]
    a1,a2,a3,a4,a6=[int(v) for v in pari(ainvs)]; b2=a1*a1+4*a2; b4=2*a4+a1*a3; b6=a3*a3+4*a6
    # twist generator mapped to short model
    pari(f"E0=ellinit({ainvs}); Es=ellinit([0,0,0,-27*E0.c4,-54*E0.c6]); CH=ellisisom(E0,Es); Et=ellinit(elltwist(E0,{D})); Qt=ellsaturation(Et,ellrank(Et)[4],1000)[1]; Edm=ellinit([0,0,0,Es.a4*({D})^2,Es.a6*({D})^3]); CHI=ellisisom(Et,Edm); QS=ellchangepoint(Qt,CHI); if(!ellisoncurve(Edm,QS),error(\"twist\"));")
    # prediction
    ker=[(i,j) for i in range(p) for j in range(p) if (i,j)!=(0,0) and int(pari(f"my(R=elladd(E0,ellmul(E0,{P1},{i}),ellmul(E0,{P2},{j}))); if(denominator(R[1])%{l}==0,1,inpE(ellinit({ainvs},{l}),R,{p}))"))]
    qd=int(pari(f"my(Em=ellminimalmodel(Et),Qm=ellchangepoint(Qt,ellisisom(Et,Em))); if(denominator(Qm[1])%{l}==0,1,inpE(ellinit(Em[1..5],{l}),Qm,{p}))"))
    pred='kappa=0' if (qd==1 or len(ker)==p*p-1) else f'line {ker}'
    log(f'[{tag}] PREDICTION: {pred}   (ker={ker}, twist gen {p}-divisible at l: {qd})')
    t0=time.time()
    r=pari(f'my(E=ellinit({ainvs}),Dl={l*l*D},M=heegforms({N},Dl),keys=Mat(M)[,1]~,F,maxA,ymin,NT,an,pts); F=vector(#keys,i,mapget(M,keys[i])); maxA=vecmax(vector(#F,i,F[i][1])); ymin=sqrt(-Dl)/(2*maxA); NT=ceil({PREC}*log(10)/(2*Pi*ymin))+60; an=ellan(E,NT); pts=vector(#F,i,ellztopoint(E,heegz(E,an,F[i][1],F[i][2],F[i][3],Dl,{PREC}))); [pts,F,keys,maxA,NT]')
    pts,F,keys=r[0],r[1],[str(k) for k in r[2]]; h=len(keys); n=l+1; assert h==n,(h,n); log(f'  {h} points, maxA={r[3]}, terms={r[4]}, {round(time.time()-t0)}s')
    T=pari('T'); f=pari(1)
    for pt in pts: f*=(T-pt[0])
    cv=pari.Vec(f); co=[pari.bestappr(pari.real(c),10**(PREC//3)) for c in cv]
    err=max(float(abs(c-rc)) for c,rc in zip(cv,co)); log(f'  x-poly recognition err {err:.1e}, max log10 den {max(float(pari.log(pari.denominator(c))/2.302585) for c in co):.1f}')
    if err>1e-30: raise Exception('recognition failed')
    f=sum(co[i]*T**(len(co)-1-i) for i in range(len(co))); f=f*pari.denominator(pari.content(f)); f=f/pari.content(f)
    t=time.time(); pr=pari.polredbest(f,1); g=pr[0]; A=pari.lift(pr[1]); B=pari.lift(pari.modreverse(pr[1])); log(f'  polredbest {round(time.time()-t)}s')
    th=[pari.subst(B,'T',pt[0]) for pt in pts]; idx={k:i for i,k in enumerate(keys)}
    def comp(a,b): return str(pari.qfbred(pari.qfbcomp(pari(a),pari(b))))
    Dl=l*l*D; ident=str(pari.qfbred(pari.Qfb(1,Dl%2,(Dl%2-Dl)//4)))
    def order(x):
        y=x;k=1
        while y!=ident: y=comp(y,x); k+=1
        return k
    gG=[x for x in keys if order(x)==n][0]; ginv=gG
    for _ in range(n-2): ginv=comp(ginv,gG)
    R={}; Rd={}; sDr=pari(f'sqrt({-D})')
    for c in [1,2,3]:
        P_=pari(1)
        for k in keys: P_*=(T-(th[idx[k]]+c*(th[idx[comp(k,gG)]]+th[idx[comp(k,ginv)]])))
        cv=pari.Vec(P_); rec=[pari.round(pari.real(x)) for x in cv]; e=max(float(abs(x-rr)) for x,rr in zip(cv,rec)); R[c]=sum(rec[i]*T**(len(rec)-1-i) for i in range(len(rec)))
        if e>1e-20: raise Exception('resolvent failed')
    for c in [2,3]:
        P_=pari(1)
        for k in keys: P_*=(T-(th[idx[k]]+c*th[idx[comp(k,gG)]]))
        cv=pari.Vec(P_); ra=[pari.bestappr(pari.real(x),10**6) for x in cv]; rb=[pari.bestappr(pari.imag(x)/sDr,10**6) for x in cv]
        e=max(max(float(abs(pari.real(x)-a)),float(abs(pari.imag(x)/sDr-b))) for x,a,b in zip(cv,ra,rb))
        if e>1e-20: raise Exception('directed resolvent failed')
        Rd[c]=(sum(ra[i]*T**(len(ra)-1-i) for i in range(len(ra))), sum(rb[i]*T**(len(rb)-1-i) for i in range(len(rb))))
    pari(f"g={g}; A={A}; nf=nfinit([g,10^5]);")
    zr=pari(f"my(xx=Mod(A,g), z2=(4*xx^3+{b2}*xx^2+{2*b4}*xx+{b6})/(4*({D})), pol=Pol([1,0,-lift(z2)],'x), r=nfroots(nf,pol)); [#r, if(#r,lift(nfbasistoalg(nf,r[1])),0)]"); mode='sqrtD'
    if int(zr[0])==0:
        zr=pari(f"my(xx=Mod(A,g), z2=(4*xx^3+{b2}*xx^2+{2*b4}*xx+{b6})/4, pol=Pol([1,0,-lift(z2)],'x), r=nfroots(nf,pol)); [#r, if(#r,lift(nfbasistoalg(nf,r[1])),0)]"); mode='plain'
    if int(zr[0])==0: raise Exception('eta not in L')
    ZE=zr[1]; sfac='sD' if mode=='sqrtD' else '1'; log(f'  y mode {mode}')
    pari(f"ZE={ZE}; Ra2={Rd[2][0]}; Rb2={Rd[2][1]}; Ra3={Rd[3][0]}; Rb3={Rd[3][1]};")
    code=f'''
    kfull(qmin,qmax,maxq)={{my(cnt=0,P1={P1},P2={P2},first=1,cands=0,hh={h},nn={n},pp={p},zeroW=0);
     forprime(q=qmin,qmax, if(cnt>=maxq,break);
      if(kronecker({D},q)!=1 || q%pp!=1 || q=={l} || {N}%q==0, next);
      if(denominator(P1[1])%q==0 || denominator(P2[1])%q==0 || denominator(QS[1])%q==0, next);
      my(rts=polrootsmod(g,q)); if(#rts!=hh, next);
      my(Eq=ellinit({ainvs},q), card=ellcard(Eq)); if(card%pp!=0, next);
      my(sD=sqrt(Mod({D},q)), S2=Set(polrootsmod(Ra2+sD*Rb2,q)), S3=Set(polrootsmod(Ra3+sD*Rb3,q)), succ=vector(hh), Y=vector(hh), ok=1);
      for(j=1,hh, my(r=rts[j], xx=subst(A,'T,r), zz=subst(ZE,'T,r)); Y[j]=[xx, (2*zz*{sfac} - {a1}*xx - {a3})/2]; if(!ellisoncurve(Eq,Y[j]),ok=0));
      if(!ok, print("q=",q," not on curve"); next);
      for(j=1,hh, my(r=rts[j], found=List());
        for(a=1,hh, if(a==j,next); if(setsearch(S2,r+2*rts[a]) && setsearch(S3,r+3*rts[a]), listput(found,a)));
        if(#found!=1, ok=0; break); succ[j]=found[1]);
      if(!ok, print("q=",q," successor ambiguous"); next);
      my(cur=1, cyc=[1]); for(i=1,nn-1, cur=succ[cur]; cyc=concat(cyc,[cur])); if(succ[cur]!=1 || #Set(cyc)!=nn, print("q=",q," cycle problem"); next);
      my(W=[0],Tr=[0]); for(i=1,nn-1, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i))); for(i=1,nn, Tr=elladd(Eq,Tr,Y[cyc[i]]));
      my(trok=inpE(Eq,Tr,pp));
      my(Esq=ellinit([0,0,0,-27*Eq.c4,-54*Eq.c6],q), Ps=[Mod(QS[1],q)/{D}, Mod(QS[2],q)/({D}*sD)], QE); if(!ellisoncurve(Esq,Ps), error("tw")); QE=ellchangepointinv(Ps,CH); if(!ellisoncurve(Eq,QE),error("tw2"));
      my(gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%pp==0, listput(G3,ellmul(Eq,gens[j],o/pp))));
      my(ex=(q-1)/pp, z=Mod(2,q)); while(z^ex==1, z=random(Mod(1,q))); z=z^ex; my(zt=vector(pp,i,z^(i-1)));
      my(dl=(t)->for(i=1,pp, if(t==zt[i], return(i-1))), tp=(Z)->vector(#G3,j,dl(elltatepairing(Eq,G3[j],Z,pp)^ex)));
      my(w=tp(W),p1=tp(P1),p2=tp(P2),qq=tp(QE), S=Set(), zero=vector(#G3,j,0));
      for(a=0,pp-1,for(b=0,pp-1,for(c=0,pp-1, if(a==0&&b==0&&c==0,next); my(v=vector(#G3,j,(a*p1[j]+b*p2[j]+c*qq[j])%pp), hit=0);
         if(w==zero, hit=(v==zero), for(k=1,pp-1, if(vector(#G3,j,(k*v[j])%pp)==w, hit=1)));
         if(hit, S=setunion(S,Set([[a,b,c]]))))));
      if(w==zero, S=setunion(S,Set([[0,0,0]])));
      cnt++; cands=if(first,S,setintersect(cands,S)); first=0; if(w==zero, zeroW++);
      print("q=",q," dim=",#G3," Tr ok:",trok," w=",w," p1=",p1," p2=",p2," qtw=",qq,"  #cands:",#cands));
     [cands, cnt, zeroW]}}'''
    open(f'kfull_{abs(hash(tag))%10**8}.gp','w').write(code); pari(f'read("kfull_{abs(hash(tag))%10**8}.gp")'); t=time.time(); res=pari(f'kfull(1000,{qmax},{maxq})')
    log(f'[{tag}] FINAL: prediction {pred} | computed candidates (a,b,c) up to scalar: {res[0]} | primes {res[1]}, W trivial in {res[2]} | {round(time.time()-t0)}s')
    return pred,str(res[0])
if __name__=='__main__':
    ainvs=sys.argv[1]; D=int(sys.argv[2]); l=int(sys.argv[3]); p=int(sys.argv[4]); PREC=int(sys.argv[5]); qmax=int(sys.argv[6]); maxq=int(sys.argv[7])
    pred,res=run(ainvs,D,l,p,PREC,qmax,maxq)
    with open('kgen3_results.jsonl','a') as fh: fh.write(json.dumps({'ainvs':ainvs,'D':D,'l':l,'p':p,'pred':pred,'res':res})+'\n')
