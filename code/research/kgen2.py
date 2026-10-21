import cypari2, math, sys, pickle, time, json
pari=cypari2.Pari(); pari.allocatemem(3*10**9)
def log(*a): print(*a, flush=True)
def run(D, l, PREC, qmax, maxq):
    ainvs='[0,1,1,-2,0]'; N=389; tag=f'389a1_D{D}_l{l}'
    pari.default('realprecision',PREC); pari('read("heeg.gp")'); pari('read("twistg.gp")'); pari('read("divtest.gp")')
    t0=time.time()
    r=pari(f'my(E=ellinit({ainvs}),Dl={l*l*D},M=heegforms({N},Dl),keys=Mat(M)[,1]~,F,maxA,ymin,NT,an,pts); F=vector(#keys,i,mapget(M,keys[i])); maxA=vecmax(vector(#F,i,F[i][1])); ymin=sqrt(-Dl)/(2*maxA); NT=ceil({PREC}*log(10)/(2*Pi*ymin))+60; an=ellan(E,NT); pts=vector(#F,i,ellztopoint(E,heegz(E,an,F[i][1],F[i][2],F[i][3],Dl,{PREC}))); [pts,F,keys,maxA,NT]')
    pts,F,keys=r[0],r[1],[str(k) for k in r[2]]; h=len(keys); n=l+1; assert h==n,(h,n); log(f'[{tag}] {h} points, maxA={r[3]}, terms={r[4]}, {round(time.time()-t0)}s')
    T=pari('T'); f=pari(1)
    for p in pts: f*=(T-p[0])
    cv=pari.Vec(f); co=[pari.bestappr(pari.real(c),10**(PREC//3)) for c in cv]
    err=max(float(abs(c-rc)) for c,rc in zip(cv,co)); maxden=max(float(pari.log(pari.denominator(c))/2.302585) for c in co)
    log(f'  x-poly recognition err {err:.1e}, max log10 den {maxden:.1f}')
    if err>1e-30: raise Exception('recognition failed; increase PREC')
    f=sum(co[i]*T**(len(co)-1-i) for i in range(len(co))); f=f*pari.denominator(pari.content(f)); f=f/pari.content(f)
    t=time.time(); pr=pari.polredbest(f,1); g=pr[0]; A=pari.lift(pr[1]); B=pari.lift(pari.modreverse(pr[1])); log(f'  polredbest {round(time.time()-t)}s; g coeff max log10 {max(float(pari.log(abs(c)+1)/2.302585) for c in pari.Vec(g)):.1f}')
    th=[pari.subst(B,'T',p[0]) for p in pts]; idx={k:i for i,k in enumerate(keys)}
    def comp(a,b): return str(pari.qfbred(pari.qfbcomp(pari(a),pari(b))))
    Dl=l*l*D; ident=str(pari.qfbred(pari.Qfb(1,Dl%2,(Dl%2-Dl)//4)))
    def order(x):
        y=x;k=1
        while y!=ident: y=comp(y,x); k+=1
        return k
    gG=[x for x in keys if order(x)==n][0]; ginv=gG
    for _ in range(n-2): ginv=comp(ginv,gG)
    assert comp(ginv,gG)==ident
    R={}; Rd={}; sDr=pari(f'sqrt({-D})')
    for c in [1,2,3]:
        P_=pari(1)
        for k in keys: P_*=(T-(th[idx[k]]+c*(th[idx[comp(k,gG)]]+th[idx[comp(k,ginv)]])))
        cv=pari.Vec(P_); rec=[pari.round(pari.real(x)) for x in cv]; e=max(float(abs(x-rr)) for x,rr in zip(cv,rec)); R[c]=sum(rec[i]*T**(len(rec)-1-i) for i in range(len(rec)))
        if e>1e-20: raise Exception(f'resolvent c={c} recognition failed {e}')
    for c in [2,3]:
        P_=pari(1)
        for k in keys: P_*=(T-(th[idx[k]]+c*th[idx[comp(k,gG)]]))
        cv=pari.Vec(P_); 
        # coefficients in O_K: a + b*sqrt(D) with a,b rational (half-integers possible if D≡1 mod 4)
        ra=[pari.bestappr(pari.real(x),10**6) for x in cv]; rb=[pari.bestappr(pari.imag(x)/sDr,10**6) for x in cv]
        e=max(max(float(abs(pari.real(x)-a)),float(abs(pari.imag(x)/sDr-b))) for x,a,b in zip(cv,ra,rb))
        if e>1e-20: raise Exception(f'directed resolvent c={c} failed {e}')
        Rd[c]=(sum(ra[i]*T**(len(ra)-1-i) for i in range(len(ra))), sum(rb[i]*T**(len(rb)-1-i) for i in range(len(rb))))
    log('  resolvents ok')
    pari(f"g={g}; A={A}; nf=nfinit([g,10^5]);")
    zr=pari(f"my(xx=Mod(A,g), z2=(4*xx^3+4*xx^2-8*xx+1)/(4*({D})), pol=Pol([1,0,-lift(z2)],'x), r=nfroots(nf,pol)); [#r, if(#r,lift(nfbasistoalg(nf,r[1])),0)]"); mode='sqrtD'
    if int(zr[0])==0:
        zr=pari("my(xx=Mod(A,g), z2=(4*xx^3+4*xx^2-8*xx+1)/4, pol=Pol([1,0,-lift(z2)],'x), r=nfroots(nf,pol)); [#r, if(#r,lift(nfbasistoalg(nf,r[1])),0)]"); mode='plain'
    if int(zr[0])==0: raise Exception('eta not in L in either mode')
    log(f'  y-coordinate mode: {mode}'); ZE=zr[1]; sfac='sD' if mode=='sqrtD' else '1'
    pickle.dump({'g':str(g),'A':str(A),'ZE':str(ZE),'mode':mode,'R':{c:str(R[c]) for c in R},'Rd':{c:(str(Rd[c][0]),str(Rd[c][1])) for c in Rd}},open(f'kgen_{tag}.pkl','wb'))
    pari(f"ZE={ZE}; Ra2={Rd[2][0]}; Rb2={Rd[2][1]}; Ra3={Rd[3][0]}; Rb3={Rd[3][1]}; TW=twistpointD({D});")
    code=f'''
    kfull(qmin,qmax,maxq)={{my(cnt=0,P1=[-1,1],P2=[0,0],first=1,cands=0,hh={h},nn={n},zeroW=0);
     forprime(q=qmin,qmax, if(cnt>=maxq,break);
      if(kronecker({D},q)!=1 || q%3!=1 || q=={l}, next);
      my(rts=polrootsmod(g,q)); if(#rts!=hh, next);
      my(Eq=ellinit({ainvs},q), card=ellcard(Eq)); if(card%3!=0, next);
      my(sD=sqrt(Mod({D},q)), S2=Set(polrootsmod(Ra2+sD*Rb2,q)), S3=Set(polrootsmod(Ra3+sD*Rb3,q)), succ=vector(hh), Y=vector(hh), ok=1);
      for(j=1,hh, my(r=rts[j]); Y[j]=[subst(A,'T,r), -1/2+subst(ZE,'T,r)*{sfac}]; if(!ellisoncurve(Eq,Y[j]),ok=0));
      if(!ok, print("q=",q," not on curve"); next);
      for(j=1,hh, my(r=rts[j], found=List());
        for(a=1,hh, if(a==j,next); if(setsearch(S2,r+2*rts[a]) && setsearch(S3,r+3*rts[a]), listput(found,a)));
        if(#found!=1, ok=0; break); succ[j]=found[1]);
      if(!ok, print("q=",q," successor ambiguous"); next);
      my(cur=1, cyc=[1]); for(i=1,nn-1, cur=succ[cur]; cyc=concat(cyc,[cur])); if(succ[cur]!=1 || #Set(cyc)!=nn, print("q=",q," cycle problem"); next);
      my(W=[0],Tr=[0]); for(i=1,nn-1, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i))); for(i=1,nn, Tr=elladd(Eq,Tr,Y[cyc[i]]));
      my(trok=in3E(Eq,Tr,q));
      my(QE=twistredD(TW,{D},q,sD), gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%3==0, listput(G3,ellmul(Eq,gens[j],o/3))));
      my(ex=(q-1)/3, z=Mod(2,q)); while(z^ex==1, z=random(Mod(1,q))); z=z^ex;
      my(dl=(t)->if(t==1,0,if(t==z,1,2)), tp=(Z)->vector(#G3,j,dl(elltatepairing(Eq,G3[j],Z,3)^ex)));
      my(w=tp(W),p1=tp(P1),p2=tp(P2),qq=tp(QE), S=Set());
      for(a=0,2,for(b=0,2,for(c=0,2, my(v=vector(#G3,j,(a*p1[j]+b*p2[j]+c*qq[j])%3), v2=vector(#G3,j,(a*p1[j]+b*p2[j]-c*qq[j])%3)); if(v==w || v2==w, S=setunion(S,Set([[a,b,c]]))))));
      cnt++; cands=if(first,S,setintersect(cands,S)); first=0; if(w==vector(#G3,j,0), zeroW++);
      print("q=",q," dim=",#G3," Tr in 3E:",trok," w=",w," p1=",p1," p2=",p2," qtw=",qq,"  intersection:",cands));
     [cands, cnt, zeroW]}}'''
    pari(code); t=time.time(); res=pari(f'kfull(1000,{qmax},{maxq})'); log(f'[{tag}] FINAL candidates (a,b,c): {res[0]}  primes used {res[1]}, W trivial in {res[2]}; local {round(time.time()-t)}s total {round(time.time()-t0)}s')
    return str(res)
if __name__=='__main__':
    D=int(sys.argv[1]); l=int(sys.argv[2]); PREC=int(sys.argv[3]); qmax=int(sys.argv[4]); maxq=int(sys.argv[5])
    out=run(D,l,PREC,qmax,maxq); json.dump({'D':D,'l':l,'result':out},open(f'result_D{D}_l{l}.json','w'))
