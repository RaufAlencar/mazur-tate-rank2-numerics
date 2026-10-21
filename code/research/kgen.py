import cypari2, math, sys, pickle, time, json
pari=cypari2.Pari(); pari.allocatemem(14*10**9)
def log(*a):
    print(*a, flush=True)
def run(l, PREC, qmax, maxq, tag):
    ainvs='[0,1,1,-2,0]'; N=389; D=-20
    pari.default('realprecision',PREC); pari('read("heeg.gp")'); pari('read("twist.gp")'); pari('read("divtest.gp")')
    t0=time.time()
    r=pari(f'my(E=ellinit({ainvs}),Dl={l*l*D},M=heegforms({N},Dl),keys=Mat(M)[,1]~,F,maxA,ymin,NT,an,pts); F=vector(#keys,i,mapget(M,keys[i])); maxA=vecmax(vector(#F,i,F[i][1])); ymin=sqrt(-Dl)/(2*maxA); NT=ceil({PREC}*log(10)/(2*Pi*ymin))+60; an=ellan(E,NT); pts=vector(#F,i,ellztopoint(E,heegz(E,an,F[i][1],F[i][2],F[i][3],Dl,{PREC}))); [pts,F,keys,maxA,NT]')
    pts,F,keys=r[0],r[1],[str(k) for k in r[2]]; h=len(keys); log(f'[{tag}] {h} points, maxA={r[3]}, terms={r[4]}, {round(time.time()-t0)}s')
    T=pari('T'); f=pari(1)
    for p in pts: f*=(T-p[0])
    cv=pari.Vec(f); co=[pari.bestappr(pari.real(c),10**(PREC//3)) for c in cv]
    err=max(float(abs(c-rc)) for c,rc in zip(cv,co)); maxden=max(float(pari.log(pari.denominator(c))/2.302585) for c in co)
    log(f'  x-poly recognition err {err:.1e}, max log10 den {maxden:.1f}, max imag {max(float(abs(pari.imag(c))) for c in cv):.1e}')
    if err>1e-30: raise Exception('recognition failed; increase PREC')
    f=sum(co[i]*T**(len(co)-1-i) for i in range(len(co))); f=f*pari.denominator(pari.content(f)); f=f/pari.content(f)
    t=time.time(); pr=pari.polredbest(f,1); g=pr[0]; A=pari.lift(pr[1]); B=pari.lift(pari.modreverse(pr[1])); log(f'  polredbest {round(time.time()-t)}s')
    th=[pari.subst(B,'T',p[0]) for p in pts]; idx={k:i for i,k in enumerate(keys)}
    def comp(a,b): return str(pari.qfbred(pari.qfbcomp(pari(a),pari(b))))
    Dl=l*l*D; ident=str(pari.qfbred(pari.Qfb(1,Dl%2,(Dl%2-Dl)//4)))
    def genus(k):
        fq=pari(k); a,b,c=int(fq[0]),int(fq[1]),int(fq[2])
        for v in (a,c,a+b+c):
            if v%2: return int(pari.kronecker(-4,v))
    ker=[k for k in keys if genus(k)==1]; n=l+1; assert len(ker)==n, (len(ker),n)
    def order(x):
        y=x;k=1
        while y!=ident: y=comp(y,x); k+=1
        return k
    gG=[x for x in ker if order(x)==n][0]
    # resolvents: undirected c=1,2,3 (Z) and directed c=2,3 (O_K)
    ginv=gG
    for _ in range(n-2): ginv=comp(ginv,gG)
    assert comp(ginv,gG)==ident
    R={}; Rd={}
    s5r=pari('sqrt(5)')
    for c in [1,2,3]:
        P_=pari(1)
        for k in keys: P_*=(T-(th[idx[k]]+c*(th[idx[comp(k,gG)]]+th[idx[comp(k,ginv)]])))
        cv=pari.Vec(P_); rec=[pari.round(pari.real(x)) for x in cv]; e=max(float(abs(x-rr)) for x,rr in zip(cv,rec)); R[c]=sum(rec[i]*T**(len(rec)-1-i) for i in range(len(rec))); log(f'  resolvent c={c} err {e:.1e}')
        if e>1e-20: raise Exception('resolvent recognition failed')
    for c in [2,3]:
        P_=pari(1)
        for k in keys: P_*=(T-(th[idx[k]]+c*th[idx[comp(k,gG)]]))
        cv=pari.Vec(P_); ra=[pari.round(pari.real(x)) for x in cv]; rb=[pari.round(pari.imag(x)/s5r) for x in cv]
        e=max(max(float(abs(pari.real(x)-a)),float(abs(pari.imag(x)/s5r-b))) for x,a,b in zip(cv,ra,rb)); log(f'  directed resolvent c={c} err {e:.1e}')
        if e>1e-20: raise Exception('directed resolvent recognition failed')
        Rd[c]=(sum(ra[i]*T**(len(ra)-1-i) for i in range(len(ra))), sum(rb[i]*T**(len(rb)-1-i) for i in range(len(rb))))
    pari(f"g={g}; A={A}; nf=nfinit([g,10^5]);")
    t=time.time(); zr=pari("my(xx=Mod(A,g), z2=-(4*xx^3+4*xx^2-8*xx+1)/20, pol=Pol([1,0,-lift(z2)],'x), r=nfroots(nf,pol)); [#r, if(#r,lift(nfbasistoalg(nf,r[1])),0)]")
    log(f'  zeta in L: {int(zr[0])>0} ({round(time.time()-t)}s)'); ZE=zr[1]
    pickle.dump({'g':str(g),'A':str(A),'ZE':str(ZE),'R':{c:str(R[c]) for c in R},'Rd':{c:(str(Rd[c][0]),str(Rd[c][1])) for c in Rd},'keys':keys},open(f'kgen_{tag}.pkl','wb'))
    pari(f"ZE={ZE}; R1={R[1]}; R2={R[2]}; R3={R[3]}; Ra2={Rd[2][0]}; Rb2={Rd[2][1]}; Ra3={Rd[3][0]}; Rb3={Rd[3][1]};")
    code=f'''
    kfull(qmin,qmax,maxq)={{my(cnt=0,P1=[-1,1],P2=[0,0],TW=twistpoint(),first=1,cands=0,hh={h},nn={n});
     forprime(q=qmin,qmax, if(cnt>=maxq,break);
      if(kronecker(-5,q)!=1 || q%3!=1 || q=={l}, next);
      my(rts=polrootsmod(g,q)); if(#rts!=hh, next);
      my(Eq=ellinit({ainvs},q), card=ellcard(Eq)); if(card%3!=0, next);
      my(s5=sqrt(Mod(-5,q)), S2=Set(polrootsmod(Ra2+s5*Rb2,q)), S3=Set(polrootsmod(Ra3+s5*Rb3,q)), succ=vector(hh), Y=vector(hh), ok=1);
      for(j=1,hh, my(r=rts[j]); Y[j]=[subst(A,'T,r), -1/2+subst(ZE,'T,r)*s5]; if(!ellisoncurve(Eq,Y[j]),ok=0));
      if(!ok, print("q=",q," not on curve"); next);
      my(rset=Map()); for(a=1,hh, mapput(rset,rts[a],a));
      for(j=1,hh, my(r=rts[j], found=List());
        for(a=1,hh, if(a==j,next); if(setsearch(S2,r+2*rts[a]) && setsearch(S3,r+3*rts[a]), listput(found,a)));
        if(#found!=1, ok=0; break); succ[j]=found[1]);
      if(!ok, print("q=",q," successor ambiguous"); next);
      my(used=vector(hh), Wtot=[0], ncyc=0, Trs=[]);
      for(st=1,hh, if(used[st],next); my(cur=st, cyc=[st]); used[st]=1; for(i=1,nn-1, cur=succ[cur]; used[cur]=1; cyc=concat(cyc,[cur])); if(succ[cur]!=st, ok=0); ncyc++;
        my(W=[0],Tr=[0]); for(i=1,nn-1, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i))); for(i=1,nn, Tr=elladd(Eq,Tr,Y[cyc[i]])); Wtot=elladd(Eq,Wtot,W); Trs=concat(Trs,[in3E(Eq,Tr,q)]));
      if(!ok || ncyc!=2, print("q=",q," cycle problem"); next);
      my(QE=twistred(TW[3],TW[1],q,s5), gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%3==0, listput(G3,ellmul(Eq,gens[j],o/3))));
      my(ex=(q-1)/3, z=Mod(2,q)); while(z^ex==1, z=random(Mod(1,q))); z=z^ex;
      my(dl=(t)->if(t==1,0,if(t==z,1,2)), tp=(Z)->vector(#G3,j,dl(elltatepairing(Eq,G3[j],Z,3)^ex)));
      my(w=tp(Wtot),p1=tp(P1),p2=tp(P2),qq=tp(QE), S=Set());
      for(a=0,2,for(b=0,2,for(c=0,2, my(v=vector(#G3,j,(a*p1[j]+b*p2[j]+c*qq[j])%3), v2=vector(#G3,j,(a*p1[j]+b*p2[j]-c*qq[j])%3)); if(v==w || v2==w, S=setunion(S,Set([[a,b,c]]))))));
      cnt++; cands=if(first,S,setintersect(cands,S)); first=0;
      print("q=",q," dim=",#G3," Tr in 3E:",Trs," w=",w," p1=",p1," p2=",p2," qtw=",qq,"  intersection:",cands));
     cands}}'''
    pari(code); t=time.time(); res=pari(f'kfull(1000,{qmax},{maxq})'); log(f'[{tag}] FINAL candidates (a,b,c): {res}   local phase {round(time.time()-t)}s, total {round(time.time()-t0)}s')
    return str(res)
if __name__=='__main__':
    l=int(sys.argv[1]); PREC=int(sys.argv[2]); qmax=int(sys.argv[3]); maxq=int(sys.argv[4])
    out=run(l,PREC,qmax,maxq,f'389a1_l{l}')
    json.dump({'l':l,'result':out},open(f'result_l{l}.json','w'))
