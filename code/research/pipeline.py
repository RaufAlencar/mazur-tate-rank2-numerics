import cypari2, math, sys, pickle, time
pari=cypari2.Pari(); pari.allocatemem(10*10**9)
def run_pipeline(ainvs, N, D, l, PREC, qmin, qmax, maxq, tag):
    pari.default('realprecision',PREC); pari('read("heeg.gp")')
    E=pari(f'ellinit({ainvs})')
    t=time.time(); r=pari(f'my(E=ellinit({ainvs}),Dl={l*l*D},M=heegforms({N},Dl),keys=Mat(M)[,1]~,F,maxA,ymin,NT,an,pts); F=vector(#keys,i,mapget(M,keys[i])); maxA=vecmax(vector(#F,i,F[i][1])); ymin=sqrt(-Dl)/(2*maxA); NT=ceil({PREC}*log(10)/(2*Pi*ymin))+60; an=ellan(E,NT); pts=vector(#F,i,ellztopoint(E,heegz(E,an,F[i][1],F[i][2],F[i][3],Dl,{PREC}))); [pts,F,keys]')
    pts,F,keys=r[0],r[1],[str(k) for k in r[2]]; h=len(keys); print(f'[{tag}] {h} Heegner points of conductor {l} computed in {round(time.time()-t)}s')
    T=pari('T'); f=pari(1)
    for p in pts: f*=(T-p[0])
    co=[pari.bestappr(pari.real(c),10**(PREC//3)) for c in pari.Vec(f)]
    err=max(float(abs(c-rc)) for c,rc in zip(pari.Vec(f),co)); print(f'  x-poly rational recognition error {err:.1e}')
    f=sum(co[i]*T**(len(co)-1-i) for i in range(len(co))); f=f*pari.denominator(pari.content(f)); f=f/pari.content(f)
    pr=pari.polredbest(f,1); g=pr[0]; A=pari.lift(pr[1]); B=pari.lift(pari.modreverse(pr[1]))
    th=[pari.subst(B,'T',p[0]) for p in pts]
    # group structure: G = kernel of genus map to Pic(O_K); if h(K)=1, G=all
    idx={k:i for i,k in enumerate(keys)}
    def comp(a,b): return str(pari.qfbred(pari.qfbcomp(pari(a),pari(b))))
    hK=int(pari.qfbclassno(D)); n=l+1
    ident=str(pari.qfbred(pari.Qfb(1,(l*l*D)%2,((l*l*D)%2-l*l*D)//4)))
    def order(x):
        y=x;k=1
        while y!=ident: y=comp(y,x); k+=1
        return k
    if hK==1: ker=keys
    else:
        # genus characters: use the one trivial on principal genus of O_K ... for D=-20 it's (-4/.), general: find character via norms
        raise Exception('hK>1 not handled generically here')
    gens=[x for x in ker if order(x)==n]; gG=gens[0]; ginv=comp(*([gG]*0+[gG,gG])) if False else None
    G=[ident]
    for i in range(n-1): G.append(comp(G[-1],gG))
    ginv=G[n-1]
    out={'g':str(g),'A':str(A)}
    for c in [1,2,3]:
        R=pari(1)
        for k in keys: R*=(T-(th[idx[k]]+c*(th[idx[comp(k,gG)]]+th[idx[comp(k,ginv)]])))
        cv=pari.Vec(R); rec=[pari.round(pari.real(x)) for x in cv]; e2=max(float(abs(x-rr)) for x,rr in zip(cv,rec))
        out[c]=str(sum(rec[i]*T**(len(rec)-1-i) for i in range(len(rec)))); print(f'  resolvent c={c}: integer recognition error {e2:.1e}')
    # zeta in L: (2y+a1 x+a3)/(2 sqrt(D)) ... general curve y^2+a1xy+a3y = ... : eta = 2y + a1 x + a3, eta^2 = 4x^3+b2x^2+2b4x+b6
    a1,a2,a3,a4,a6=[int(v) for v in pari(ainvs)]
    b2=a1*a1+4*a2; b4=2*a4+a1*a3; b6=a3*a3+4*a6
    pari(f"g={g}; A={A}; nf=nfinit([g,10^5]);")
    zr=pari(f"my(xx=Mod(A,g), z2=(4*xx^3+{b2}*xx^2+{2*b4}*xx+{b6})/(4*({D})), pol=Pol([1,0,-lift(z2)],'x), r=nfroots(nf,pol)); [#r, if(#r,lift(nfbasistoalg(nf,r[1])),0)]")
    mode='sqrtD'
    if int(zr[0])==0:
        zr=pari(f"my(xx=Mod(A,g), z2=(4*xx^3+{b2}*xx^2+{2*b4}*xx+{b6})/4, pol=Pol([1,0,-lift(z2)],'x), r=nfroots(nf,pol)); [#r, if(#r,lift(nfbasistoalg(nf,r[1])),0)]"); mode='plain'
    print('  eta/2 in L via', mode, ':', int(zr[0])>0); ZE=zr[1]; sDfac = 'sD' if mode=='sqrtD' else '1' 
    pickle.dump({**out,'ZE':str(ZE),'ainvs':ainvs,'D':D,'l':l},open(f'pipe_{tag}.pkl','wb'))
    # local test
    code=f'''
    loc(qmin,qmax,maxq)={{my(cnt=0,res=List(),hh={h},nn={n});
     forprime(q=qmin,qmax, if(cnt>=maxq,break);
      if(kronecker({D},q)!=1 || q%3!=1 || q=={l}, next);
      my(rts=polrootsmod(g,q)); if(#rts!=hh, next);
      my(Eq=ellinit({ainvs},q), card=ellcard(Eq)); if(card%3!=0, next);
      my(sD=sqrt(Mod({D},q)), S1=Set(polrootsmod(R1,q)), S2=Set(polrootsmod(R2,q)), S3=Set(polrootsmod(R3,q)), nb=vector(hh), Y=vector(hh), ok=1);
      for(j=1,hh, my(r=rts[j], xx=subst(A,'T,r), zz=subst(ZE,'T,r)); Y[j]=[xx, (zz*{sDfac}*2 - {a1}*xx - {a3})/2]; if(!ellisoncurve(Eq,Y[j]), ok=0));
      if(!ok, print("q=",q," not on curve"); next);
      for(j=1,hh, my(r=rts[j], found=List());
        for(a=1,hh, for(b=a+1,hh, if(a==j||b==j,next); my(s=rts[a]+rts[b]);
          if(setsearch(S1,r+s) && setsearch(S2,r+2*s) && setsearch(S3,r+3*s), listput(found,[a,b]))));
        if(#found!=1, ok=0; break); nb[j]=found[1]);
      if(!ok, next);
      my(cyc=[1],prev=0,cur=1,nxt); for(i=1,nn-1, nxt=if(nb[cur][1]==prev, nb[cur][2], nb[cur][1]); prev=cur; cur=nxt; cyc=concat(cyc,[cur]));
      if(#Set(cyc)!=nn, print("q=",q," bad cycle"); next);
      my(Tr=[0]); for(i=1,nn, Tr=elladd(Eq,Tr,Y[cyc[i]]));
      my(W=[0]); for(i=1,nn-1, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i)));
      cnt++; print("q=",q," group=",ellgroup(Eq)," Tr in 3E:",in3E(Eq,Tr,q),"  W in 3E: ",in3E(Eq,W,q)); listput(res,[q,in3E(Eq,W,q)]));
     Vec(res)}}'''
    pari(f"R1={out[1]}; R2={out[2]}; R3={out[3]}; ZE={ZE};"); pari('read("divtest.gp")'); pari(code)
    return pari(f'loc({qmin},{qmax},{maxq})')
if __name__=='__main__':
    print(run_pipeline('[0,1,1,0,0]',43,-8,29,200,1000,400000,15,'43a1_D-8_l29'))
