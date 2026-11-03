# Two-leg Kolyvagin class kappa_{l1 l2} in E(Q)/p for rank-3 E (tower method, no polredbest at level n).
# usage: python kgen2p.py ainvs D l1 l2 p PREC qmax maxq
import cypari2, sys, time, json, itertools, math
from hreps import best_reps, divisors_exact
pari=cypari2.Pari(); pari.allocatemem(3*10**9)
def log(*a): print(*a, flush=True)

def rec_K(v, D, bound):
    """complex v = a + b*sqrt(D), sqrt(D)=I*sqrt(-D) -> (a,b) rationals, error, log10 den"""
    a=pari.real(v); b=pari.imag(v)/pari(f'sqrt({-D})')
    ra=pari.bestappr(a,bound); rb=pari.bestappr(b,bound)
    err=max(float(abs(a-ra)),float(abs(b-rb)))
    den=max(int(pari.denominator(ra)),int(pari.denominator(rb)))
    return ra,rb,err,den

LASTD=[1]
def solve_K(rows, vals, D, bound, what, sol=None):
    if sol is None:
        M=pari.matrix(len(rows),len(rows[0]),[x for r in rows for x in r])
        sol=pari.matsolve(M,pari.Col(vals))
        if DUMP:  # numerical solution kept for offline recognition
            open(DUMP+'_'+''.join(ch for ch in what if ch.isalnum())+'.txt','w').write(str(sol))
    out=[]; emax=0; dmax=1
    for v in sol:
        ra,rb,e,d=rec_K(v,D,bound); out.append((ra,rb)); emax=max(emax,e); dmax=max(dmax,d)
    log(f'  {what}: {len(out)} coeffs in K, max rec err {emax:.1e}, max log10 den {len(str(dmax))-1}')
    if len(str(dmax))-1 > len(str(bound))-1-30:
        # fallback: common denominator by LLL (needs ~1.15x the digits of d instead of 2x)
        P=int(pari.default('realprecision'))
        RR=[pari.real(v) for v in sol]+[pari.imag(v)/pari(f'sqrt({-D})') for v in sol]
        res=pari('reccommon')(pari.Vec(RR),P,8); pre=1
        if int(res[0])==0 and LASTD[0]!=1:
            # retry with the previous common denominator as a known factor (y after x)
            pre=LASTD[0]; Pe=P-(len(str(pre))-1)
            log(f'  {what}: LLL failed ({res}); retry with known factor 10^{len(str(pre))-1}, effective precision {Pe}')
            res=pari('reccommon')(pari.Vec([pre*v for v in RR]),Pe,8)
        if int(res[0])==0: raise Exception(f'{what}: denominators reach the bestappr bound and LLL common denominator failed ({res})')
        d=res[0]*pre; nums=res[1]; m=len(sol); LASTD[0]=d
        out=[(nums[k]/d, nums[m+k]/d) for k in range(m)]
        log(f'  {what}: LLL common denominator, log10 d = {len(str(d))-1} (precision {P})')
        return out,0.0
    return out,emax

LEGPREC=500
DUMP=None
def leg(N, D, l, beta1, PREC, bound, ws, K):
    LP=min(PREC,LEGPREC); pari.default('realprecision',LP)
    keys,pts,_=level_points(N,D,l,beta1,1,LP,ws,K)
    pari.default('realprecision',PREC)
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]; h=len(keys)
    T=pari('T'); Dl=D*l*l
    if Dl%8==1 and Dl%3!=0:
        # small generator: Weber class polynomial; root <-> class matched through j = (f^24-16)^3/f^24
        g=pari(f'polclass({Dl},1,T)'); rts=pari(f'polroots({g})')
        jr=[pari(f'my(r={r}); (r^24-16)^3/r^24') for r in rts]; th=[]
        for k in keys:
            jk=pari(f'my(q=Vec({k})); ellj((-q[2]+sqrt({Dl}))/(2*q[1]))')
            d=[float(abs(j-jk)/(1+abs(jk))) for j in jr]; j0=min(range(len(d)),key=lambda i:d[i]); srt=sorted(d)
            if not (srt[0]<1e-30 and srt[1]>1e-10): raise Exception(f'Weber matching failed {srt[:2]}')
            th.append(rts[j0])
        log(f'  leg l={l}: h={h}, Weber class polynomial {g}')
    elif Dl%3!=0:
        # gamma2 class polynomial, reduced by polredbest (partial factorization); root <-> class via A(theta)^3 = j
        W=pari(f'polclass({Dl},5,T)'); pr=pari(f'polredbest([{W},10^5],1)'); g=pr[0]; Apol=pari.lift(pr[1])
        rts=pari(f'polroots({g})'); gr=[pari(f'subst({Apol},T,{r})^3') for r in rts]; th=[]
        for k in keys:
            jk=pari(f'my(q=Vec({k})); ellj((-q[2]+sqrt({Dl}))/(2*q[1]))')
            d=[float(abs(j-jk)/(1+abs(jk))) for j in gr]; j0=min(range(len(d)),key=lambda i:d[i]); srt=sorted(d)
            if not (srt[0]<1e-8 and srt[1]>1e6*srt[0]): raise Exception(f'gamma2 matching failed {srt[:2]}')
            th.append(rts[j0])
        log(f'  leg l={l}: h={h}, gamma2 class polynomial reduced to {str(g)[:80]}...')
    else:
        f=pari(1)
        for x in xs: f*=(T-x)
        cv=pari.Vec(f); co=[pari.bestappr(pari.real(c),bound) for c in cv]
        err=max(float(abs(c-rc)) for c,rc in zip(cv,co))
        lden=max(len(str(pari.denominator(c))) for c in co)-1
        log(f'  leg l={l}: h={h}, x-poly rec err {err:.1e}, max log10 den {lden}')
        if err>10**(-PREC//4): raise Exception('x-poly recognition failed')
        fq=sum(co[i]*T**(len(co)-1-i) for i in range(len(co))); fq=fq/pari.content(fq)
        pr=pari(f'polredbest([{fq},10^5],1)'); g=pr[0]; Apol=pari.lift(pr[1])
        rts=pari(f'polroots({g})'); th=[]
        for x in xs:
            dists=[float(abs(pari.subst(Apol,'T',r)-x)) for r in rts]; j=min(range(len(rts)),key=lambda i:dists[i])
            srt=sorted(dists)
            if not (srt[0]<1e-6 and srt[1]>1e6*srt[0]): raise Exception(f'theta matching failed {srt[:2]}')
            th.append(rts[j])
    if len(set(str(t) for t in th))!=h: raise Exception('theta matching not injective')
    pari(f'DLR=dlog({[k for k in keys]},{Dl})'.replace("'",'"'))
    gk=str(pari('DLR[1]')); s_of={k:int(pari(f'mapget(DLR[2],"{k}")')) for k in keys}
    idx={k:i for i,k in enumerate(keys)}
    succ=[idx[str(pari(f'qcomp("{k}","{gk}")'))] for k in keys]
    nfz=pari(f'nfinit([{g},10^6]).zk'); zk=[pari.lift(w) if False else w for w in nfz]
    V=[[pari.subst(w,'T',th[i]) for w in zk] for i in range(h)]
    # directed O_K-resolvents prod (T - (th_k + c th_{k g})), c=2,3  (as in kgen3)
    Rd={}; sDr=pari(f'sqrt({-D})')
    for c in [2,3]:
        P_=pari(1)
        for i in range(h): P_*=(T-(th[i]+c*th[succ[i]]))
        cv=pari.Vec(P_); ra=[pari.round(2*pari.real(x))/2 for x in cv]; rb=[pari.round(2*pari.imag(x)/sDr)/2 for x in cv]
        e=max(max(float(abs(pari.real(x)-a)),float(abs(pari.imag(x)/sDr-b))) for x,a,b in zip(cv,ra,rb))
        dd=max(max(int(pari.denominator(a)),int(pari.denominator(b))) for a,b in zip(ra,rb))
        log(f'  resolvent c={c} l={l}: err {e:.1e}, den {dd}')
        if e>10**(-PREC//4) or dd>4: raise Exception('directed resolvent failed')
        Rd[c]=(sum(ra[i]*T**(len(ra)-1-i) for i in range(len(ra))), sum(rb[i]*T**(len(rb)-1-i) for i in range(len(rb))))
    lb=pari(f'10^{LP//2}')
    X,_=solve_K(V,xs,D,lb,f'x in zk_{l}')
    Y,_=solve_K(V,ys,D,lb,f'y in zk_{l}')
    return dict(keys=keys,th=th,g=g,Rd=Rd,zk=zk,X=X,Y=Y,s_of=s_of,h=h,xs=xs,ys=ys)


def rep_coprime(key, n):
    a,b,c=[int(v) for v in pari(key)[:3]] if False else [int(v) for v in pari('Vec('+key+')')[:3]]
    Qv=lambda x,y: a*x*x+b*x*y+c*y*y
    for r in range(1,60):
        for x in range(-r,r+1):
            for y in (range(0,r+1)):
                if max(abs(x),abs(y))!=r and r>1: continue
                if math.gcd(x,y)!=1: continue
                v=Qv(x,y)
                if math.gcd(v,n)!=1: continue
                # find (r0,s0) with x*s0 - r0*y = 1
                g,u,w=egcd(x,-y)  # x*u + (-y)*w = g = +-1
                s0,r0=u*g,w*g
                assert x*s0-r0*y==1
                A=v; B=2*a*x*r0+b*(x*s0+r0*y)+2*c*y*s0; C=Qv(r0,s0)
                assert B*B-4*A*C==b*b-4*a*c
                return [A,B,C]
    raise Exception('no coprime rep')

def egcd(a,b):
    if b==0: return (a if a>=0 else -a, 1 if a>=0 else -1, 0)
    g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)

def level_points(N,D,m,beta1,cop,PREC,ws,K):
    br=best_reps(N,D,m,beta1,cop,pari,K=K,wsign=ws)
    keys=sorted(br); pts=[]
    for k in keys:
        A,B,C,Q,eps,im=br[k]
        pts.append(pari(f'ellztopoint(E,{eps}*hz2(an,{A},{B},{D*m*m},{PREC}))'))
    return keys,pts,br

def Kstr(c):  # (a,b) -> GP expression in sD
    return f'({c[0]}+({c[1]})*sD)'

KMAX=3000
def run(ainvs,D,l1,l2,p,PREC,qmax,maxq):
    n=l1*l2
    pari.default('realprecision',PREC); pari('read("two.gp")'); pari('read("divp.gp")'); pari('read("recc.gp")')
    pari(f'E=ellinit({ainvs})'); N=int(pari('ellglobalred(E)[1]'))
    G=pari('G=ellsaturation(E,ellrank(E)[4],1000)'); assert len(G)==3
    tag=f'{N}_D{D}_n{l1}x{l2}_p{p}'
    global DUMP; DUMP=f'runs/dump_{tag}_P{PREC}'
    beta1=int(pari(f'my(b=-1);for(x=0,2*{N}-1,if((x^2-({D}))%(4*{N})==0,b=x;break));b'))
    # prediction (recomputed here; registered separately beforehand)
    def ker(l): return set(v for v in itertools.product(range(p),repeat=3) if v!=(0,0,0) and int(pari(f'my(R=elladd(E,elladd(E,ellmul(E,G[1],{v[0]}),ellmul(E,G[2],{v[1]})),ellmul(E,G[3],{v[2]}))); if(R==[0]||denominator(R[1])%{l}==0,1,inpE(ellinit(E,{l}),R,{p}))')))
    I=sorted(ker(l1)&ker(l2))
    pred = 'kappa=0 (forced)' if len(I)!=p-1 else f'line {I} or 0'
    log(f'[{tag}] gens {G}  PREDICTION: {pred}')
    t0=time.time()
    ws={Q:int(pari(f'my(s=1); fordiv({Q},q,if(isprime(q),s*=-ellap(E,q))); s')) for Q in divisors_exact(N,pari)}
    log(f'  Atkin-Lehner signs {ws}')
    ymin=min(min(v[5] for v in best_reps(N,D,m,beta1,(n if m==n else 1),pari,K=KMAX,wsign=ws).values()) for m in (l1,l2,n))
    NT=int(PREC*2.302585/(6.283185*ymin))+100
    log(f'  terms {NT}'); pari(f'an=ellan(E,{NT})')
    # self-check of the Atkin-Lehner signs: same class computed from two orientations, low precision
    for Q in ws:
        if Q==1: continue
        from hreps import al_pairs
        pr=al_pairs(N,D,l1,beta1,Q,pari)[0]; (A,B,C),Fp=pr[1],pr[2]
        old=pari.default('realprecision'); pari.default('realprecision',40)
        z1=pari(f'hz(an,{A},{B},{D*l1*l1},30)'); z2=pari(f'{ws[Q]}*hz(an,{Fp[0]},{Fp[1]},{D*l1*l1},30)') if Fp[0]>0 else None
        d=pari(f'my(P1=ellztopoint(E,{z1}),P2=ellztopoint(E,{z2})); abs(P1[1]-P2[1])+abs(P1[2]-P2[2])')
        pari.default('realprecision',old)
        log(f'  AL check Q={Q}: |phi(tau)-w_Q phi(W_Q^-1 tau)| = {float(d):.1e}')
        if float(d)>1e-15: raise Exception('Atkin-Lehner sign check failed')
    bound=pari(f'10^{PREC//2}')
    L1=leg(N,D,l1,beta1,PREC,bound,ws,KMAX); L2=leg(N,D,l2,beta1,PREC,bound,ws,KMAX)
    log(f'  legs done {round(time.time()-t0)}s')
    import os
    dx=DUMP+'_xynintheta1itheta2j.txt'; dy=DUMP+'_yynintheta1itheta2j.txt'
    if os.environ.get('FROMDUMP') and os.path.exists(dx) and os.path.exists(dy):
        keys=sorted(best_reps(N,D,n,beta1,n,pari,K=KMAX,wsign=ws)); pts=[[0,0]]*len(keys)
        SOLX=pari(open(dx).read()); SOLY=pari(open(dy).read()); log(f'  level-n solutions read from dump {dx}')
    else:
        SOLX=SOLY=None
        t=time.time(); keys,pts,br=level_points(N,D,n,beta1,n,PREC,ws,KMAX); log(f'  {len(pts)} level-n points {round(time.time()-t)}s')
    rows=[];xs=[];ys=[];ST=[]
    for key,pnt in zip(keys,pts):
        F=rep_coprime(key,n); k1=str(pari(f'proj({F},{D},{l1},{l2})')); k2=str(pari(f'proj({F},{D},{l2},{l1})'))
        i1=L1['keys'].index(k1); i2=L2['keys'].index(k2)
        if SOLX is None: rows.append([pari.subst(w1,'T',L1['th'][i1])*pari.subst(w2,'T',L2['th'][i2]) for w1 in L1['zk'] for w2 in L2['zk']])
        xs.append(pnt[0]); ys.append(pnt[1]); ST.append((L1['s_of'][k1],L2['s_of'][k2]))
    assert len(set(ST))==len(pts)==L1['h']*L2['h'], 'class labels not a bijection'
    CX,ex=solve_K(rows,xs,D,bound,'x(y_n) in theta1^i theta2^j',SOLX)
    CY,ey=solve_K(rows,ys,D,bound,'y(y_n) in theta1^i theta2^j',SOLY)
    if max(ex,ey)>10**(-PREC//4): raise Exception('level-n recognition failed')
    h1,h2=L1['h'],L2['h']
    def poly(C): return '['+','.join(Kstr(c) for c in C)+']'
    pari(f'g1={L1["g"]}; g2={L2["g"]}; RA12={L1["Rd"][2][0]}; RB12={L1["Rd"][2][1]}; RA13={L1["Rd"][3][0]}; RB13={L1["Rd"][3][1]}; RA22={L2["Rd"][2][0]}; RB22={L2["Rd"][2][1]}; RA23={L2["Rd"][3][0]}; RB23={L2["Rd"][3][1]};')
    code=f'''
    ev(C,Z,r)=sum(k=1,#C,C[k]*subst(Z[k],'T,r));
    ev2(C,Z1,Z2,r1,r2)=my(s=0,a=vector(#Z1,i,subst(Z1[i],'T,r1)),b=vector(#Z2,j,subst(Z2[j],'T,r2)));for(i=1,#Z1,for(j=1,#Z2,s+=C[(i-1)*#Z2+j]*a[i]*b[j]));s;
    succv(rts,S2,S3)={{my(h=#rts,sc=vector(h));for(j=1,h,my(f=List());for(a=1,h,if(a!=j && setsearch(S2,rts[j]+2*rts[a]) && setsearch(S3,rts[j]+3*rts[a]),listput(f,a)));if(#f!=1,return(0));sc[j]=f[1]);sc}};
    cyc(rts,sc)={{my(v=vector(#rts),c=1);for(s=1,#rts,v[s]=rts[c];c=sc[c]);if(c!=1||#Set(v)!=#rts,return(0));v}};
    kq(qmin,qmax,maxq)={{my(cnt=0,first=1,cands=0,pp={p},zeroW=0,zeroW1=0,zeroW2=0,ST={[list(x) for x in ST]},al1=ellap(E,{l1}),al2=ellap(E,{l2}));
     forprime(q=qmin,qmax, if(cnt>=maxq,break);
      if(kronecker({D},q)!=1 || q%pp!=1 || {N}%q==0 || q=={n}, next);
      if(vecsum(vector(3,i,denominator(G[i][1])%q==0)), next);
      if(BADDEN%q==0, next);
      my(r1s=polrootsmod(g1,q)); if(#r1s!={h1}, next); my(r2s=polrootsmod(g2,q)); if(#r2s!={h2}, next);
      my(Eq=ellinit(E,q), card=ellcard(Eq)); if(card%pp!=0, next);
      my(sD=sqrt(Mod({D},q)), Z1={L1['zk']}, Z2={L2['zk']}, X1={poly(L1['X'])}, Y1={poly(L1['Y'])}, X2={poly(L2['X'])}, Y2={poly(L2['Y'])}, CX={poly(CX)}, CY={poly(CY)});
      my(sc1=succv(r1s,Set(polrootsmod(RA12+sD*RB12,q)),Set(polrootsmod(RA13+sD*RB13,q))), sc2=succv(r2s,Set(polrootsmod(RA22+sD*RB22,q)),Set(polrootsmod(RA23+sD*RB23,q))));
      if(sc1==0||sc2==0, print("q=",q," successor ambiguous"); next);
      my(v1=cyc(r1s,sc1),v2=cyc(r2s,sc2)); if(v1==0||v2==0, print("q=",q," Galois cycle problem"); next);
      my(P1=vector({h1},s,[ev(X1,Z1,v1[s]),ev(Y1,Z1,v1[s])]), P2=vector({h2},s,[ev(X2,Z2,v2[s]),ev(Y2,Z2,v2[s])]), ok=1);
      for(s=1,{h1},if(!ellisoncurve(Eq,P1[s]),ok=0)); for(s=1,{h2},if(!ellisoncurve(Eq,P2[s]),ok=0));
      if(!ok, print("q=",q," leg point not on curve"); next);
      my(Yn=vector(#ST), Tr1=vector({h1},s,[0]), Tr2=vector({h2},s,[0]));
      for(k=1,#ST, my(s=ST[k][1],t=ST[k][2]); Yn[k]=[ev2(CX,Z1,Z2,v1[s+1],v2[t+1]), ev2(CY,Z1,Z2,v1[s+1],v2[t+1])]; if(!ellisoncurve(Eq,Yn[k]),ok=0);
         if(ok, Tr1[s+1]=elladd(Eq,Tr1[s+1],Yn[k]); Tr2[t+1]=elladd(Eq,Tr2[t+1],Yn[k])));
      if(!ok, print("q=",q," level-n point not on curve"); next);
      my(trok=1); for(s=1,{h1}, if(Tr1[s]!=ellmul(Eq,P1[s],al2),trok=0)); for(t=1,{h2}, if(Tr2[t]!=ellmul(Eq,P2[t],al1),trok=0));
      my(T1=[0],T2=[0]); for(s=1,{h1},T1=elladd(Eq,T1,P1[s])); for(t=1,{h2},T2=elladd(Eq,T2,P2[t]));
      my(W=[0],W1=[0],W2=[0]);
      for(k=1,#ST, W=elladd(Eq,W,ellmul(Eq,Yn[k],ST[k][1]*ST[k][2])));
      for(s=1,{h1},W1=elladd(Eq,W1,ellmul(Eq,P1[s],s-1))); for(t=1,{h2},W2=elladd(Eq,W2,ellmul(Eq,P2[t],t-1)));
      my(gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%pp==0, listput(G3,ellmul(Eq,gens[j],o/pp))));
      my(ex=(q-1)/pp, z=Mod(2,q)); while(z^ex==1, z=random(Mod(1,q))); z=z^ex; my(zt=vector(pp,i,z^(i-1)));
      my(dl=(u)->for(i=1,pp, if(u==zt[i], return(i-1))), tp=(Z)->if(Z==[0],vector(#G3,j,0),vector(#G3,j,dl(elltatepairing(Eq,G3[j],Z,pp)^ex))));
      my(w=tp(W),w1=tp(W1),w2=tp(W2),pg=vector(3,i,tp(Mod(G[i],q))), S=Set(), zero=vector(#G3,j,0));
      for(a=0,pp-1,for(b=0,pp-1,for(c=0,pp-1, my(v=vector(#G3,j,(a*pg[1][j]+b*pg[2][j]+c*pg[3][j])%pp), hit=0);
         if(w==zero, hit=(v==zero), for(k=1,pp-1, if(vector(#G3,j,(k*v[j])%pp)==w, hit=1)));
         if(hit && (a||b||c||w==zero), S=setunion(S,Set([[a,b,c]]))))));
      cnt++; cands=if(first,S,setintersect(cands,S)); first=0; if(w==zero, zeroW++); if(w1==zero, zeroW1++); if(w2==zero, zeroW2++);
      print("q=",q," dim=",#G3," trace rel ok:",trok," Tr(y_l1)=0:",T1==[0]," Tr(y_l2)=0:",T2==[0]," w=",w," w_l1=",w1," w_l2=",w2," gens=",pg,"  #cands:",#cands));
     [cands, cnt, zeroW, zeroW1, zeroW2]}}'''
    dens=[c for C in (L1['X'],L1['Y'],L2['X'],L2['Y'],CX,CY) for ab in C for c in ab]
    pari('BADDEN=1'); 
    for c in dens: pari(f'BADDEN=lcm(BADDEN,denominator({c}))')
    for L in (L1,L2):
        for w in L['zk']: pari(f'BADDEN=lcm(BADDEN,denominator(content({w})))')
    log(f'  log10 lcm of denominators: {len(str(pari("BADDEN")))-1}')
    fn=f'kq_{tag}.gp'; open(fn,'w').write(code); pari(f'read("{fn}")')
    t=time.time(); res=pari(f'kq(1000,{qmax},{maxq})')
    log(f'[{tag}] FINAL: prediction {pred} | kappa_n candidates (a,b,c) up to scalar: {res[0]} | primes {res[1]}, W_n trivial in {res[2]}, W_l1 trivial in {res[3]}, W_l2 trivial in {res[4]} | {round(time.time()-t0)}s')
    return pred,str(res[0]),[int(x) for x in res[1:]]

if __name__=='__main__':
    ainvs=sys.argv[1]; D=int(sys.argv[2]); l1=int(sys.argv[3]); l2=int(sys.argv[4]); p=int(sys.argv[5]); PREC=int(sys.argv[6]); qmax=int(sys.argv[7]); maxq=int(sys.argv[8])
    if len(sys.argv)>9: KMAX=int(sys.argv[9])
    pred,res,cnts=run(ainvs,D,l1,l2,p,PREC,qmax,maxq)
    with open('kgen2p_results.jsonl','a') as fh: fh.write(json.dumps({'ainvs':ainvs,'D':D,'l1':l1,'l2':l2,'p':p,'PREC':PREC,'pred':pred,'res':res,'counts':cnts})+'\n')
