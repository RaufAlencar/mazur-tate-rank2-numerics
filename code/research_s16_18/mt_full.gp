hval(Q,m,g,n)={my(x=Q[1],y=Q[2],d=sqrtint(denominator(x)),r=-x/y/d); znlog(Mod(numerator(r),m)/Mod(denominator(r),m),g)%n}
imgsize(E,m,P1,P2)={my(Eq=ellinit(E,m), Q1=[Mod(P1[1],m),Mod(P1[2],m)], Q2=[Mod(P2[1],m),Mod(P2[2],m)], o1, o2, inter=0);
  o1=ellorder(Eq,Q1); o2=ellorder(Eq,Q2);
  for(k=0,o2-1, my(T=ellmul(Eq,Q2,k)); if(T==[0], inter++, my(l=elllog(Eq,T,Q1,o1)); if(ellmul(Eq,Q1,l)==T, inter++)));
  o1*o2/inter}
mtfull(E,m,P1,P2)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], g=znprimroot(m), n=(m-1)/2, c=vector(n), den, a0, a1, a2, c2cls, N, R1, R2, h11, h22, hs, I, fourdet);
  for(a=1,m-1, my(t=znlog(Mod(a,m),g)%n); c[t+1]+=mseval(M,xp,[a/m,oo]));
  den=denominator(content(c)); c=c*den;
  a0=sum(t=0,n-1,c[t+1]); a1=sum(t=0,n-1,c[t+1]*t); a2=sum(t=0,n-1,c[t+1]*binomial(t,2));
  if(a0%n!=0, return([-1,"a0",a0%n])); if(a1%n!=0, return([-2,"a1",a1%n]));
  c2cls=(a2 - (a1/n)*(n*(n-1)/2))%n;
  N=ellcard(E,m); R1=ellmul(E,P1,N); R2=ellmul(E,P2,N);
  h11=hval(R1,m,g,n); h22=hval(R2,m,g,n); hs=hval(elladd(E,R1,R2),m,g,n);
  fourdet=(4*h11*h22-(hs-h11-h22)^2)%n;
  I=imgsize(E,m,P1,P2);
  [n, c2cls, N%n, I%n, fourdet, den]}
hvalk(Q,m,g,n,k)={my(x=Q[1],y=Q[2],d=sqrtint(denominator(x)),r=(-x/y)*d^k); znlog(Mod(numerator(r),m)/Mod(denominator(r),m),g)%n}
mtfullk(E,m,P1,P2,k)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], g=znprimroot(m), n=(m-1)/2, c=vector(n), den, a0, a1, a2, c2cls, N, R1, R2, h11, h22, hs, fourdet);
  for(a=1,m-1, my(t=znlog(Mod(a,m),g)%n); c[t+1]+=mseval(M,xp,[a/m,oo]));
  den=denominator(content(c)); c=c*den;
  a0=sum(t=0,n-1,c[t+1]); a1=sum(t=0,n-1,c[t+1]*t); a2=sum(t=0,n-1,c[t+1]*binomial(t,2));
  if(a0%n!=0 || a1%n!=0, return([-1]));
  c2cls=(a2 - (a1/n)*(n*(n-1)/2))%n;
  N=ellcard(E,m); R1=ellmul(E,P1,N); R2=ellmul(E,P2,N);
  my(S=elladd(E,R1,R2), out=vector(#k));
  for(i=1,#k, h11=hvalk(R1,m,g,n,k[i]); h22=hvalk(R2,m,g,n,k[i]); hs=hvalk(S,m,g,n,k[i]); out[i]=(4*h11*h22-(hs-h11-h22)^2)%n);
  [n, c2cls, N, out]}
hvalu(Q,m,g,n,a)={my(x=Q[1],y=Q[2],d=sqrtint(denominator(x)),t=-x/y,vt=valuation(t,m),vd=valuation(d,m),tu=t/m^vt,d0=d/m^vd,r=tu*d0^a); znlog(Mod(numerator(r),m)/Mod(denominator(r),m),g)%n}
mtfullu(E,m,P1,P2,as)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], g=znprimroot(m), n=(m-1)/2, c=vector(n), den, a0, a1, a2, c2cls, N, R1, R2, h11, h22, hs);
  for(a=1,m-1, my(t=znlog(Mod(a,m),g)%n); c[t+1]+=mseval(M,xp,[a/m,oo]));
  den=denominator(content(c)); c=c*den;
  a0=sum(t=0,n-1,c[t+1]); a1=sum(t=0,n-1,c[t+1]*t); a2=sum(t=0,n-1,c[t+1]*binomial(t,2));
  if(a0%n!=0 || a1%n!=0, return([-1]));
  c2cls=(a2 - (a1/n)*(n*(n-1)/2))%n;
  N=ellcard(E,m); R1=ellmul(E,P1,N); R2=ellmul(E,P2,N);
  my(S=elladd(E,R1,R2), out=vector(#as));
  for(i=1,#as, h11=hvalu(R1,m,g,n,as[i]); h22=hvalu(R2,m,g,n,as[i]); hs=hvalu(S,m,g,n,as[i]); out[i]=(4*h11*h22-(hs-h11-h22)^2)%n);
  [n, c2cls, N, out, den]}
mtfullc(E,m,P1,P2,cmul)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], g=znprimroot(m), n=(m-1)/2, c=vector(n), a0, a1, a2, c2cls, N, Mm, R1, R2, h11, h22, hs);
  for(a=1,m-1, my(t=znlog(Mod(a,m),g)%n); c[t+1]+=mseval(M,xp,[a/m,oo]));
  c=c*denominator(content(c));
  a0=sum(t=0,n-1,c[t+1]); a1=sum(t=0,n-1,c[t+1]*t); a2=sum(t=0,n-1,c[t+1]*binomial(t,2));
  if(a0%n!=0 || a1%n!=0, return([-1]));
  c2cls=(a2 - (a1/n)*(n*(n-1)/2))%n;
  N=ellcard(E,m); Mm=cmul*N; R1=ellmul(E,P1,Mm); R2=ellmul(E,P2,Mm);
  h11=hvalu(R1,m,g,n,-1); h22=hvalu(R2,m,g,n,-1); hs=hvalu(elladd(E,R1,R2),m,g,n,-1);
  [n, c2cls, Mm, (4*h11*h22-(hs-h11-h22)^2)%n]}
