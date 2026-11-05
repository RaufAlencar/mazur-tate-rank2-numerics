dd(T)={sqrtint(denominator(T[1]))}
\\ Portillo-Bobadilla eq. (18): g(P,Q,P',p) = d(P'+P) d(P'+Q) / ( d(P') d(P'+P+Q) ) mod p  (Q in formal group at p)
gPB(E,P,Q,Pp,m)={my(a=dd(elladd(E,Pp,P)),b=dd(elladd(E,Pp,Q)),c=dd(Pp),ee=dd(elladd(E,Pp,elladd(E,P,Q)))); Mod(a*b,m)/Mod(c*ee,m)}
\\ ours: <R,S> = gamma/A(R+S) with gamma = constant term of primitive line through R,S; diagonal [A(R)]^-1 ; recipe (t/d)^2
ours(E,R,S,m)={my(lam,c,den,L,RS,A3); if(R==S, lam=(3*R[1]^2+2*E.a2*R[1]+E.a4-E.a1*R[2])/(2*R[2]+E.a1*R[1]+E.a3), lam=(S[2]-R[2])/(S[1]-R[1])); c=R[2]-lam*R[1];
  den=lcm(denominator(lam),denominator(c)); L=[den,-den*lam,-den*c]; L=L/content(L); RS=elladd(E,R,S); A3=numerator(RS[1]); Mod(L[3],m)/Mod(A3,m)}
cmpPB(E,P,m,Pp)={my(N=ellcard(E,m),R=ellmul(E,P[1],N),S=ellmul(E,P[2],N),g1,g2,g12,o1,o2,o12);
  g1=gPB(E,R,R,Pp,m); g2=gPB(E,S,S,Pp,m); g12=gPB(E,R,S,Pp,m);
  o1=ours(E,R,R,m); o2=ours(E,S,S,m); o12=ours(E,R,S,m);
  printf("m=%d  PB diag: %s %s  off: %s | ours diag: %s %s off: %s | A(R)^-1: %s | PB*ours: %s %s %s | PB/ours: %s %s %s\n",m,g1,g2,g12,o1,o2,o12,Mod(numerator(R[1]),m)^-1, g1*o1,g2*o2,g12*o12, g1/o1,g2/o2,g12/o12)}
