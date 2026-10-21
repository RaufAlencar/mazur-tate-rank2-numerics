\\ brute-force: is Z in 3E(F_q)? enumerate 3*E[3^inf]
in3E(Eq,Z,q)={my(n=ellcard(Eq), e=valuation(n,3), m=n/3^e, Zp=ellmul(Eq,Z,m), gens=ellgenerators(Eq), T=List(), S);
  \\ 3-primary generators
  for(j=1,#gens, my(o=ellorder(Eq,gens[j]), v=valuation(o,3)); if(v>0, listput(T,[ellmul(Eq,gens[j],o/3^v),3^v])));
  \\ enumerate 3*E[3^inf]
  S=Set(); my(idx=vector(#T));
  my(tot=prod(j=1,#T,T[j][2]/3));
  for(c=0,tot-1, my(cc=c, R=[0]); for(j=1,#T, my(k=cc%(T[j][2]/3)); cc\=(T[j][2]/3); R=elladd(Eq,R,ellmul(Eq,T[j][1],3*k))); S=setunion(S,Set([Str(R)])));
  setsearch(S,Str(Zp))>0}
