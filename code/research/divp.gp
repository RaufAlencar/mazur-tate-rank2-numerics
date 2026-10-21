inpE(El,Z,p)={my(n=ellcard(El), e=valuation(n,p)); if(e==0, return(1));
  my(m=n/p^e, Zp=ellmul(El,Z,m), gens=ellgenerators(El), T=List(), S=Set());
  for(j=1,#gens, my(o=ellorder(El,gens[j]), v=valuation(o,p)); if(v>0, listput(T,[ellmul(El,gens[j],o/p^v),p^v])));
  my(tot=prod(j=1,#T,T[j][2]/p));
  for(c=0,tot-1, my(cc=c,R=[0]); for(j=1,#T, my(k=cc%(T[j][2]/p)); cc=cc\(T[j][2]/p); R=elladd(El,R,ellmul(El,T[j][1],p*k))); S=setunion(S,Set([Str(R)])));
  setsearch(S,Str(Zp))>0}
