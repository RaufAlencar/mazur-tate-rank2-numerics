in3Eff(El,Z)={my(n=ellcard(El), e=valuation(n,3)); if(e==0, return(1));
  my(m=n/3^e, Zp=ellmul(El,Z,m), gens=ellgenerators(El), T=List(), S=Set());
  for(j=1,#gens, my(o=ellorder(El,gens[j]), v=valuation(o,3)); if(v>0, listput(T,[ellmul(El,gens[j],o/3^v),3^v])));
  my(tot=prod(j=1,#T,T[j][2]/3));
  for(c=0,tot-1, my(cc=c,R=[0]); for(j=1,#T, my(k=cc%(T[j][2]/3)); cc=cc\(T[j][2]/3); R=elladd(El,R,ellmul(El,T[j][1],3*k))); S=setunion(S,Set([Str(R)])));
  setsearch(S,Str(Zp))>0}
