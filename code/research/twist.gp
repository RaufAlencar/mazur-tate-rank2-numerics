twistpoint()={my(E=ellinit([0,1,1,-2,0]), D=-20, Es, ch, Et, Edm, chi, Qt, Qs, A, B);
  Es=ellinit([0,0,0,-27*E.c4,-54*E.c6]); ch=ellisisom(E,Es);
  Et=ellinit(elltwist(E,D)); Qt=ellsaturation(Et,ellrank(Et)[4],1000)[1];
  A=Es.a4; B=Es.a6; Edm=ellinit([0,0,0,A*D^2,B*D^3]); chi=ellisisom(Et,Edm); Qs=ellchangepoint(Qt,chi);
  [ch, Qt, Qs, ellisoncurve(Edm,Qs), E.c4, E.c6]}
\\ reduce twist point into E(F_q) given s5 = sqrt(-5) mod q: on Es: (x/D, y/(D*sqrtD)), sqrtD=2*s5; then inverse change ch back to E
twistred(Qs,ch,q,s5)={my(D=-20, sD=2*s5, Ps=[Mod(Qs[1],q)/D, Mod(Qs[2],q)/(D*sD)], E=ellinit([0,1,1,-2,0],q), Es=ellinit([0,0,0,-27*E.c4,-54*E.c6],q), P);
  if(!ellisoncurve(Es,Ps), error("twist point not on Es mod q"));
  P=ellchangepointinv(Ps,ch); if(!ellisoncurve(E,P), error("not on E")); P}
