twistpointD(D)={my(E=ellinit([0,1,1,-2,0]), Es, ch, Et, Edm, chi, Qt, Qs, A, B);
  Es=ellinit([0,0,0,-27*E.c4,-54*E.c6]); ch=ellisisom(E,Es);
  Et=ellinit(elltwist(E,D)); Qt=ellsaturation(Et,ellrank(Et)[4],1000)[1];
  A=Es.a4; B=Es.a6; Edm=ellinit([0,0,0,A*D^2,B*D^3]); chi=ellisisom(Et,Edm); Qs=ellchangepoint(Qt,chi);
  if(!ellisoncurve(Edm,Qs),error("twist")); [ch, Qs]}
twistredD(TW,D,q,sD)={my(Ps=[Mod(TW[2][1],q)/D, Mod(TW[2][2],q)/(D*sD)], E=ellinit([0,1,1,-2,0],q), Es=ellinit([0,0,0,-27*E.c4,-54*E.c6],q), P);
  if(!ellisoncurve(Es,Ps), error("twist point not on Es mod q")); P=ellchangepointinv(Ps,TW[1]); if(!ellisoncurve(E,P), error("not on E")); P}
