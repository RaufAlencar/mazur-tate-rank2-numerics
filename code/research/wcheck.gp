wcheck(g,A,ZE,R1,R2,R3,q)={
  my(rts=polrootsmod(g,q), Eq=ellinit([0,1,1,-2,0],q), s5=sqrt(Mod(-5,q)), S1=Set(polrootsmod(R1,q)), S2=Set(polrootsmod(R2,q)), S3=Set(polrootsmod(R3,q)), nb=vector(36), Y=vector(36), P1=[-1,1],P2=[0,0]);
  for(j=1,36, my(r=rts[j]); Y[j]=[subst(A,'T,r), -1/2+subst(ZE,'T,r)*s5]);
  for(j=1,36, my(r=rts[j], found=List());
    for(a=1,36, for(b=a+1,36, if(a==j||b==j,next); my(s=rts[a]+rts[b]);
      if(setsearch(S1,r+s) && setsearch(S2,r+2*s) && setsearch(S3,r+3*s), listput(found,[a,b]))));
    nb[j]=Vec(found));
  print("neighbour counts:",vector(36,j,#nb[j]));
  my(cyc=[1],prev=0,cur=1,nxt); for(i=1,17, nxt=if(nb[cur][1][1]==prev, nb[cur][1][2], nb[cur][1][1]); prev=cur; cur=nxt; cyc=concat(cyc,[cur]));
  my(W=[0]); for(i=1,17, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i)));
  my(D=ellsub(Eq,P1,P2));
  print("q=",q," group ",ellgroup(Eq)," W in 3E:",in3E(Eq,W,q)," W-(P1-P2) in 3E:",in3E(Eq,ellsub(Eq,W,D),q)," W+(P1-P2) in 3E:",in3E(Eq,elladd(Eq,W,D),q),
        " W-P1:",in3E(Eq,ellsub(Eq,W,P1),q)," W+P1:",in3E(Eq,elladd(Eq,W,P1),q)," W-P2:",in3E(Eq,ellsub(Eq,W,P2),q)," W+P2:",in3E(Eq,elladd(Eq,W,P2),q)," W-(P1+P2):",in3E(Eq,ellsub(Eq,W,elladd(Eq,P1,P2)),q)," W+(P1+P2):",in3E(Eq,elladd(Eq,W,elladd(Eq,P1,P2)),q));
  \\ also: order of W and size of 3-part
  print("   #E=",ellcard(Eq)," 3-adic val:",valuation(ellcard(Eq),3), " cycle:",cyc);}
