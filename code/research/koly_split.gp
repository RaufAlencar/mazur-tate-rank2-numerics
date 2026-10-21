kolysplit(g,A,ZE,R1,R2,R3,qmin,qmax,maxq)={my(cnt=0,res=List(),P1=[-1,1],P2=[0,0]);
 forprime(q=qmin,qmax, if(cnt>=maxq,break);
  if(kronecker(-5,q)!=1 || q%3!=1, next);
  my(den=lcm(lcm(denominator(content(A)),denominator(content(ZE))),denominator(content(g)))); if(den%q==0,next);
  my(rts=polrootsmod(g,q)); if(#rts!=36, next);
  my(Eq=ellinit([0,1,1,-2,0],q), card=ellcard(Eq)); if(card%3!=0, next);
  my(s5=sqrt(Mod(-5,q)), S1=Set(polrootsmod(R1,q)), S2=Set(polrootsmod(R2,q)), S3=Set(polrootsmod(R3,q)), nb=vector(36), Y=vector(36), ok=1);
  \\ points
  for(j=1,36, my(r=rts[j], xx=subst(A,'T,r), yy=-1/2+subst(ZE,'T,r)*s5); Y[j]=[xx,yy]; if(!ellisoncurve(Eq,Y[j]), ok=0));
  if(!ok, print("q=",q," point not on curve"); next);
  \\ neighbours: pair {a,b} with r+c(a+b) in S_c for c=1,2,3
  for(j=1,36, my(r=rts[j], found=List());
    for(a=1,36, for(b=a+1,36, if(a==j||b==j,next); my(s=rts[a]+rts[b]);
      if(setsearch(S1,r+s) && setsearch(S2,r+2*s) && setsearch(S3,r+3*s), listput(found,[a,b])))); 
    if(#found!=1, ok=0; break); nb[j]=found[1]);
  if(!ok, print("q=",q," neighbour structure ambiguous"); next);
  \\ walk cycle from root 1
  my(cyc=[1], prev=0, cur=1, nxt);
  for(i=1,17, nxt=if(nb[cur][1]==prev, nb[cur][2], nb[cur][1]); if(i==1, nxt=nb[cur][1]); prev=cur; cur=nxt; cyc=concat(cyc,[cur]));
  \\ closure check
  my(nxt2=if(nb[cur][1]==prev,nb[cur][2],nb[cur][1])); if(nxt2!=1, print("q=",q," cycle not length 18"); next);
  \\ trace check: sum over cycle = -6*Y1 for some sign
  my(Tr=[0]); for(i=1,18, Tr=elladd(Eq,Tr,Y[cyc[i]]));
  my(ii=sqrt(Mod(-1,q)), x1=Mod(-13,q)/4, Y1a=[x1,-1/2+33/8*ii], Y1b=[x1,-1/2-33/8*ii], trok= (Tr==ellmul(Eq,Y1a,-6)) || (Tr==ellmul(Eq,Y1b,-6)));
  \\ Kolyvagin derivative W = sum_{i=1}^{17} i * Y[cyc[i+1]]
  my(W=[0]); for(i=1,17, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i)));
  \\ 3-torsion and Tate pairing
  my(gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%3==0, listput(G3,ellmul(Eq,gens[j],o/3))));
  my(TW=twistpoint(), QE=twistred(TW[3],TW[1],q,s5), ex=(q-1)/3, tp=(Z)->vector(#G3,j,elltatepairing(Eq,G3[j],Z,3)^ex), vW=tp(W), v1=tp(P1), v2=tp(P2), vQ=tp(QE), one=vector(#G3,j,1), indep=1, sol=[]);
  if(v1==one||v2==one||v1==v2||vector(#G3,j,v1[j]*v2[j])==one, indep=0);
  for(a=0,2,for(b=0,2, if(vector(#G3,j,v1[j]^a*v2[j]^b)==vW, sol=concat(sol,[[a,b]]))));
  cnt++; my(aa=[]); for(a=0,2, if(vector(#G3,j,v1[j]^a)==vW, aa=concat(aa,[a]))); my(bb=[]); for(b=0,2, if(vector(#G3,j,v2[j]^b)==vW, bb=concat(bb,[b]))); my(cc=[]); for(c=0,2, if(vector(#G3,j,vQ[j]^c)==vW, cc=concat(cc,[c]))); print("q=",q," trace ok:",trok," Qtw triv:",vQ==one," W=c*Qtwist:",cc," dimE[3]=",#G3," indep:",indep," P1 triv:",v1==one," P2 triv:",v2==one," W triv:",vW==one," W=aP1:",aa," W=bP2:",bb," W=aP1+bP2:",sol);
  listput(res,[q,trok,#G3,indep,sol]));
 Vec(res)}
