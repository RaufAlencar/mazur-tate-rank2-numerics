dlog3(t,z)={if(t==1,0,if(t==z,1,2))}
kolysolve(g,A,ZE,R1,R2,R3,qmin,qmax,maxq)={my(cnt=0,cands=0,first=1,P1=[-1,1],P2=[0,0],TW=twistpoint());
 forprime(q=qmin,qmax, if(cnt>=maxq,break);
  if(kronecker(-5,q)!=1 || q%3!=1, next);
  my(rts=polrootsmod(g,q)); if(#rts!=36, next);
  my(Eq=ellinit([0,1,1,-2,0],q), card=ellcard(Eq)); if(card%3!=0, next);
  my(s5=sqrt(Mod(-5,q)), S1=Set(polrootsmod(R1,q)), S2=Set(polrootsmod(R2,q)), S3=Set(polrootsmod(R3,q)), nb=vector(36), Y=vector(36), ok=1);
  for(j=1,36, my(r=rts[j]); Y[j]=[subst(A,'T,r), -1/2+subst(ZE,'T,r)*s5]);
  for(j=1,36, my(r=rts[j], found=List());
    for(a=1,36, for(b=a+1,36, if(a==j||b==j,next); my(s=rts[a]+rts[b]);
      if(setsearch(S1,r+s) && setsearch(S2,r+2*s) && setsearch(S3,r+3*s), listput(found,[a,b]))));
    if(#found!=1, ok=0; break); nb[j]=found[1]);
  if(!ok, next);
  my(cyc=[1],prev=0,cur=1,nxt); for(i=1,17, nxt=if(nb[cur][1]==prev, nb[cur][2], nb[cur][1]); prev=cur; cur=nxt; cyc=concat(cyc,[cur]));
  my(W=[0]); for(i=1,17, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i)));
  my(QE=twistred(TW[3],TW[1],q,s5));
  my(gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%3==0, listput(G3,ellmul(Eq,gens[j],o/3))));
  my(ex=(q-1)/3, z=Mod(2,q)); while(z^ex==1, z=random(Mod(1,q))); z=z^ex;
  my(tp=(Z)->vector(#G3,j,dlog3(elltatepairing(Eq,G3[j],Z,3)^ex,z)));
  my(w=tp(W),p1=tp(P1),p2=tp(P2),qq=tp(QE), S=Set());
  for(a=0,2,for(b=0,2,for(c=0,2, my(v=vector(#G3,j,(a*p1[j]+b*p2[j]+c*qq[j])%3), v2=vector(#G3,j,(a*p1[j]+b*p2[j]-c*qq[j])%3), mw=vector(#G3,j,(-w[j])%3)); if(v==w || v==mw || v2==w || v2==mw, S=setunion(S,Set([[a,b,c]]))))));
  cnt++; cands=if(first, S, setintersect(cands,S)); first=0;
  print("q=",q," dim=",#G3," w=",w," p1=",p1," p2=",p2," q=",qq,"  local candidates:",#S,"  running intersection:",cands));
 cands}
