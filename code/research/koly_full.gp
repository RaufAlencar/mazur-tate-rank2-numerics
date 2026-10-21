\\ directed resolvents: Rd_c(T) = sum (a_i + b_i*sqrt(-5)) T^i ; mod q with s5
kolyfull(g,A,ZE,Ra2,Rb2,Ra3,Rb3,qmin,qmax,maxq)={my(cnt=0,P1=[-1,1],P2=[0,0],TW=twistpoint(),first=1,cands=0);
 forprime(q=qmin,qmax, if(cnt>=maxq,break);
  if(kronecker(-5,q)!=1 || q%3!=1, next);
  my(rts=polrootsmod(g,q)); if(#rts!=36, next);
  my(Eq=ellinit([0,1,1,-2,0],q), card=ellcard(Eq)); if(card%3!=0, next);
  my(s5=sqrt(Mod(-5,q)), S2=Set(polrootsmod(Ra2+s5*Rb2,q)), S3=Set(polrootsmod(Ra3+s5*Rb3,q)), succ=vector(36), Y=vector(36), ok=1);
  for(j=1,36, my(r=rts[j]); Y[j]=[subst(A,'T,r), -1/2+subst(ZE,'T,r)*s5]);
  for(j=1,36, my(r=rts[j], found=List());
    for(a=1,36, if(a==j,next); if(setsearch(S2,r+2*rts[a]) && setsearch(S3,r+3*rts[a]), listput(found,a)));
    if(#found!=1, ok=0; break); succ[j]=found[1]);
  if(!ok, print("q=",q," successor ambiguous"); next);
  \\ two oriented cycles
  my(used=vector(36), Wtot=[0], ncyc=0, Trs=[]);
  for(st=1,36, if(used[st],next); my(cur=st, cyc=[st]); used[st]=1; for(i=1,17, cur=succ[cur]; used[cur]=1; cyc=concat(cyc,[cur])); if(succ[cur]!=st, ok=0); ncyc++;
    my(W=[0],Tr=[0]); for(i=1,17, W=elladd(Eq,W,ellmul(Eq,Y[cyc[i+1]],i))); for(i=1,18, Tr=elladd(Eq,Tr,Y[cyc[i]])); Wtot=elladd(Eq,Wtot,W); Trs=concat(Trs,[in3E(Eq,Tr,q)]));
  if(!ok || ncyc!=2, print("q=",q," cycle problem"); next);
  my(QE=twistred(TW[3],TW[1],q,s5), gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%3==0, listput(G3,ellmul(Eq,gens[j],o/3))));
  my(ex=(q-1)/3, z=Mod(2,q)); while(z^ex==1, z=random(Mod(1,q))); z=z^ex;
  my(dl=(t)->if(t==1,0,if(t==z,1,2)), tp=(Z)->vector(#G3,j,dl(elltatepairing(Eq,G3[j],Z,3)^ex)));
  my(w=tp(Wtot),p1=tp(P1),p2=tp(P2),qq=tp(QE), S=Set());
  for(a=0,2,for(b=0,2,for(c=0,2, my(v=vector(#G3,j,(a*p1[j]+b*p2[j]+c*qq[j])%3), v2=vector(#G3,j,(a*p1[j]+b*p2[j]-c*qq[j])%3)); if(v==w || v2==w, S=setunion(S,Set([[a,b,c]]))))));
  cnt++; cands=if(first,S,setintersect(cands,S)); first=0;
  print("q=",q," dim=",#G3," traces in 3E:",Trs," w=",w," p1=",p1," p2=",p2," qtw=",qq,"  candidates (a,b,c):",S,"  intersection:",cands));
 cands}
