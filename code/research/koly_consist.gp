kolycons(g,A,ZE,R1,R2,R3,qmin,qmax,maxq)={my(cnt=0);
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
  \\ build the two cycles
  my(cycles=List(), used=vector(36));
  for(st=1,36, if(used[st],next); my(cyc=[st],prev=0,cur=st,nxt); used[st]=1;
    for(i=1,17, nxt=if(nb[cur][1]==prev, nb[cur][2], nb[cur][1]); prev=cur; cur=nxt; used[cur]=1; cyc=concat(cyc,[cur]));
    listput(cycles,cyc));
  if(#cycles!=2, print("q=",q," cycles:",#cycles); next);
  my(gens=ellgenerators(Eq), G3=List()); for(j=1,#gens, my(o=ellorder(Eq,gens[j])); if(o%3==0, listput(G3,ellmul(Eq,gens[j],o/3))));
  my(ex=(q-1)/3, tp=(Z)->vector(#G3,j,elltatepairing(Eq,G3[j],Z,3)^ex));
  my(classes=[]);
  for(c=1,2, my(cyc=cycles[c]);
    for(sh=0,2, for(dir=0,1,   \\ start shift and direction
      my(W=[0]); for(i=1,17, my(idx=if(dir==0, (sh+i)%18, (sh-i)%18)+1); W=elladd(Eq,W,ellmul(Eq,Y[cyc[idx]],i)));
      classes=concat(classes,[[c,sh,dir,tp(W)]]))));
  cnt++; print("q=",q," dimE[3]=",#G3); for(i=1,#classes, print("   cycle ",classes[i][1]," shift ",classes[i][2]," dir ",classes[i][3],"  pairing(W)=",classes[i][4])));}
