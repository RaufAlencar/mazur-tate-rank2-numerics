kolytest(f,l,D,qmin,qmax)={my(res=List(),P1=[-1,1],P2=[0,0],n=l+1,lc=pollead(f));
  forprime(q=qmin,qmax,
    if(q==17||q==389||q==3||q==5||lc%q==0, next);
    if(kronecker(D,q)!=1 || kronecker(-4,q)!=1, next);
    my(fq=factormod(f,q,1), degs=fq[,1]); if(#degs!=2 || degs[1]!=n || degs[2]!=n, next);
    my(h=factormod(f,q)[1,1], g=ffgen(h,'g), Eg=ellinit([0,1,1,-2,0],g), xr=g, yr, W=[0], card, gens3, G3=List());
    yr=(-1+sqrt(1+4*(xr^3+xr^2-2*xr)))/2;
    if(!ellisoncurve(Eg,[xr,yr]), error("not on curve"));
    for(i=1,l, W=elladd(Eg,W, ellmul(Eg,[xr^(q^i),yr^(q^i)],i)));
    card=ellcard(Eg); if(card%3!=0, print("q=",q," 3 does not divide #E"); next);
    gens3=ellgenerators(Eg);
    for(j=1,#gens3, my(o=ellorder(Eg,gens3[j])); if(o%3==0, listput(G3, ellmul(Eg,gens3[j],o/3))));
    if(#G3==0, next);
    my(ex=(q^n-1)/3, tp=(Z)->vector(#G3,j,elltatepairing(Eg,Z,G3[j],3)^ex));
    my(vW=tp(W), v1=tp(P1), v2=tp(P2), one=vector(#G3,j,1), ok=1, sol=[]);
    my(vWF=tp([W[1]^q,W[2]^q]));  \\ Frobenius of W: class must be Frob-invariant
    if(v1==one || v2==one || v1==v2 || vector(#G3,j,v1[j]*v2[j])==one, ok=0);
    for(a=0,2, for(b=0,2, if(vector(#G3,j,v1[j]^a*v2[j]^b)==vW, sol=concat(sol,[[a,b]]))));
    print("q=",q,"  dim E[3](F_q^18)=",#G3,"  P1,P2 indep mod 3E: ",ok,"  Frob-invariant: ",vW==vWF,"  W = a*P1+b*P2 mod 3E for (a,b) in ",sol);
    listput(res,[q,#G3,ok,vW==vWF,sol]));
  Vec(res)}
