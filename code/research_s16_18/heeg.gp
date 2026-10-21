heegforms(N,Dl)={my(h=qfbclassno(Dl),beta,M=Map(),k=1,A,C,key);
  beta=0; for(b=0,2*N-1, if((b^2-Dl)%(4*N)==0, beta=b; break));
  while(#M<h && k<20000, A=N*k;
    forstep(B=-A+1,A,1, if((B-beta)%(2*N)==0 && (B^2-Dl)%(4*A)==0, C=(B^2-Dl)/(4*A);
      if(gcd(gcd(A,B),C)==1, key=Str(qfbred(Qfb(A,B,C))); if(!mapisdefined(M,key) || mapget(M,key)[1]>A, mapput(M,key,[A,B,C])))));
    k++);
  if(#M!=h, error("classes ",#M," != ",h)); M}
heegz(E,an,A,B,C,Dl,digits)={my(tau=(-B+sqrt(Dl))/(2*A),y=imag(tau),nt=ceil(digits*log(10)/(2*Pi*y))+50,q=exp(2*Pi*I*tau),s=0,qn=1);
  if(nt>#an,error("need ",nt," terms"));
  for(n=1,nt, qn*=q; s+=an[n]*qn/n); s}
run(N,D,l,digits)={my(Dl=l^2*D,E=ellinit([0,1,1,-2,0]),M=heegforms(N,Dl),keys=Mat(M)[,1]~,F,maxA,ymin,NT,an,xs,P,co,ok=1,maxmag=0,dens=[],a,b,ra,rb,bound,errm=0);
  F=vector(#keys,i,mapget(M,keys[i])); maxA=vecmax(vector(#F,i,F[i][1])); ymin=sqrt(-Dl)/(2*maxA); NT=ceil(digits*log(10)/(2*Pi*ymin))+60;
  an=ellan(E,NT); print("l=",l," h=",#keys," maxA=",maxA," terms=",NT);
  xs=vector(#F,i,ellztopoint(E,heegz(E,an,F[i][1],F[i][2],F[i][3],Dl,digits))[1]);
  P=prod(i=1,#xs,'T-xs[i]); co=Vec(P);
  for(i=2,#co, maxmag=max(maxmag,log(abs(co[i])+1)/log(10)));
  bound=10^floor(digits/2-maxmag-5);
  for(i=2,#co, a=real(co[i]); b=imag(co[i])/sqrt(5); ra=bestappr(a,bound); rb=bestappr(b,bound);
    errm=max(errm,max(abs(a-ra),abs(b-rb))); dens=concat(dens,[max(denominator(ra),denominator(rb))]));
  print("max log10|coef|=",maxmag,"  max recognition error=",errm,"  denominators(log10): ",vector(min(8,#dens),i,round(log(dens[i])/log(10))));
  [co,xs,F,keys]}
run2(N,D,l,digits)={my(Dl=l^2*D,E=ellinit([0,1,1,-2,0]),M=heegforms(N,Dl),keys=Mat(M)[,1]~,F,maxA,ymin,NT,an,pts);
  F=vector(#keys,i,mapget(M,keys[i])); maxA=vecmax(vector(#F,i,F[i][1])); ymin=sqrt(-Dl)/(2*maxA); NT=ceil(digits*log(10)/(2*Pi*ymin))+60;
  an=ellan(E,NT);
  pts=vector(#F,i,ellztopoint(E,heegz(E,an,F[i][1],F[i][2],F[i][3],Dl,digits)));
  [pts,F,keys]}
