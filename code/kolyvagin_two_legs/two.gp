\\ Heegner points of conductor m with a FIXED orientation: B == m*beta1 (mod 2N), beta1^2 == D (mod 4N).
\\ Forms [A,B,C], A = N*k with gcd(k,cop)==1; keeps the smallest A per class (largest Im tau).
hforms(N,D,m,beta1,cop)={my(Dm=D*m^2,h=qfbclassno(Dm),bm=(m*beta1)%(2*N),M=Map(),k=1,A,C,key);
  while(#M<h && k<200000, if(gcd(k,cop)==1, A=N*k;
    forstep(B=bm-2*N*ceil((A+bm)/(2*N)),A,2*N, if(B>-A && (B^2-Dm)%(4*A)==0, C=(B^2-Dm)/(4*A);
      if(gcd(gcd(A,B),C)==1, key=Str(qfbred(Qfb(A,B,C))); if(!mapisdefined(M,key) || mapget(M,key)[1]>A, mapput(M,key,[A,B,C]))))));
    k++);
  if(#M!=h, error("classes ",#M," != ",h)); M}
\\ modular parametrization value z(tau) = sum a_n/n q^n
hz(an,A,B,Dm,digits)={my(tau=(-B+sqrt(Dm))/(2*A),nt=ceil(digits*log(10)/(2*Pi*imag(tau)))+50,q=exp(2*Pi*I*tau),s=0.,qn=1);
  if(nt>#an,error("need ",nt," terms"));
  for(n=1,nt, qn*=q; if(an[n], s+=an[n]*qn/n)); s}
\\ list of [key, form, point] for conductor m
hpoints(E,an,N,D,m,beta1,cop,digits)={my(M=hforms(N,D,m,beta1,cop),keys=Mat(M)[,1]~);
  vector(#keys,i,my(F=mapget(M,keys[i])); [keys[i],F,ellztopoint(E,hz(an,F[1],F[2],D*m^2,digits))])}
\\ projection of a form of disc D*n^2 (gcd(A,n)=1) to disc D*l^2, n = l*l2
proj(F,D,l,l2)={my(A=F[1],u=lift(Mod(l2,2*A)^-1),B=(u*F[2])%(2*A),C);
  if(B>A,B-=2*A); C=(B^2-D*l^2)/(4*A); if(denominator(C)!=1,error("proj")); Str(qfbred(Qfb(A,B,C)))}
qcomp(a,b)=Str(qfbred(qfbcomp(eval(a),eval(b))));
\\ discrete log table: Map key -> s with key = g^s, g of order h (cyclic group assumed)
dlog(keys,Dl)={my(id=Str(qfbred(Qfb(1,Dl%2,(Dl%2-Dl)/4))),h=#keys,g,x,M=Map(),s);
  for(i=1,h, g=keys[i]; x=id; M=Map(); for(s=0,h-1, mapput(M,x,s); x=qcomp(x,g)); if(#M==h, return([g,M])));
  error("not cyclic")}
\\ same value as hz, baby-step/giant-step: inner blocks as C-level dot products
hz2(an,A,B,Dm,digits)={my(tau=(-B+sqrt(Dm))/(2*A),nt=ceil(digits*log(10)/(2*Pi*imag(tau)))+50,q=exp(2*Pi*I*tau),bs,J,Qp,qB,s=0);
  if(nt>#an,error("need ",nt," terms"));
  bs=ceil(sqrt(nt)); J=ceil(nt/bs); Qp=vector(bs); Qp[1]=q; for(i=2,bs,Qp[i]=Qp[i-1]*q); qB=Qp[bs];
  forstep(j=J-1,0,-1, my(o=j*bs, c=vector(bs,i,my(n=o+i); if(n<=nt,an[n]/n,0))); s=s*qB+c*Qp~);
  s}
