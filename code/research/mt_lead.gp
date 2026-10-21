\\ theta_m in F_3[x]/(x^9) using primitive root g; returns [ord, c2 = coefficient of x^ord... we return full coeff vector]
thetacoefs(E,m,g,p,e)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], n=p^e, vals=vector(n), pol);
  for(a=1,m-1, my(t=znlog(Mod(a,m),g)%n, v=mseval(M,xp,[a/m,oo])); vals[t+1]+=v);
  if(denominator(content(vals))%p==0, return(0));
  pol=sum(t=0,n-1, vals[t+1]*(1+'x)^t)*Mod(1,p); vector(n,i,polcoeff(pol,i-1))}
\\ Mazur-Tate pairing matrix via t/d of N*P (reconstructed formula), cube classes w.r.t. g
cub(h,m,g)={my(z=Mod(h,m)^((m-1)/3), w=Mod(g,m)^((m-1)/3)); if(z==1,0,if(z==w,1,2))}
aval(E,Q,m,g)={my(x=Q[1],y=Q[2],d=sqrtint(denominator(x)),t=-x/y, r=t/d); cub(lift(Mod(numerator(r),m)/Mod(denominator(r),m)),m,g)}
mtreg(E,m,g,P1,P2)={my(N=ellcard(E,m),R1,R2,A11,A22,A12);
  if(N%3==0, return([-1]));
  R1=ellmul(E,P1,N); R2=ellmul(E,P2,N);
  A11=aval(E,R1,m,g); A22=aval(E,R2,m,g); A12=((aval(E,elladd(E,R1,R2),m,g)-A11-A22)*2)%3;
  [(N*(A11*A22-A12^2))%3, [A11,A12,A22], N%3]}
