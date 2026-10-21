\\ class of theta_m in I^r/I^{r+1} for Z[Z/n], n=(m-1)/2, with relations from (1+x)^n=1, assuming theta in I^r
thetaclass(E,m,r)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], g=znprimroot(m), n=(m-1)/2, c=vector(n), a, k, cls);
  for(b=1,m-1, my(t=znlog(Mod(b,m),g)%n); c[t+1]+=mseval(M,xp,[b/m,oo]));
  c=c*denominator(content(c));
  a=vector(r+1,j,sum(t=0,n-1,c[t+1]*binomial(t,j-1)));   \\ a[j] = coefficient of x^(j-1)
  if(a[1]%n!=0, return([-1,0]));
  \\ successive reduction: for j=1..r-1, a[j] = n*k_j, substitute n x^(j-1) ... careful: use relation n*x^(j) = -sum_{i>=2} C(n,i) x^(j-1+i)
  for(j=2,r, if(a[j]%n!=0, return([-j,a[j]%n])); k=a[j]/n; for(i=2,r+2-j, a[j-1+i]-=k*binomial(n,i)));
  [n, a[r+1]%n]}
pairdet(E,m,P,Mm,n,g)={my(r=#P, R=vector(r,i,ellmul(E,P[i],Mm)), H=matrix(r,r), hv=vector(r,i,hvalu(R[i],m,g,n,-1)));
  for(i=1,r, H[i,i]=2*hv[i]; for(j=i+1,r, my(s=hvalu(elladd(E,R[i],R[j]),m,g,n,-1)); H[i,j]=H[j,i]=s-hv[i]-hv[j]));
  matdet(H)%n}   \\ = 2^r det(h)
