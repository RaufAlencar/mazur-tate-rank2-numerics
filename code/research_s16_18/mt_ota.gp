\\ Composite-level Mazur-Tate / Ota test (rank 2), S = m1*m2, pushed to Gbar = Z/n1 x Z/n2, n_i=(m_i-1)/2.
\\ theta_code = Den * sum_{a in (Z/S)^x} [a/S]^+ sigma_a  (Den fixed per curve: 1 if c_inf=2, 2 if c_inf=1)
\\ Returns [n1,n2,g, a20 mod n1, a11 mod g, a02 mod n2, M, D1 mod n1, D12 mod g, D2 mod n2, a_{m1}, a_{m2}, prime-level c2(m1), c2(m2)]
\\ with H^k_{ij} = polarised pairing matrix B at level m_k (B_ii = 2 h_ii), R_i = M P_i, M = Tam*N1*N2.

hv2(Q,m,g,n)={my(x=Q[1],y=Q[2],d=sqrtint(denominator(x)),t=-x/y,vt=valuation(t,m),vd=valuation(d,m),r=(t/m^vt)/(d/m^vd)); znlog(Mod(numerator(r),m)/Mod(denominator(r),m),g)%n}
Bmat(E,m,g,n,R)={my(r=#R,H=matrix(r,r),hv=vector(r,i,hv2(R[i],m,g,n)));
  for(i=1,r, H[i,i]=2*hv[i]; for(j=i+1,r, my(s=hv2(elladd(E,R[i],R[j]),m,g,n)); H[i,j]=H[j,i]=s-hv[i]-hv[j])); H}

\\ prime-level leading coefficient c2(m) (same reduction as thetaclass), with fixed Den
c2prime(Ms,m,Den)={my(M=Ms[1],xp=Ms[2],g=znprimroot(m),n=(m-1)/2,c=vector(n),a0,a1,a2);
  for(b=1,m-1, my(t=znlog(Mod(b,m),g)%n); c[t+1]+=mseval(M,xp,[b/m,oo]));
  c=c*Den; if(denominator(content(c))!=1, return([-9]));
  a0=sum(t=0,n-1,c[t+1]); a1=sum(t=0,n-1,c[t+1]*t); a2=sum(t=0,n-1,c[t+1]*binomial(t,2));
  if(a0!=0 || a1%n!=0, return([-1,a0,a1%n]));
  [n,(a2-(a1/n)*binomial(n,2))%n]}

otacomp(E,Ms,m1,m2,P,Den)={my(M=Ms[1],xp=Ms[2],S=m1*m2,g1=znprimroot(m1),g2=znprimroot(m2),n1=(m1-1)/2,n2=(m2-1)/2,gg=gcd(n1,n2),
   c=matrix(n1,n2),a00,a10,a01,a20,a11,a02,k1,k2,cT,N1,N2,Mm,R,B1,B2,D1,D2,D12,p1,p2);
  for(a=1,S-1, if(gcd(a,S)==1, my(t1=znlog(Mod(a,m1),g1)%n1, t2=znlog(Mod(a,m2),g2)%n2); c[t1+1,t2+1]+=mseval(M,xp,[a/S,oo])));
  c=c*Den; if(denominator(content(c))!=1, return([-9,content(c)]));
  a00=sum(i=1,n1,sum(j=1,n2,c[i,j]));
  a10=sum(i=1,n1,sum(j=1,n2,c[i,j]*(i-1)));  a01=sum(i=1,n1,sum(j=1,n2,c[i,j]*(j-1)));
  a20=sum(i=1,n1,sum(j=1,n2,c[i,j]*binomial(i-1,2))); a02=sum(i=1,n1,sum(j=1,n2,c[i,j]*binomial(j-1,2)));
  a11=sum(i=1,n1,sum(j=1,n2,c[i,j]*(i-1)*(j-1)));
  if(a00!=0, return([-1,"a00",a00]));
  if(a10%n1!=0 || a01%n2!=0, return([-2,"not in I^2",a10%n1,a01%n2]));
  k1=a10/n1; k2=a01/n2; a20-=k1*binomial(n1,2); a02-=k2*binomial(n2,2);
  cT=prod(i=1,#ellglobalred(E)[5],ellglobalred(E)[5][i][4]); N1=ellcard(E,m1); N2=ellcard(E,m2); Mm=cT*N1*N2;
  R=vector(#P,i,ellmul(E,P[i],Mm));
  B1=Bmat(E,m1,g1,n1,R); B2=Bmat(E,m2,g2,n2,R);
  D1=matdet(B1)%n1; D2=matdet(B2)%n2;
  D12=(B1[1,1]*B2[2,2]+B2[1,1]*B1[2,2]-2*B1[1,2]*B2[1,2])%gg;
  p1=c2prime(Ms,m1,Den); p2=c2prime(Ms,m2,Den);
  [n1,n2,gg, a20%n1, a11%gg, a02%n2, Mm, D1, D12, D2, ellap(E,m1), ellap(E,m2), p1, p2]}
