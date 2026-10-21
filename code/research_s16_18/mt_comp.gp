\\ composite level S=m1*m2: theta in Z[G1 x G2], Gi = Z/ni; second-order class and pairing determinants
mtcomp(E,m1,m2,P1,P2,cT)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], g1=znprimroot(m1), g2=znprimroot(m2), n1=(m1-1)/2, n2=(m2-1)/2, S=m1*m2, C=matrix(n1,n2), a00=0,a10=0,a01=0,a20=0,a11=0,a02=0, den, N1,N2,Mm,R1,R2,R12, h1,h2, D20,D02,D11, k1,k2, gg=gcd(n1,n2));
  for(a=1,S-1, if(gcd(a,S)>1,next); my(t1=znlog(Mod(a,m1),g1)%n1, t2=znlog(Mod(a,m2),g2)%n2); C[t1+1,t2+1]+=mseval(M,xp,[a/S,oo]));
  den=denominator(content(C)); C=C*den;
  for(i=0,n1-1,for(j=0,n2-1, my(c=C[i+1,j+1]); if(c==0,next); a00+=c; a10+=c*i; a01+=c*j; a20+=c*binomial(i,2); a02+=c*binomial(j,2); a11+=c*i*j));
  if(a00%n1!=0 || a00%n2!=0, return([-1,"a00"]));
  if(a10%n1!=0, return([-2,"a10",a10%n1])); if(a01%n2!=0, return([-3,"a01",a01%n2]));
  k1=a10/n1; k2=a01/n2;
  my(c20=(a20-k1*binomial(n1,2))%n1, c02=(a02-k2*binomial(n2,2))%n2, c11=a11%gg);
  N1=ellcard(E,m1); N2=ellcard(E,m2); Mm=cT*N1*N2; R1=ellmul(E,P1,Mm); R2=ellmul(E,P2,Mm); R12=elladd(E,R1,R2);
  h1=[hvalu(R1,m1,g1,n1,-1), hvalu(R2,m1,g1,n1,-1), hvalu(R12,m1,g1,n1,-1)];
  h2=[hvalu(R1,m2,g2,n2,-1), hvalu(R2,m2,g2,n2,-1), hvalu(R12,m2,g2,n2,-1)];
  \\ 4*det for pure terms; mixed: 4*(h1_11 h2_22 + h2_11 h1_22 - 2 h1_12 h2_12) with 2h12 = hs - h11 - h22
  D20=(4*h1[1]*h1[2]-(h1[3]-h1[1]-h1[2])^2)%n1; D02=(4*h2[1]*h2[2]-(h2[3]-h2[1]-h2[2])^2)%n2;
  D11=(4*h1[1]*h2[2]+4*h2[1]*h1[2]-2*(h1[3]-h1[1]-h1[2])*(h2[3]-h2[1]-h2[2]))%gg;
  [n1,n2,gg, c20,c02,c11, D20,D02,D11, N1,N2,Mm, den]}
