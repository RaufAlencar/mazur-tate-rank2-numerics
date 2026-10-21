tunit(Q,m)={my(t=-Q[1]/Q[2],d=sqrtint(denominator(Q[1])),r=(t/m^valuation(t,m))/(d/m^valuation(d,m))); Mod(numerator(r),m)/Mod(denominator(r),m)}
offd(E,R,S,m)={my(lam,c,den,L,A3,RS,u,pol);
  lam=(S[2]-R[2])/(S[1]-R[1]); c=R[2]-lam*R[1];
  den=lcm(denominator(lam),denominator(c));
  L=[den, -den*lam, -den*c]; L=L/content(L);
  RS=elladd(E,R,S); A3=numerator(RS[1]);
  if(L[1]%m!=0 || L[2]%m!=0, return([-1, L[1]%m, L[2]%m]));
  u=Mod(-A3,m)/Mod(L[3],m);
  pol=tunit(RS,m)/(tunit(R,m)*tunit(S,m));
  [if(u^-1==pol || u^-1==-pol,1,0), u^-1, pol]}
