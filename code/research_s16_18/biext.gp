\\ Biextension computation: f = (x - x(2R)) / (y - lam*x - c), tangent line at R; globally primitive normalization F;
\\ psi_m(F) = reduction of F mod m = A_{2R} / (E*c)  ;  returns [ (E c / A_{2R}) mod m , (t/d) mod m, (d/t) mod m ]
biext(E,R,m)={my(x=R[1],y=R[2],a1=E.a1,a2=E.a2,a3=E.a3,a4=E.a4,lam,c,R2,A2,d2,den,Ec,u,t,d,tu,du);
  lam=(3*x^2+2*a2*x+a4-a1*y)/(2*y+a1*x+a3); c=y-lam*x;
  R2=ellmul(E,R,2); A2=numerator(R2[1]); d2=sqrtint(denominator(R2[1]));
  den=lcm(denominator(lam),denominator(c));  \\ E
  Ec=den*c;  \\ integer
  if(valuation(Ec,m)!=0 || valuation(A2,m)!=0, print("warning: unexpected valuations ",valuation(Ec,m)," ",valuation(A2,m)));
  u=Mod(Ec,m)/Mod(A2,m);
  t=-x/y; d=sqrtint(denominator(x)); tu=t/m^valuation(t,m); du=d/m^valuation(d,m);
  [u, Mod(numerator(tu),m)/Mod(denominator(tu),m)*Mod(du,m)^-1, Mod(du,m)*(Mod(numerator(tu),m)/Mod(denominator(tu),m))^-1, valuation(d,m), valuation(den,m)]}
