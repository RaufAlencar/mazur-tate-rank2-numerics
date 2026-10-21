pairf(E,R,m)={my(x=R[1],y=R[2],A=numerator(x),d=sqrtint(denominator(x)),B=numerator(y),a1=E.a1,a2=E.a2,a3=E.a3,a4=E.a4,W,Lam,R2,A2,N4,g,cL,val,td);
  W=2*B+a1*A*d+a3*d^3; Lam=3*A^2+2*a2*A*d^2+a4*d^4-a1*B*d;
  R2=ellmul(E,R,2); A2=numerator(R2[1]);
  N4=A^4-E.b4*A^2*d^4-2*E.b6*A*d^6-E.b8*d^8; g=sqrtint(N4/A2);
  cL=gcd(gcd(d^3*W,d^2*Lam),B*W-Lam*A);
  val=Mod(-A^3,m)/Mod(A2*cL,m);          \\ <R,R> from the biextension, in F_m^x (up to sign)
  td=(Mod(-A,m)/Mod(B,m))^2;              \\ (t/d)^2 = (A/B)^2
  [val==td || val==-td, val, td, Mod(A,m)^-1, g, cL]}
