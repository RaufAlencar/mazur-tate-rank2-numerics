gc(E,R)={my(x=R[1],y=R[2],A=numerator(x),d=sqrtint(denominator(x)),B=numerator(y),W,Lam,R2,A2,N4,g,c);
  W=2*B+E.a1*A*d+E.a3*d^3; Lam=3*A^2+2*E.a2*A*d^2+E.a4*d^4-E.a1*B*d;
  R2=ellmul(E,R,2); A2=numerator(R2[1]); N4=A^4-E.b4*A^2*d^4-2*E.b6*A*d^6-E.b8*d^8; g=sqrtint(N4/A2);
  c=gcd(gcd(d^3*W,d^2*Lam),B*W-Lam*A); [g,c,W]}
