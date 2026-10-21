chk(E,R)={my(x=R[1],y=R[2],A=numerator(x),d=sqrtint(denominator(x)),B=numerator(y),a1=E.a1,a2=E.a2,a3=E.a3,a4=E.a4,W,Lam,R2,A2,d2,N4,g,cL);
  if(denominator(y)!=d^3, return("y denominator not d^3"));
  W=2*B+a1*A*d+a3*d^3; Lam=3*A^2+2*a2*A*d^2+a4*d^4-a1*B*d;
  R2=ellmul(E,R,2); A2=numerator(R2[1]); d2=sqrtint(denominator(R2[1]));
  N4=A^4-E.b4*A^2*d^4-2*E.b6*A*d^6-E.b8*d^8;
  g=sqrtint(N4/A2);
  cL=gcd(gcd(d^3*W,d^2*Lam),B*W-Lam*A);
  [N4==g^2*A2, d2==d*W/g, g^2==cL, g, cL]}
