\\ Mazur-Tate theta element of E at prime m, projected to Z/p^e quotient of (Z/m)^x, then to F_p[x]/(x^(p^e)) via sigma -> 1+x. Returns ord_x and leading coeff.
thetaord(E,m,p,e)={my(Ms=msfromell(E,1), M=Ms[1], xp=Ms[2], n=p^e, g=znprimroot(m), vals=vector(n), den=1, pol, co, ordv=-1);
  if((m-1)%n!=0, error("need m = 1 mod p^e"));
  for(a=1,m-1, my(t=znlog(Mod(a,m),g)%n, v=mseval(M,xp,[a/m,oo])); vals[t+1]+=v);
  den=denominator(content(vals)); if(den%p==0, return([-2,den]));  \\ p in denominator: skip
  \\ polynomial sum_t vals[t+1]*(1+x)^t mod p, truncated mod x^n
  pol=sum(t=0,n-1, vals[t+1]*(1+'x)^t); pol=pol*Mod(1,p);
  co=Vec(pol); \\ highest degree first
  co=vector(n,i,polcoeff(pol,i-1));
  for(i=1,n, if(co[i]!=0, ordv=i-1; break));
  [ordv, if(ordv>=0,co[ordv+1],0), den]}
