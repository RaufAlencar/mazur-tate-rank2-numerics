\\ common-denominator recognition: r = vector of reals sharing a denominator d (|d| <~ 10^L, L < (P-30)*t/(t+1));
\\ LLL on t of them gives d; then every r_k*d must be within 10^-20 of an integer (small extra denominators,
\\ up to 10^6, are absorbed). Returns [d, integers] or [0, k] on failure at coordinate k.
reccommon(r,P,t)={my(m=#r,C=10^(P-30),M,v,d,dd=1,nums);
  t=min(t,m);
  M=matrix(t+1,t+1,i,j, if(j==1, if(i==1,1,round(C*r[i-1])), if(i==j,C,0)));
  v=(M*qflll(M))[,1]; d=abs(v[1]); if(d==0,return([0,0]));
  for(k=1,m, my(x=r[k]*d*dd, y=round(x));
    if(abs(x-y)>10^(-20), y=bestappr(x,10^6); if(type(y)!="t_INT" && type(y)!="t_FRAC", return([0,k])); if(abs(x-y)*denominator(y)>10^(-20), return([0,k])); dd*=denominator(y)));
  d*=dd; nums=vector(m,k,round(r[k]*d)); for(k=1,m, if(abs(r[k]*d-nums[k])>10^(-20),return([0,-k]))); [d,nums]}
