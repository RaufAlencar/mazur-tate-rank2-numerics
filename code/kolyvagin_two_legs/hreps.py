# Best Heegner representatives per class using Atkin-Lehner involutions W_Q (Q || N).
# For orientation beta (B == m*beta1 mod 2N) and each class c of disc D m^2, returns
#   key -> (A', B', C', Q, eps)  with  phi(tau_c) = eps * phi(tau')  (tau' root of [A',B',C'], eps = w_Q),
# choosing the representative tau' with the largest imaginary part.
import math, itertools

def divisors_exact(N, pari):
    ps=[int(p) for p in pari(f'factor({N})[,1]')]
    out=[]
    for r in range(len(ps)+1):
        for S in itertools.combinations(ps,r):
            Q=1
            for p in S: Q*=p**int(pari(f'valuation({N},{p})'))
            out.append(Q)
    return out

def wmat(N, Q):
    # W_Q = [Q x, y; N z, Q w], det Q:  Q x w - (N/Q) y z = 1
    if Q==1: return (1,0,0,1)
    M=N//Q
    for z in range(1,10**6):
        g=math.gcd(Q, M*z)
        if g==1:
            # solve Q*x*w - M*z*y = 1 with w=1: Q x - M z y = 1
            y=(-pow(M*z,-1,Q))%Q  # Q x = 1 + M z y
            x=(1+M*z*y)//Q
            assert Q*x*1-M*y*z==1
            return (Q*x, y, N*z, Q*1)
    raise Exception

def transform(F, g):
    # form of tau = g tau', where tau' is root of F; returns primitive [A,B,C] with A>0
    A,B,C=F; a,b,c,d=g
    A2=A*d*d - B*d*c + C*c*c
    B2=-2*A*d*b + B*(d*a+b*c) - 2*C*c*a
    C2=A*b*b - B*a*b + C*a*a
    g0=math.gcd(math.gcd(A2,B2),C2)
    A2//=g0;B2//=g0;C2//=g0
    if A2<0: A2,B2,C2=-A2,-B2,-C2
    return (A2,B2,C2)

def orient_forms(N, Dm, b, K, cop, pari):
    # forms [N k, B, C] with B == b mod 2N, minimal k per class (k <= K, gcd(k,cop)=1)
    S={}
    for k in range(1,K+1):
        if math.gcd(k,cop)!=1: continue
        A=N*k
        B0=b-2*N*((A+b)//(2*N)+1)
        for B in range(B0, A+1, 2*N):
            if B>-A and (B*B-Dm)%(4*A)==0:
                C=(B*B-Dm)//(4*A)
                if math.gcd(math.gcd(A,B),C)==1:
                    key=str(pari.qfbred(pari.Qfb(A,B,C)))
                    if key not in S: S[key]=(A,B,C)
    return S

def best_reps(N, D, m, beta1, cop, pari, K=3000, wsign=None):
    Dm=D*m*m; bm=(m*beta1)%(2*N)
    h=int(pari.qfbclassno(Dm))
    Qs=divisors_exact(N,pari)
    best={}
    for Q in Qs:
        g=wmat(N,Q)
        # orientation b' of tau' such that W_Q tau' has orientation bm: take any form of orientation bm,
        # apply W_Q^{-1} (= W_Q mod Gamma0(N) up to scalar), read its B mod 2N
        F0=(N, bm, (bm*bm-Dm)//(4*N)) if (bm*bm-Dm)%(4*N)==0 else None
        a,b_,c,d=g
        ginv=(d,-b_,-c,a)
        Fp=transform(F0, ginv)
        assert Fp[0]%N==0, (Q,Fp)
        bp=Fp[1]%(2*N)
        S=orient_forms(N,Dm,bp,K,cop,pari)
        eps=1 if Q==1 else wsign[Q]
        for key,Fq in S.items():
            Ft=transform(Fq,g)
            assert Ft[0]%N==0 and (Ft[1]-bm)%(2*N)==0, ('orientation',Q,Fq,Ft)
            k2=str(pari.qfbred(pari.Qfb(*Ft)))
            im=math.sqrt(-Dm)/(2*Fq[0])
            if k2 not in best or best[k2][5]<im: best[k2]=(Fq[0],Fq[1],Fq[2],Q,eps,im)
    if len(best)!=h: raise Exception(f'classes {len(best)} != {h}')
    return best

def al_pairs(N, D, m, beta1, Q, pari, K=300):
    """pairs (F1, Fq): F1 orientation bm (direct), Fq orientation of W_Q^{-1}; same class after W_Q. Smallest A's first."""
    Dm=D*m*m; bm=(m*beta1)%(2*N); g=wmat(N,Q); a,b_,c,d=g
    Fp=transform((N,bm,(bm*bm-Dm)//(4*N)),(d,-b_,-c,a)); bp=Fp[1]%(2*N)
    S1=orient_forms(N,Dm,bm,K,1,pari); Sq=orient_forms(N,Dm,bp,K,1,pari)
    out=[]
    for key,Fq in Sq.items():
        k2=str(pari.qfbred(pari.Qfb(*transform(Fq,g))))
        if k2 in S1: out.append((max(S1[k2][0],Fq[0]),S1[k2],Fq))
    out.sort(); return out
