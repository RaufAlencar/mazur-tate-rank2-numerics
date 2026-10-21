import cypari2, math, sys
pari=cypari2.Pari(); pari.allocatemem(6*10**9); PREC=int(sys.argv[2]); pari.default('realprecision',PREC)
N=389; D=-20; l=int(sys.argv[1]); Dl=l*l*D; E=pari.ellinit([0,1,1,-2,0])
h=int(pari.qfbclassno(Dl)); beta=[b for b in range(2*N) if (b*b-Dl)%(4*N)==0][0]
cls={}; k=1
while len(cls)<h and k<5000:
    A=N*k
    for B in range(-A+1,A+1):
        if (B-beta)%(2*N) or (B*B-Dl)%(4*A): continue
        C=(B*B-Dl)//(4*A)
        if math.gcd(math.gcd(A,abs(B)),C)!=1: continue
        key=str(pari.qfbred(pari.Qfb(A,B,C)))
        if key not in cls or A<cls[key][0]: cls[key]=(A,B,C)
    k+=1
print('l=',l,'h=',h,'classes found',len(cls),'maxA',max(v[0] for v in cls.values()))
ymin=math.sqrt(-Dl)/(2*max(v[0] for v in cls.values())); NT=int(PREC*math.log(10)/(2*math.pi*ymin))+50
an=pari.ellan(E,NT+10)
def z_of(A,B,C):
    tau=(-B+pari.sqrt(pari(Dl)))/(2*A); y=float(pari.imag(tau)); nterms=int(PREC*math.log(10)/(2*math.pi*y))+50
    q=pari.exp(2*pari.Pi()*pari('I')*tau); s=pari(0); qn=pari(1)
    for n in range(1,nterms+1): qn*=q; s+=an[n-1]*qn/n
    return s
xs=[pari.ellztopoint(E,z_of(*v))[0] for v in cls.values()]
T=pari('T'); Pk=pari(1)
for x in xs: Pk*=(T-x)
coef=pari.Vec(Pk); mags=[float(pari.log(abs(c)+pari(10)**-PREC))/math.log(10) for c in coef]
print('log10|coef| max',round(max(mags),1))
bound=10**max(5,int(PREC/2-max(mags)-5)); ok=True; dens=[]
for c in coef[1:]:
    a=pari.real(c); b=pari.imag(c)/pari.sqrt(pari(5)); ra=pari.bestappr(a,bound); rb=pari.bestappr(b,bound)
    if float(max(abs(a-ra),abs(b-rb)))>10**(-(PREC/2-max(mags)-10)): ok=False
    dens.append(math.log10(max(1,int(pari.denominator(ra)))))
print('recognized?',ok,' log10 denominators:',[round(x,1) for x in dens[:8]])
# alternative hypothesis: coefficients in Q(sqrt(Dl))? same field. Try: is the polynomial over Q of degree 2h with real coefs?
import sys; sys.exit()
P2=pari(1)
for x in xs: P2*=(T-x)*(T-pari.conj(x))
c2=pari.Vec(P2); ok2=True; d2=[]
for c in c2[1:]:
    a=pari.real(c); ra=pari.bestappr(a,bound)
    if float(abs(a-ra))>10**(-(PREC/2-2*max(mags)-10)): ok2=False
    d2.append(math.log10(max(1,int(pari.denominator(ra)))))
print('degree-2h real polynomial recognized over Q?',ok2,[round(x,1) for x in d2[:8]])
