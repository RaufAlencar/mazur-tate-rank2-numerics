import cypari2, math, sys, pickle
pari=cypari2.Pari(); pari.allocatemem(6*10**9); PREC=int(sys.argv[1]); pari.default('realprecision',PREC)
N=389; D=-20; l=17; Dl=l*l*D; E=pari.ellinit([0,1,1,-2,0])
d=pickle.load(open('koly17.pkl','rb')); keys,cls,G=d['keys'],d['cls'],d['G']
an=pari.ellan(E,3000000)
def z_of(k_):
    A,B,C=cls[k_]; tau=(-B+pari.sqrt(pari(Dl)))/(2*A); y=float(pari.imag(tau))
    nterms=int(PREC*math.log(10)/(2*math.pi*y))+50; assert nterms<=len(an),(nterms,A)
    q=pari.exp(2*pari.Pi()*pari('I')*tau); s=pari(0); qn=pari(1)
    for n in range(1,nterms+1):
        qn*=q; s+=an[n-1]*qn/n
    return s
Z={k_:z_of(k_) for k_ in keys}
X={k_:pari.ellztopoint(E,Z[k_]) for k_ in keys}
xs=[X[k_][0] for k_ in keys]
# polynomial over K: coefficients should be in Q(sqrt-5)
Pk=pari(1); T=pari('T')
for x in xs: Pk*=(T-x)
coef=pari.Vec(Pk)
print('degree',len(coef)-1)
# check whether root set closed under complex conjugation
conjclosed=all(min(abs(pari.conj(x)-y) for y in xs) < pari(10)**(-PREC//2) for x in xs[:5])
print('x-coordinates closed under conjugation (=> Q(x) degree 36 over Q?)', conjclosed)
rec=[]; ok=True; maxden=0
for i,c in enumerate(coef):
    a=pari.real(c); b=pari.imag(c)/pari.sqrt(pari(5))   # c = a + b*sqrt(-5)
    ra=pari.bestappr(a,10**(PREC//3)); rb=pari.bestappr(b,10**(PREC//3))
    err=max(abs(a-ra),abs(b-rb))
    if err>pari(10)**(-(PREC//3)): ok=False
    maxden=max(maxden,int(pari.denominator(ra)),int(pari.denominator(rb)))
    rec.append((ra,rb))
print('all coefficients recognized in Q(sqrt-5)?',ok,' max denominator',maxden)
print('first coeffs:',rec[:3]); print('leading-ish coeffs:',rec[-4:])
pickle.dump({'rec':[(str(a),str(b)) for a,b in rec],'ok':ok},open('koly17_poly.pkl','wb'))
