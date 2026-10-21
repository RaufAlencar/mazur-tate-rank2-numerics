import cypari2, math, sys, pickle
pari=cypari2.Pari(); pari.allocatemem(8*10**9); PREC=int(sys.argv[1]); pari.default('realprecision',PREC)
N=389; D=-20; l=17; Dl=l*l*D; E=pari.ellinit([0,1,1,-2,0])
d=pickle.load(open('koly17.pkl','rb')); keys,cls,G=d['keys'],d['cls'],d['G']
maxA=max(v[0] for v in cls.values()); ymin=math.sqrt(-Dl)/(2*maxA)
NT=int(PREC*math.log(10)/(2*math.pi*ymin))+50; print('terms needed',NT); an=pari.ellan(E,NT+10)
def z_of(k_):
    A,B,C=cls[k_]; tau=(-B+pari.sqrt(pari(Dl)))/(2*A); y=float(pari.imag(tau))
    nterms=int(PREC*math.log(10)/(2*math.pi*y))+50
    q=pari.exp(2*pari.Pi()*pari('I')*tau); s=pari(0); qn=pari(1)
    for n in range(1,nterms+1):
        qn*=q; s+=an[n-1]*qn/n
    return s
Z={k_:z_of(k_) for k_ in keys}; X={k_:pari.ellztopoint(E,Z[k_]) for k_ in keys}
xs=[X[k_][0] for k_ in keys]
print('|x| range:',min(float(abs(x)) for x in xs),max(float(abs(x)) for x in xs))
T=pari('T'); Pk=pari(1)
for x in xs: Pk*=(T-x)
coef=pari.Vec(Pk)
mags=[float(pari.log(abs(c)+pari(10)**-PREC))/math.log(10) for c in coef]
print('log10|coef| :',[round(m,1) for m in mags])
pickle.dump({'xs':[(str(pari.real(x)),str(pari.imag(x))) for x in xs],'ys':[(str(pari.real(X[k][1])),str(pari.imag(X[k][1]))) for k in keys],'coef':[(str(pari.real(c)),str(pari.imag(c))) for c in coef],'PREC':PREC},open(f'koly17_x_{PREC}.pkl','wb'))
# try recognition with denominators bounded by 10^(PREC/2 - maxmag - 5)
import itertools
ok=True; dens=[]
for c in coef[1:]:
    a=pari.real(c); b=pari.imag(c)/pari.sqrt(pari(5)); bound=10**max(5,int(PREC/2-max(mags)-5))
    ra=pari.bestappr(a,bound); rb=pari.bestappr(b,bound)
    err=float(max(abs(a-ra),abs(b-rb))); dens.append(int(pari.denominator(ra)))
    if err>10**(-(PREC/2-max(mags)-10)): ok=False
print('recognized?',ok,' denominators (first 6):',dens[:6], ' max log10 den', max(math.log10(x) for x in dens))
