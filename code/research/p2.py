import cypari2,sys
pari=cypari2.Pari(); pari.allocatemem(2*10**9)
def analyze(ainvs,label,primes,n=10,verbose=False):
    E=pari.ellinit(pari(ainvs)); gens=pari.ellsaturation(E,pari.ellrank(E)[3],1000)
    if len(gens)!=2: return None
    out=[]
    for p in primes:
        N=int(pari.ellglobalred(E)[0])
        if N%p==0 or int(pari.ellap(E,p))%p==0: continue
        A,B=pari.ellpadicheightmatrix(E,p,n,gens); s2=pari.ellpadics2(E,p,n)
        H=A-s2*B; R=pari.ellpadicregulator(E,p,n,gens)
        ok = pari.valuation(pari.matdet(H)-R,p) >= n-3
        card=int(pari.ellcard(pari.ellinit(E,p))); L=pari.truncate(pari.ellformallog(E,40,'t'))
        def plog(P):
            Q=pari.ellmul(E,P,card); t=-Q[0]/Q[1]+pari(f'O({p}^{n+2})'); return pari.subst(L,'t',t)/card
        l1,l2=plog(gens[0]),plog(gens[1])
        # P0 = l1*P2 - l2*P1 ; h(P0,P0)= l1^2 H22 - 2 l1 l2 H12 + l2^2 H11
        hP0 = l1**2*H[1,1]-2*l1*l2*H[0,1]+l2**2*H[0,0]
        Nval = hP0/R
        vN=int(pari.valuation(Nval,p)) if Nval!=0 else 'INF'
        out.append((p,vN,int(pari.valuation(R,p)),int(pari.valuation(l1,p)),int(pari.valuation(l2,p)),ok))
        if verbose: print(label,'p',p,'det(H)==R?',ok,' N=h(P0,P0)/R =',Nval)
    return out
if __name__=='__main__':
    analyze('[0,1,1,-2,0]','389a1',[5,7,11,13],verbose=True)

if __name__=='__main__' and len(sys.argv)>1:
    rows=[l.split() for l in open('data/allbsd.00000-09999')]
    sel=[r for r in rows if r[4]=='2' and r[2]=='1' and int(r[0])<=int(sys.argv[1])]
    from collections import Counter
    vals=Counter(); okc=Counter(); zero=[]
    for r in sel:
        res=analyze(r[3],r[0]+r[1],[3,5,7,11,13],n=8)
        if not res: continue
        for p,vN,vR,v1,v2,ok in res:
            vals[vN]+=1; okc[ok]+=1
            if vN=='INF' or vN>=4: zero.append((r[0]+r[1],p,vN))
        print(r[0]+r[1],[(p,vN,vR,v1,v2) for p,vN,vR,v1,v2,ok in res],flush=True)
    print('det(H)==R check:',dict(okc)); print('valuation of N=h(P0,P0)/R:',dict(vals)); print('high/zero:',zero)
