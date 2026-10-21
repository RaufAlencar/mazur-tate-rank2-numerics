import cypari2,sys; pari=cypari2.Pari(); pari.allocatemem(8*10**9); pari('read("mt_lead.gp")')
from collections import Counter
def test(label,ainvs,mmax):
    pari(f'E=ellinit({ainvs}); G=ellsaturation(E,ellrank(E)[4],1000); P1=G[1]; P2=G[2]; N0=ellglobalred(E)[1];')
    N0=int(pari('N0')); ratios=Counter(); rows=[]; skipped=0
    for m in [int(q) for q in pari.primes(600) if int(q)%9==1 and N0%int(q) and int(q)<=mmax]:
        co=pari(f'thetacoefs(E,{m},znprimroot({m}),3,2)')
        if co==0: continue
        ordv=next((i for i in range(9) if int(pari.lift(co[i]))!=0),-1); c2=int(pari.lift(co[2]))
        reg=pari(f'mtreg(E,{m},znprimroot({m}),P1,P2)')
        if len(reg)==1: skipped+=1; continue
        R=int(reg[0]); ratio='both0' if (c2==0 and R==0) else ('MISMATCH' if (c2==0)!=(R==0) else (c2*pow(R,-1,3))%3)
        ratios[ratio]+=1; rows.append((m,ordv,c2,R,ratio))
    print(f'{label}: tested {sum(ratios.values())} primes (skipped {skipped} with 3|#E(F_m)); outcomes {dict(ratios)}',flush=True)
    bad=[r for r in rows if r[4]=='MISMATCH' or (r[4] not in ('both0','MISMATCH') and r[4]!=rows[0][4] and rows[0][4]!='both0')]
    print('   rows:',rows,flush=True)
test('389a1',sys.argv[1] if len(sys.argv)>1 else '[0,1,1,-2,0]', int(sys.argv[2]) if len(sys.argv)>2 else 3000)
