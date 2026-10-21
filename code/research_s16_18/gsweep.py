import cypari2,sys; pari=cypari2.Pari(); pari.allocatemem(3*10**9); pari('read("gsweep.gp")')
from collections import Counter
cnt=Counter(); bad=[]; skipped=0
lines=[l.split() for l in open('data/allcurves.00000-09999') if int(l.split()[0])<=int(sys.argv[1]) and int(l.split()[4])>=1]
for N,iso,num,ainvs,rank,tors in lines:
    pari(f'E=ellinit({ainvs}); G=ellrank(E)[4]; cT=prod(i=1,#ellglobalred(E)[5],ellglobalred(E)[5][i][4]);')
    for i in range(int(pari('#G'))):
        for k in range(1,9):
            try: r=pari(f'gc(E, ellmul(E,G[{i+1}],cT*{k}))')
            except Exception: skipped+=1; continue
            g,c,W=int(r[0]),int(r[1]),int(r[2]); cnt[(g,c)]+=1
            if g!=1 or c!=1: bad.append((N+iso+num,ainvs,i,k,g,c,W%2)); print('EXC',bad[-1],flush=True)
print(len(lines),'curves; (g,c) counts',dict(cnt),'skipped',skipped,flush=True); print('exceptions:',len(bad))
