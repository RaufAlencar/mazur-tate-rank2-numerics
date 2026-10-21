import cypari2,sys,time
pari=cypari2.Pari(); pari.allocatemem(2*10**9)
import glob
rows=[l.split() for f in sorted(glob.glob('data/allbsd.*')) for l in open(f)]
sel=[r for r in rows if r[4]=='2' and int(sys.argv[2])<int(r[0])<=int(sys.argv[1])]
print(len(sel),'rank-2 curves (all in isogeny classes)',flush=True)
bad=0;t=time.time()
for r in sel:
    E=pari.ellinit(pari(r[3])); rk=pari.ellrank(E,0); lo,hi=int(rk[0]),int(rk[1])
    t2=int(pari.elltors(E)[0]); dim2=len([x for x in pari.elltors(E)[1] if int(x)%2==0])
    flag='' if hi==2 else ' <<< Sel_2 EXCESS (upper=%d)'%hi
    if hi!=2: bad+=1
    print(r[0]+r[1]+r[2],'rank bounds',[lo,hi],'dimE[2]=',dim2,flag,flush=True)
print('done',len(sel),'curves, excess Selmer in',bad,'curves,',round(time.time()-t),'s')
