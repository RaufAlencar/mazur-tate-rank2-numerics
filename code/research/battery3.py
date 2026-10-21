import cypari2, subprocess, sys, json
pari=cypari2.Pari(); pari.allocatemem(10**9)
jobs=[]
def add(ainvs,p,Ds,lmax):
    E=pari.ellinit(pari(ainvs)); N=int(pari.ellglobalred(E)[0])
    for D in Ds:
        for l in pari.primes(200):
            l=int(l)
            if l<5 or l>lmax or N%l==0 or (l+1)%p or int(pari.kronecker(D,l))!=-1 or int(pari.ellap(E,l))%p: continue
            deg=l+1; PREC=150 if deg<=20 else (220 if deg<=36 else (300 if deg<=48 else 380))
            jobs.append((ainvs,D,l,p,PREC))
add('[0,1,1,-2,0]',5,[-7,-11],20)          # 389a1, p=5, l=19
add('[0,1,1,-4,2]',3,[-7,-8,-19,-67],60)   # 571b1
add('[0,-1,1,-2,0]',3,[-7,-11,-19],60)     # 709a1
print(len(jobs),'jobs:',[(j[0],j[1],j[2],j[3]) for j in jobs],flush=True)
for ainvs,D,l,p,PREC in jobs:
    print(f'=== {ainvs} D={D} l={l} p={p} PREC={PREC} ===',flush=True)
    r=subprocess.run(['python3','kgen3.py',ainvs,str(D),str(l),str(p),str(PREC),'6000000','20'],capture_output=True,text=True)
    for line in r.stdout.split('\n')+r.stderr.split('\n'):
        if ('PREDICTION' in line or 'FINAL' in line or 'Error' in line or 'Exception' in line) and 'Warn' not in line: print(line,flush=True)
print('BATTERY3_DONE',flush=True)
