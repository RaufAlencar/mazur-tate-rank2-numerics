# Independent prediction: ker(red_l1) ∩ ker(red_l2) in E(Q)/3 (saturated ellrank basis), p=3.
import cypari2, itertools, json, sys
pari=cypari2.Pari(); pari.allocatemem(10**9); pari('read("divp.gp")')
CASES=[('11642a1','[1,-1,0,-16,28]',-7,5,17),('5077a1','[0,0,1,-7,6]',-163,11,23),('18097a1','[1,1,1,-10,6]',-163,11,29)]
out=[]
for name,a,D,l1,l2 in CASES:
    pari(f'E=ellinit({a}); G=ellsaturation(E,ellrank(E)[4],1000)')
    def ker(l): return set(v for v in itertools.product(range(3),repeat=3) if v!=(0,0,0) and int(pari(f'my(R=elladd(E,elladd(E,ellmul(E,G[1],{v[0]}),ellmul(E,G[2],{v[1]})),ellmul(E,G[3],{v[2]}))); if(R==[0]||denominator(R[1])%{l}==0,1,inpE(ellinit(E,{l}),R,3))')))
    I=sorted(ker(l1)&ker(l2))
    if len(I)==2: pred={'tipo':'reta ou zero','reta':[list(v) for v in I]}
    else: pred={'tipo':'zero forcado','dim_intersecao':{8:2,26:3}.get(len(I),'?'),'intersecao':[list(v) for v in I]}
    out.append({'curva':name,'ainvs':a,'D':D,'l1':l1,'l2':l2,'p':3,'base':str(pari('G')),'ker_l1':sorted(map(list,ker(l1))),'ker_l2':sorted(map(list,ker(l2))),'previsao_kappa_n':pred,
                'previsao_kappa_l1':'zero forcado','previsao_kappa_l2':'zero forcado'})
json.dump(out,open('prereg/previsoes.json','w'),indent=1)
print(json.dumps([(o['curva'],o['previsao_kappa_n']) for o in out]))
