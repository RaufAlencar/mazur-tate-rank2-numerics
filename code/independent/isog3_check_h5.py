# Conferencia independente (registro) do teste H5 a partir dos .jsonl brutos da sessao mt-estatistica.
# Subconjunto primario: m = 1 mod 3 (ja no arquivo), 3 | n, 3 nao divide M = c_Tam*#E(F_m); curvas sem 3-torcao; >= 10 niveis.
import json, sys, math, collections
D = r"C:\Users\Administrativo\Desktop\BSD_Estatistica_MT\dados"
def rates(fn):
    per = collections.defaultdict(list); tors3 = set()
    for line in open(f"{D}\{fn}", encoding="utf-8"):
        o = json.loads(line)
        if o["tors"] % 3 == 0: tors3.add(o["label"]); continue
        M = o["ctam"] * o["Np"]
        if o["n"] % 3 or M % 3 == 0: continue
        per[o["label"]].append(1 if o["det"] % 3 == 0 else 0)
    r = [sum(v)/len(v) for v in per.values() if len(v) >= 10]
    return r, len(per), len(tors3)
def ms(r):
    m = sum(r)/len(r); s = math.sqrt(sum((x-m)**2 for x in r)/(len(r)-1)/len(r)); return m, s
for rk, model in (("r2", 1/3), ("r3", 29/81)):
    T, nT, t3 = rates(f"isog3_{rk}_T.jsonl"); C, nC, c3 = rates(f"isog3_{rk}_C.jsonl")
    mT, sT = ms(T); mC, sC = ms(C)
    print(f"{rk}: isog3 {len(T)} curvas (>=10 niveis; {t3} com 3-torcao excluidas) taxa {mT:.4f} +- {sT:.4f} | controle {len(C)} taxa {mC:.4f} +- {sC:.4f} | modelo {model:.4f}")
    print(f"     z(isog3-modelo) = {(mT-model)/sT:.2f}   z(isog3-controle) = {(mT-mC)/math.hypot(sT,sC):.2f}")
