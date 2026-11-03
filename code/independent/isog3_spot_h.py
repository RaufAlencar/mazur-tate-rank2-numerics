# Amostra de linhas brutas (isog3 posto 2 e 3): recalcula det H mod n com mt_indep.gp (pair_det), geradores do ellrank saturados.
import cypari2, json, random, pathlib
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
here = pathlib.Path(__file__).resolve().parent
pari('read("%s")' % (here.parent / "mt_indep.gp").as_posix())
D = "/mnt/c/Users/Administrativo/Desktop/BSD_Estatistica_MT/dados/"
random.seed(1711)
import re, glob
AG = {}
for f in glob.glob("/mnt/c/Users/Administrativo/Desktop/BSD_Estatistica_MT/ecdata_raw/allgens.*"):
    for line in open(f):
        t = line.split()
        if len(t) < 5: continue
        lab = t[0] + t[1] + t[2]; r = int(t[4])
        pts = re.findall(r"\[(-?\d+):(-?\d+):(\d+)\]", line)[:r]
        AG[lab] = (t[3], r, "[" + ",".join(f"[{x}/{z},{y}/{z}]" for x,y,z in pts) + "]")
for fn, k in (("isog3_r2_T.jsonl", 16), ("isog3_r2_C.jsonl", 8), ("isog3_r3_T.jsonl", 5)):
    rows = [json.loads(l) for l in open(D + fn)]
    rows = [o for o in rows if o["tors"] % 3]
    labs = sorted({o["label"] for o in rows}); pick = random.sample(labs, k)
    ok = bad = 0
    for lab in pick:
        rs = [o for o in rows if o["label"] == lab and o["n"] % 3 == 0 and (o["ctam"]*o["Np"]) % 3 and o["ctam"]*o["Np"] < 1500 and o["n"] >= 30]; o = random.choice(rs) if rs else None
        if o is None: continue
        if lab not in AG: print(lab, "fora dos arquivos allgens; pulo"); continue
        ai, r, gens = AG[lab]
        E = pari.ellinit(pari(ai)); g = pari(gens)
        m, n = o["m"], o["n"]; M = o["ctam"] * o["Np"]
        pari('pair_det_py(E,P,M,m)=pair_det(E,P,M,m,znprimroot(m))[1]')
        d = pari('pair_det_py')(E, g, M, m)
        d = int(pari.lift(pari.Mod(d, n))); mine0 = d % 3 == 0; theirs0 = o["det"] % 3 == 0
        same = (d - o["det"]) % n == 0 or (d + o["det"]) % n == 0
        print(f"{lab} M={M} m={m} n={n} det_registro={o['det']} det_meu={d} igual(+-)={same} zero_mod3: {theirs0}/{mine0}", flush=True)
        ok += same; bad += not same
    print(fn, "iguais", ok, "diferentes", bad)
