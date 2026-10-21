# PASSO 1 do teste cego: registra a constante u_E PREVISTA pela fórmula da nota (§4),
#   u_E = (-2)^r * c_inf / #tors^2,
# para curvas que a nota não testou. Usa só invariantes da curva; não calcula nenhum theta nem pareamento.
import cypari2, datetime, json
pari = cypari2.Pari()
pari.allocatemem(10**9)

CURVES = {  # Cremona: a-invariantes (tabela allbsd)
    "563a1": "[1,1,1,-15,16]",
    "643a1": "[1,0,0,-4,3]",
    "664a1": "[0,0,0,-7,10]",
    "681c1": "[0,-1,1,0,2]",
    "707a1": "[0,1,1,-12,12]",
    "997b1": "[0,-1,1,-5,-3]",
    "11197a1": "[1,-1,1,-6,0]",
}
out = {"registrado_em": datetime.datetime.now().isoformat(timespec="seconds"),
       "formula": "u_E = (-2)^r * c_inf / tors^2", "curvas": {}}
for label, a in CURVES.items():
    E = pari(f"ellinit({a})")
    r = int(pari.ellanalyticrank(E)[0])
    cinf = 2 if pari(f"ellinit({a}).disc") > 0 else 1
    tors = int(pari.elltors(E)[0])
    ctam = int(pari.ellglobalred(E)[2])
    gens = pari.ellsaturation(E, pari.ellrank(E)[3], 1000)
    assert len(gens) == r, (label, gens)
    u = (-2) ** r * cinf // tors ** 2
    out["curvas"][label] = {"ainvs": a, "r": r, "c_inf": cinf, "tors": tors, "c_Tam": ctam,
                            "base": str(gens), "u_E_previsto": u}
    print(label, out["curvas"][label])
json.dump(out, open("blind_predictions.json", "w"), indent=1, ensure_ascii=False)
print("gravado blind_predictions.json em", out["registrado_em"])
