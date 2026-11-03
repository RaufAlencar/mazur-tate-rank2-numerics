import json, pathlib, re
here = pathlib.Path(__file__).resolve().parent
lv = json.loads((here.parent.parent / "results/independent/site_levels_389a1_5077a1.json").read_text())
src = (here / "src.html").read_text(encoding="utf-8")
out = src.replace("__LEVELS__", json.dumps(lv, separators=(",", ":")))
bad = re.findall("[–—]", out)
assert not bad, f"travessao encontrado: {len(bad)}"
(here.parent / "index.html").write_text(out, encoding="utf-8")
print("ok", len(out))
