# Gera um PDF tipografado da nota a partir do markdown exportado do Claude Doc:
# blocos ```latex viram fórmulas de destaque renderizadas pelo KaTeX; o resto é o texto do Doc sem alterações,
import re, sys, html, subprocess, pathlib, markdown

SRC = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "note_rev105.md")
OUT_HTML = SRC.with_suffix(".html")
OUT_PDF = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "Alencar_MazurTate_rank2_3_v0.1.pdf")

md = SRC.read_text(encoding="utf-8")
# Texto da revisão congelada sem nenhuma edição (rev 111 em diante).

blocks = []
def stash(m):
    blocks.append(m.group(1).strip())
    return f"\n\nKATEXBLOCK{len(blocks)-1}\n\n"
md = re.sub(r"```latex\n(.*?)\n```", stash, md, flags=re.S)

body = markdown.markdown(md, extensions=["tables"])
for i, tex in enumerate(blocks):
    body = body.replace(f"<p>KATEXBLOCK{i}</p>", f'<div class="math">{html.escape(tex)}</div>')
assert "KATEXBLOCK" not in body

page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Numerical verification of the refined Mazur–Tate conjecture in ranks 2 and 3</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<style>
@page {{ size: A4; margin: 22mm 20mm; }}
body {{ font-family: "Cambria", "Georgia", serif; font-size: 11.5pt; line-height: 1.5; color: #111; }}
h1 {{ font-size: 17pt; line-height: 1.25; margin: 0 0 6pt; }}
h2 {{ font-size: 12.5pt; margin: 16pt 0 6pt; }}
p {{ margin: 0 0 7pt; text-align: justify; }}
table {{ border-collapse: collapse; font-size: 9pt; margin: 8pt 0; width: 100%; }}
th, td {{ border: 0.5pt solid #999; padding: 2pt 4pt; text-align: left; vertical-align: top; }}
th {{ background: #f0f0f0; }}
.math {{ margin: 8pt 0; overflow-x: hidden; }}
code {{ font-size: 9pt; }}
li {{ margin-bottom: 4pt; text-align: justify; }}
</style></head><body>
{body}
<script>
document.querySelectorAll('.math').forEach(function(el) {{
  katex.render(el.textContent, el, {{displayMode: true, throwOnError: true, strict: false}});
}});
document.title = document.title;
</script></body></html>"""
OUT_HTML.write_text(page, encoding="utf-8")

import time
if OUT_PDF.exists(): OUT_PDF.unlink()
edge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
subprocess.run([edge, "--headless=new", f"--user-data-dir={pathlib.Path('edge_profile').resolve()}", "--disable-gpu", "--no-pdf-header-footer", "--virtual-time-budget=15000",
                f"--print-to-pdf={OUT_PDF.resolve()}", OUT_HTML.resolve().as_uri()], check=True)
last = -1
for _ in range(120):  # o Edge pode devolver o controle antes de terminar de gravar
    time.sleep(1)
    if OUT_PDF.exists() and OUT_PDF.stat().st_size == last and last > 0: break
    last = OUT_PDF.stat().st_size if OUT_PDF.exists() else -1
print("ok", OUT_PDF, OUT_PDF.stat().st_size, "bytes;", len(blocks), "fórmulas")
