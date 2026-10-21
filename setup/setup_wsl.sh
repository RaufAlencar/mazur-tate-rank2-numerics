#!/bin/bash
# Prepara o ambiente Python + PARI/GP (cypari2) no WSL para reproduzir os resultados.
set -e
rm -rf ~/bsdvenv
if [ ! -x ~/.local/bin/uv ]; then
  curl -LsSf https://astral.sh/uv/install.sh | sh >/dev/null 2>&1
fi
UV=~/.local/bin/uv
$UV --version
$UV venv --python 3.12 ~/bsdvenv
$UV pip install --python ~/bsdvenv/bin/python cypari2==2.2.4 numpy
~/bsdvenv/bin/python - <<'EOF'
import cypari2
pari = cypari2.Pari()
print("PARI", pari.version())
print("posto analitico de 389a1:", pari('ellanalyticrank(ellinit([0,1,1,-2,0]))'))
EOF
