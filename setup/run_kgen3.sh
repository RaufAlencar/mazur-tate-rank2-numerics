#!/bin/bash
# Section 5 of the note: Kolyvagin class for 389a1, K=Q(sqrt -7), l=5, p=3. Runs in a scratch copy.
set -e
SRC="$(cd "$(dirname "$0")/../code/research" && pwd)"
WORK="$(mktemp -d)"
cp "$SRC"/*.py "$SRC"/*.gp "$WORK"/
cd "$WORK"
~/bsdvenv/bin/python kgen3.py "[0,1,1,-2,0]" -7 5 3 120 300000 20
