#!/bin/bash
# Section 8 (Exp. 1): strong BSD on optimal rank-2 curves with N <= 1500. Run from the repository root.
set -e
cd "$(dirname "$0")/../code/research"
mkdir -p data
[ -f data/allbsd.00000-09999 ] || curl -sSfL -o data/allbsd.00000-09999 \
  https://raw.githubusercontent.com/JohnCremona/ecdata/master/allbsd/allbsd.00000-09999
~/bsdvenv/bin/python bsd2.py 1500
