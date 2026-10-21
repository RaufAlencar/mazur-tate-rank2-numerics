#!/bin/bash
cd /tmp/claude-0/-home-claude/3e963557-071e-5378-94ef-6245b6256453/scratchpad
run(){ echo "=== D=$1 l=$2 (predicted: $5) ==="; python3 kgen2.py $1 $2 $3 $4 20 2>&1 | grep -v Warn | grep -v "^q="; }
run -7 59 400 6000000 "kappa!=0 line P2"
run -19 59 400 6000000 "kappa!=0 line P2"
run -4 59 400 6000000 "kappa=0 (D=-4, units 4: degree may be 30)"
echo BATTERY2_DONE
