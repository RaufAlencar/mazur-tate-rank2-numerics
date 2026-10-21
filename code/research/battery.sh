#!/bin/bash
cd /tmp/claude-0/-home-claude/3e963557-071e-5378-94ef-6245b6256453/scratchpad
run(){ echo "=== D=$1 l=$2 (predicted: $5) ==="; python3 kgen2.py $1 $2 $3 $4 20 2>&1 | grep -v Warn | grep -v "^q="; }
run -7 17 200 1500000 "kappa=0"
run -11 17 200 1500000 "kappa!=0 line P1-P2"
run -11 29 260 2500000 "kappa=0"
run -19 29 260 2500000 "kappa!=0 line P1-P2"
run -7 41 320 4000000 "kappa!=0 line P1-P2"
run -11 41 320 4000000 "kappa!=0 line P1-P2"
run -19 41 320 4000000 "kappa!=0 line P1-P2"
run -19 53 380 6000000 "kappa=0"
run -7 59 400 6000000 "kappa!=0 line P2"
run -19 59 400 6000000 "kappa!=0 line P2"
echo BATTERY_DONE
