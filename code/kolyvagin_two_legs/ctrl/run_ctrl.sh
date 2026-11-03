source ~/bsdvenv/bin/activate
cd "$(dirname "$0")"
python kgen3.py "[0,1,1,-2,0]" -7 5 3 200 2000000 20 > ctrl_389a1_D7_l5.log 2>&1
