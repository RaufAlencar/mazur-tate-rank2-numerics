# The runs behind the u_E stress test (note §4) and the 433a1 row; mt_rank_run.py ainvs mmax r
python3 mt_rank_run.py "[0,0,1,-13,18]" 330 2   # 655a1  c_inf=1 tam=2  -> u=4
python3 mt_rank_run.py "[0,1,0,-5,4]"   330 2   # 1132a1 c_inf=1 tam=3  -> u=4
python3 mt_rank_run.py "[0,0,0,-4,1]"   330 2   # 916c1  c_inf=2 tam=3  -> u=8
python3 mt_rank_run.py "[0,1,1,-4,2]"   330 2   # 571b1  c_inf=1 tam=1  -> u=4
python3 mt_rank_run.py "[0,1,0,-25,39]" 300 2   # 1088j1 c_inf=2 tam=4 tors=2 -> u=2
python3 mt_rank_run.py "[0,0,1,-7,6]"   300 3   # 5077a1 rank 3 -> u=-16
python3 mt_rank_run.py "[1,0,0,0,1]"    113 2   # 433a1  c_inf=1 -> u=4 (21 non-trivial levels)
python3 mt_ota_run.py  "[0,1,1,-2,0]"   61 3    # composite levels 389a1 (28 pairs)
python3 mt_ota_run.py  "[1,0,0,0,1]"    61 3    # composite levels 433a1
python3 gsweep.py 400                            # delta = c_L = 1 sweep (3631 points, 453 curves)
