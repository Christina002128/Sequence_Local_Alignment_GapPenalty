# Comparison of constant, linear, convex, affine gap penalties (local alignment)
# Using script files including constant_gap.py, linear_gap.py, convex_gap.py, affine_gap.py
# for all of the scripts, match score = 20, mismatch score = -0.5 ;
# for constant and linear gaps, gap_penalty = -5 ;
# for convex and affine gaps, open_gap_penalty = -5, extend_gap_penalty = -0.5 ;

# Compare different 5 types of sequences:
# 1. short seqs high similarity dna1 2 EMBOSS Water: match 89/98 (90.8%) gaps 2/98 (2.04%)
# TCTGCTTTGTTATTTCTAGGATTGGCTTCCAGGTCCTCCTCAGCATTGGGGAACTCAGCTTGGCCATCGTGATAGCTCTCACGTCTGGTCTCCTGACAGG
# CTTTGTTATTTCTAGGATTGGCTTCCAGGTCCTCCTCAGCATTGGGCAAATCTGCTCTGGGCCATCGTGATAGCTTTCACGTCCGGTATCCTGACAGG
# 2. short seqs low similarity dna3 4 (EMBOSS Water: match 70/103 (68.0%) gaps 24/103 (23.3%)
# TCTGCTTTGTTATTTCTAGGATTGGCTTCCAGGTCCTCCTCAGCATTGGGGAACTCAGCTTGGCCATCGTGATAGCTCTCACGTCTGGTCTCCTGACAGG
# TTTGTTCATTGGGGAACTTCCAGGTCCTCCTCAGCATTGGGCTAGGATTGGCTTCCATAGCTCTCACGTCTGGTCCCATCGTCCTCCCATTGGGGA
# 3. long seqs high similarity dna5 6  EMBOSS Water: match 553/576 (96.0%)  gaps 5/576 (0.9%)
# NM_001002027.2 (576 bp) 
# XM_008013033.2 
# 4. long seqs low similarity dna7 8  EMBOSS Water: match 254/700 (36.3%)  gaps 365/700 (52.1%)
# NM_001002027.2 (576 bp) 
# JS563783.1 (540bp) 
# 5. long seqs with large countinuous gaps dna9 10  EMBOSS Water: match 1251/3074 (40.7%) Gaps 1451/3074 (47.2%)
# NM_001689.5   2587 bp
# XM_030253440.2  2207bp


import os,time
import numpy as np
dna_1=['dna1','dna3','dna5','dna7','dna9']
dna_2=['dna2','dna4','dna6','dna8','dna10']
count_time=np.zeros(4)
total_time=[]
stats=[]
string=['constant','linear','convex','affine']
for i in range(5):
    start = time.time()
    os.system(f'python constant_gap.py {dna_1[i]} {dna_2[i]} 1 > constant{i}.result')
    end = time.time()
    count_time[0]=(end - start)

    start = time.time()
    os.system(f'python linear_gap.py {dna_1[i]} {dna_2[i]} 1 > linear{i}.result')
    end = time.time()
    count_time[1]=(end - start)

    start = time.time()
    os.system(f'python convex_gap.py {dna_1[i]} {dna_2[i]} 1 > convex{i}.result')
    end = time.time()
    count_time[2]=(end - start)

    start = time.time()
    os.system(f'python affine_gap.py {dna_1[i]} {dna_2[i]} 1 > affine{i}.result')
    end = time.time()
    count_time[3]=(end - start)

    total_time.extend(count_time)
    
    print(f'\nComparing dna',2*i+1,'and dna',2*i+2,'...')
    for x in range(4):
        print(f'\n{string[x]} gaps:')
        print('time:',count_time[x])
        os.system(f'tail -6 {string[x]}{i}.result')
        with open(f'{string[x]}{i}.result','r') as file:
            res=file.readlines()
        align=res[-3].replace('\n', '')
        max_score=float(res[-6].replace('max score: ', '').replace(' ', '').replace('\n', ''))
        # count number of match and mismatch and gap
        match=align.count('|')
        mismatch=align.count('.')
        gap=align.count(' ')
        al_len=len(align)
        # percentage
        per_mismatch=mismatch/al_len
        per_gap=gap/al_len
        identity=round(match/al_len,4)
        stats.append([identity,per_mismatch,per_gap,max_score])

print(total_time)
print(stats)

import pickle
with open('total_time.pkl', 'wb') as f:
    pickle.dump(total_time, f)
with open('stats.pkl', 'wb') as f:
    pickle.dump(stats, f)

os.system('rm *.result')