# plot the proformance of different gap penalties
import pickle
import numpy as np
import matplotlib.pylab as plt

# performace results from comparison.py
with open('total_time.pkl', 'rb') as f:
    time0 = pickle.load(f)
with open('stats.pkl', 'rb') as f:
    stats = pickle.load(f)


identity0=[stats[i][0] for i in range(20)]
mismatch0=[stats[i][1] for i in range(20)]
gap0=[stats[i][2] for i in range(20)]
score0=[stats[i][3] for i in range(20)]

# make them into matrices, each sub list is the performance of each gap penalty for 5 types of sequence-pair
time=[[time0[(4*k)+i] for k in range(5)]for i in range(4)]
identity=[[identity0[(4*k)+i] for k in range(5)]for i in range(4)]
mismatch=[[mismatch0[(4*k)+i] for k in range(5)]for i in range(4)]
gap=[[gap0[(4*k)+i] for k in range(5)]for i in range(4)]
score=[[score0[(4*k)+i] for k in range(5)]for i in range(4)]

# results of each sequence-pairs from EMBOSS Water (https://www.ebi.ac.uk/Tools/psa/emboss_water/)
emboss_match=[89,70,553,254,1251]
emboss_gap_number=[2,24,5,365,1451]
emboss_align_len=[98,103,576,700,3074]
emboss_identity=[emboss_match[i]/emboss_align_len[i] for i in range(5)]
emboss_mismatch=[1-(emboss_match[i] + emboss_gap_number[i])/emboss_align_len[i] for i in range(5)]
emboss_gap=[emboss_gap_number[i]/emboss_align_len[i] for i in range(5)]
string=['constant','linear','convex','affine']
dna=['short_high','short_low','long_high','long_low','large_gap']

# plot
a,b=2,2
plt.figure(figsize=(12,10))

# identity plot
plt.subplot(a,b,1)
for i in range(4):
    plt.plot(dna,identity[i])
plt.plot(dna,emboss_identity,'-.')
plt.title('Identity')
plt.ylabel('identity%')
plt.xlabel('sequence-pair')
plt.legend(string+['EMBOSS_Water'])

# mismatches plot
plt.subplot(a,b,2)
for i in range(4):
    plt.plot(dna,mismatch[i])
plt.plot(dna,emboss_mismatch,'-.')
plt.title('Mismatches')
plt.ylabel('mismatches%')
plt.xlabel('sequence-pair')
plt.legend(string+['EMBOSS_Water'])

# number of gaps plot
plt.subplot(a,b,3)
for i in range(4):
    plt.plot(dna,gap[i])
plt.plot(dna,emboss_gap,'-.')
plt.title('Gaps')
plt.ylabel('gaps%')
plt.xlabel('sequence-pair')
plt.legend(string+['EMBOSS_Water'])

# maximum score plot
plt.subplot(a,b,4)
for i in range(4):
    plt.plot(dna,score[i])
plt.title('Maximum Score')
plt.ylabel('scores')
plt.xlabel('sequence-pair')
plt.legend(string)


plt.savefig('Align_Performance.png')

# Time of Running plot
plt.figure(figsize=(5,10))
for i in range(4):
    plt.plot(dna,time[i])
plt.title('Time of Running')
plt.ylabel('time(s)')
plt.xlabel('sequence-pair')
plt.legend(string)

plt.savefig('Time_Performance.png')
