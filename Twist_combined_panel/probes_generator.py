#tested on python3.9
#dir = sys.argv[1] eg. /Users/chyiyin/Desktop/crow_project/00_probe/
#basename of file = sys.argv[2] eg. snp_panel_combined_184k_0based

import sys
import pandas as pd
import random

dir = sys.argv[1]
name = sys.argv[2]
f = open(f'{dir}{name}_final.txt', "r")
output = open(f'{dir}{name}_final_assignedallele.txt', 'w')
base = ["A", "T", "C", "G"]
for x in f:
    line = x.split()
    ref = line[2]
    alt = line[3]
    if ref != "C" and alt != "C":
        assign = list(set(base)-set(ref)-set(alt))
        final = line + [random.choice(assign)]
        output.write(' '.join(final) +  '\n')
    elif ref == "C" or alt == "C":
        assign = list(set(base)-set(ref)-set(alt)-set("T"))
        final = line + [random.choice(assign)]
        output.write(' '.join(final) + '\n')
    else:
        print("problem!")
output.close()

fa=pd.read_csv(f'{dir}{name}_upstream.fa.bed', header=None, sep='\t')
fb=pd.read_csv(f'{dir}{name}_final_assignedallele.txt', header=None, sep=' ')
fc=pd.read_csv(f'{dir}{name}_downstream.fa.bed', header=None, sep='\t')
df=fa[1]+fb[4]+fc[1]
result = pd.concat([fb[0],fb[1]-40,fb[1]+40,df],axis="columns")
result.to_csv(f'{dir}{name}_probeset.txt',index=False,header= False, sep='\t')

file = open(f'{dir}{name}_probeset.txt', "r")
fastaprobe = open(f'{dir}{name}_probeset.fa', 'w')
for y in file:
    probeline = y.split()
    scaffoldname = str(f'>{probeline[0]}:{probeline[1]}-{probeline[2]}')
    fastaprobe.write(scaffoldname + '\n')
    fastaprobe.write(probeline[3] + '\n')
fastaprobe.close()
