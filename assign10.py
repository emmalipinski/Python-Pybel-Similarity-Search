#!/usr/bin/env python3

import sys
from openbabel import pybel
import operator
import numpy as np
import statistics as stat

file1 = sys.argv[1]
file2 = sys.argv[2]
thresh = sys.argv[3]
thresh = float(thresh)


#60%
mols1 = list(pybel.readfile('smi', file1))
mols2 = list(pybel.readfile('smi', file2))

def lipinski(mol):
    desc = mol.calcdesc()
    return desc['HBD'] <= 5 and desc['HBA1'] <= 10 and desc['MW'] <= 500 and desc['logP'] <= 5

lipRules = []
for m in mols2:
    lipRules.append(lipinski(m))

molec2 = []
for i in pybel.readfile('smi', file2):
        molec2.append(i.write('can').rstrip())
    
lipMols = operator.itemgetter(*[i for i, mols2 in enumerate(lipRules) if mols2])(molec2)

#print(*lipMols, sep='\n')

#70%
mol1first = mols1[0]
fp = mol1first.calcfp()

tan = []
for m in mols2:
    f = m.calcfp()
    tan.append(f | fp)

lipTan = operator.itemgetter(*[i for i, tan in enumerate(lipRules) if tan])(tan)

#for i in range(len(lipMols)):
#    print(lipMols[i], '{:.4f}'.format(round(lipTan[i],4)))


#80%
tanimoto = []
for i in range(len(lipMols)):
    tanimoto.append([lipMols[i],lipTan[i]])

tanSort = sorted(tanimoto, key=operator.itemgetter(1), reverse=True)

#for i in range(0,10):
#    print(tanSort[i][0], '{:.4f}'.format(round(tanSort[i][1],4)))


#90%
fpAll = []
for i in range(len(mols1)):
    fpAll.append(mols1[i].calcfp())

f = []
for m in mols2:
    f.append(m.calcfp())

avg = []
for j in range(len(f)):
    means = []
    for i in range(len(fpAll)):
        means.append(f[j] | fpAll[i])
    avg.append(stat.mean(means))

lipAvg = operator.itemgetter(*[i for i, avg in enumerate(lipRules) if avg])(avg)

avgSim = []
for i in range(len(lipMols)):
    avgSim.append([lipMols[i],lipAvg[i]])

avgSimSort = sorted(avgSim, key=operator.itemgetter(1), reverse=True)

#for i in range(0,10):
#    print(avgSimSort[i][0], '{:.4f}'.format(round(avgSimSort[i][1],4)))


#100%
minim = []
for j in range(len(f)):
    mins = []
    for i in range(len(fpAll)):
        mins.append(f[j] | fpAll[i])
    minim.append(min(mins))

lipMins = operator.itemgetter(*[i for i, minim in enumerate(lipRules) if minim])(minim)

maxim = []
for j in range(len(f)):
    maxes = []
    for i in range(len(fpAll)):
        maxes.append(f[j] | fpAll[i])
    maxim.append(max(maxes))

lipMax = operator.itemgetter(*[i for i, maxim in enumerate(lipRules) if maxim])(maxim)

final = []
for i in range(len(lipMols)):
    final.append([lipMols[i], lipMins[i], lipAvg[i], lipMax[i]])

finalFilt = [x for x in final if x[1] > thresh]

finalSort = sorted(finalFilt, key=operator.itemgetter(2), reverse=True)

for i in range(0,10):
    print(finalSort[i][0], '{:.4f}'.format(round(finalSort[i][1],4)), '{:.4f}'.format(round(finalSort[i][2],4)), '{:.4f}'.format(round(finalSort[i][3],4)))









