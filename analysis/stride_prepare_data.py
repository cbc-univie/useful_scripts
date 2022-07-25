import numpy as np
import itertools
import glob
import re

#path = "/site/raid1/student5/flora_schoefbeck/PROA_A/charmm-gui-3413037598/analysis/stride/"
path = "."
first, last = 0,200
arrays = [np.loadtxt(f"{path}/stride{i}.dat", usecols=(0,2,3,5)) for i in range(first,last)] 
for pos, array in enumerate(arrays):
    array[:,0]+=((first+pos)*5000)



data_old = np.concatenate(arrays)
ind = np.lexsort((data_old[:,1], data_old[:,0]))
data = data_old[ind]

np.savetxt("data.txt", data)

final_values = []
number = 1
frame = 5000
last_residue = 24
for pos, val in enumerate(data):
    if frame != int(val[0]):
        if number <= last_residue:
            for i in range(number, last_residue+1):
                final_values.append([val[0], i, 0])
        frame = int(val[0])
        number = 1
    if number < int(val[1]):
        for i in range(int(number),int(val[1])):
            final_values.append([val[0], i, 0])
        number = int(val[1])
    if int(val[1]) <= number < int(val[2]):
        for i in range(number, int(val[2])):
            final_values.append([val[0], i, val[3]])
        number=int(val[2])
if number <= last_residue:
    for i in range(number, last_residue+1):
        final_values.append([val[0], i, 0])
final_values = np.asarray(final_values)
print(final_values.shape)
np.savetxt("final.txt", final_values)

