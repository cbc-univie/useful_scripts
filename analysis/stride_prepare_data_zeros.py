import numpy as np
import os
import itertools
import glob
import re

# Alternative version for stride_prepare_data.py with np.zeros
# probably less efficient than the old version due to explicit for loops,
# but easier to adapt and better at dealing with overlapping regions

#path = "/site/raid4/andras/sandbox/students/maria/3k6s_BPA_GEN/"
# Due to the nature of the STRIDE output, you should always explicitly
# check your psf to see how segments and amino acids are numbered
path = "."
files = sorted(glob.glob(f'{os.getcwd()}/stride*.dat'),
               key=lambda x: int(re.findall('\d+', os.path.basename(x))[0]))
arrays = [np.genfromtxt(f"{i}", usecols=(0,1,2,3,5)) for i in files] 
frame = int(arrays[0][-1][0])+1 # NOTE: if last frame has no defined
                                #sec. structure, this will be false

last_residue = 963
seg0 = 597
seg1 = 57

for pos, array in enumerate(arrays):
    array[:,0]+=(pos*frame)
    for index, residue in enumerate(array):
        if residue[1] == 0:
            residue[2] += seg0-seg1
            residue[3] += seg0-seg1

data_old = np.concatenate(arrays)
ind = np.lexsort((data_old[:,2], -data_old[:,1], data_old[:,0]))
data = data_old[ind]

np.savetxt("data.txt", data, fmt='%i')

final_values = np.zeros((frame*len(files)*last_residue,3), dtype='short')
for i in range(frame*len(files)):
    for j in range(last_residue):
        final_values[i*last_residue+j][0] = i
        final_values[i*last_residue+j][1] = j+1

for i in data:
    frame    = int(i[0])
    aa_begin = int(i[2])
    aa_end   = int(i[3])+1
    sec_str  = int(i[4])
    for j in range(aa_begin, aa_end):
        final_values[last_residue*frame+j-1][2] = sec_str

np.savetxt("final.txt", final_values, fmt='%i')
