import numpy as np
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import sys

frames = np.loadtxt("final.txt", usecols=(0,))
nums   = np.loadtxt("final.txt", usecols=(1,))
secstr = np.loadtxt("final.txt", usecols=(2,))

num_aa = 24

x = np.reshape(frames, (int(len(frames)/num_aa),num_aa))
y = np.reshape(nums, (int(len(nums)/num_aa),num_aa))
z = np.reshape(secstr, (int(len(secstr)/num_aa),num_aa))

fig, ax = plt.subplots()
plt.subplot(111)
plt.xlabel("Time / ps")
plt.ylabel("Amino Acid")
plt.yticks([i for i in range(1,num_aa+1)], [str(i) for i in range(1,num_aa+1)])
plt.ylim((0.5,num_aa+0.5))
cbarcolors = plt.cm.Accent.colors[0:7]
cmap = mpl.colors.ListedColormap(cbarcolors)
cs = plt.pcolormesh(x,y,z, shading='nearest', cmap=cmap)
cb = fig.colorbar(cs)
cb.set_ticks([3/7,9/7,15/7,21/7,27/7,33/7,39/7,45/7])
cb.set_ticklabels(["Coil", "Helix", "3_10 Helix", "Pi Helix", "Strand", "Bridge", "Turn"])
plt.savefig("heatmap.png")
#plt.close()
#print("DONE WITH FIG "+str(i+1))
