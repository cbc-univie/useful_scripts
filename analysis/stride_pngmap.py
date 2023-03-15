import numpy as np
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

uni_hexes = ["0063A6", "A71C49", "DD4814", "F6A800", "94C154", "11897A", "666666"]
uni_rgbs = []

for i in range(len(uni_hexes)):
    uni_rgbs.append((int(uni_hexes[i][0:2], 16)/255., int(uni_hexes[i][2:4], 16)/255., int(uni_hexes[i][4:6], 16)/255.))

uni_tuple = tuple(uni_rgbs)        

frames = np.genfromtxt("final.txt", usecols=(0,))
nums   = np.genfromtxt("final.txt", usecols=(1,))
secstr = np.genfromtxt("final.txt", usecols=(2,))

next_frame = frames[0]
num_aa = -1
while next_frame == frames[0]:
    num_aa += 1
    next_frame = frames[num_aa]


x = np.reshape(frames, (int(len(frames)/num_aa),num_aa))
y = np.reshape(nums, (int(len(nums)/num_aa),num_aa))
z = np.reshape(secstr, (int(len(secstr)/num_aa),num_aa))

fig, ax = plt.subplots()
plt.subplot(111)
plt.xlabel("Time / ns")
plt.ylabel("Amino Acid")
plt.yticks([i for i in range(100,num_aa+2,100)], [str(i) for i in range(100,num_aa+2,100)])
plt.xticks([i for i in range(0,20000,2000)],  [str(i//100) for i in range(0,20000,2000)])

plt.ylim((0.5,num_aa+0.5))
cmap = mpl.colors.ListedColormap(uni_tuple)
cs = plt.pcolormesh(x,y,z, shading='nearest', cmap=cmap)
cb = fig.colorbar(cs)
cb.set_ticks([3/7,9/7,15/7,21/7,27/7,33/7,39/7,45/7])
cb.set_ticklabels(["Coil", "Helix", "3_10 Helix", "Pi Helix", "Strand", "Bridge", "Turn"])
plt.savefig("heatmap.png")
#plt.close()
#print("DONE WITH FIG "+str(i+1))
