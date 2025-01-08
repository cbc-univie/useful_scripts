import numpy as np
import sys
from newanalysis.diffusion import msdCOM, unfoldTraj, msdMJ, msdMJdecomp, msdMJcrossterms
from pathlib import Path
import h5py
sys.path.append("/home/florian/pythonfiles")
from helpers import msd_com 

infile = "msd_1-100_skip5.hdf5"

hdf5_file = h5py.File(infile, "r")

firstfile = hdf5_file.attrs["firstfile"]
lastfile = hdf5_file.attrs["lastfile"]
skip = hdf5_file.attrs["skip"]
n = hdf5_file.attrs["n"]
dt = hdf5_file.attrs["dt"] # skip value is excluded hence indexing of array by time is dt*skip [ps], =original timestep = dt/(skip*dcd_save_freq)
boxl = hdf5_file.attrs["boxl"]

print(f"Reading file {infile}")
print(f"Used file from {firstfile} to {lastfile} with skip {skip}")
print(f"{boxl = }, {dt = }, {n = }")
print(f"Distance between two entries is dt*skip = {dt*skip} ps")

resn_cat = "IM1H"
resn_an = "OAC"
resn_cat1 = "IM1"
resn_an1 = "HOAC"


#folded_com_cat = hdf5_file["folded_IM1H"][:]
#folded_com_an = hdf5_file["folded_OAC"][:]
com_cat = hdf5_file["unfolded_im1h"][:]
com_an = hdf5_file["unfolded_oac"][:]
com_cat1 = hdf5_file["unfolded_im1"][:]
com_an1 = hdf5_file["unfolded_hoac"][:]

hdf5_file.close()

#print(com_cat.shape, com_an.shape)

#with_time = False
#if with_time:
timestart = float(sys.argv[1]) #ps
timeend = float(sys.argv[2]) #ps
indexstart = int(timestart/(dt*skip))
indexend = int((timeend/(dt*skip))-1)
if indexend >= len(com_cat):
    print(f"Indexend {indexend} was greater or equal than array com_cat {len(com_cat)}!")
    print(f"Indexend set to maximum value {len(com_cat)-1}")
    indexend = len(com_cat)-1
#else:
#    filestart = 1
#    fileend = 10
#    #frames per file = len(com_an) * skip / lastfile
#    frames_per_file = len(com_an)/lastfile
#    indexstart =  int(filestart*frames_per_file-frames_per_file)
#    indexend =  int(fileend*frames_per_file)
print(f"Input: Calculate msd from time {timestart} ps to {timeend} ps")
print(f"Array is indexed from {indexstart} to {indexend}")

com_an = com_an[indexstart:indexend]
com_cat = com_cat[indexstart:indexend]
com_an1 = com_an1[indexstart:indexend]
com_cat1 = com_cat1[indexstart:indexend]

n = len(com_an)

time_axis = np.arange(0, n*skip*dt, skip*dt)



print("calculating msd ..")
msd_cat = msd_com(com_cat)
msd_an  = msd_com(com_an )
msd_cat1 = msd_com(com_cat1)
msd_an1  = msd_com(com_an1 )

Path("msd").mkdir(parents=True, exist_ok=True)
np.savetxt(f"msd/{resn_cat.lower()}_msd_{int(timestart)}-{int(timeend)}_skip{skip}.dat", np.c_[time_axis, msd_cat])
np.savetxt(f"msd/{resn_an.lower()}_msd_{int(timestart)}-{int(timeend)}_skip{skip}.dat",  np.c_[time_axis, msd_an])
np.savetxt(f"msd/{resn_cat1.lower()}_msd_{int(timestart)}-{int(timeend)}_skip{skip}.dat", np.c_[time_axis, msd_cat1])
np.savetxt(f"msd/{resn_an1.lower()}_msd_{int(timestart)}-{int(timeend)}_skip{skip}.dat",  np.c_[time_axis, msd_an1])

print("calculating msdMJ ..")

outfolder = "./msdmj"
Path(outfolder).mkdir(parents=True, exist_ok=True)

msdmj = msdMJ(com_cat,com_an)
timeseries_msdMJ_cat = msdMJdecomp(com_cat, com_an, com_cat,  1)
timeseries_msdMJ_ani = msdMJdecomp(com_cat, com_an, com_an,  -1)

timeseries_msdMJ_catcat = msdMJcrossterms(com_cat, com_cat, 1,  1)
timeseries_msdMJ_aniani = msdMJcrossterms(com_an,  com_an, -1, -1)
timeseries_msdMJ_catani = msdMJcrossterms(com_cat, com_an,  1, -1)

np.savetxt(f"{outfolder}/msdMJ_{int(timestart)}-{int(timeend)}_total.dat",        np.c_[time_axis, msdmj]) #Change these file names to save to different files
np.savetxt(f"{outfolder}/msdMJ_{int(timestart)}-{int(timeend)}_total_cation.dat", np.c_[time_axis, timeseries_msdMJ_cat])
np.savetxt(f"{outfolder}/msdMJ_{int(timestart)}-{int(timeend)}_total_anion.dat",  np.c_[time_axis, timeseries_msdMJ_ani])
np.savetxt(f"{outfolder}/msdMJ_{int(timestart)}-{int(timeend)}_cation_cation.dat", np.c_[time_axis, timeseries_msdMJ_catcat])
np.savetxt(f"{outfolder}/msdMJ_{int(timestart)}-{int(timeend)}_anion_anion.dat",   np.c_[time_axis, timeseries_msdMJ_aniani])
np.savetxt(f"{outfolder}/msdMJ_{int(timestart)}-{int(timeend)}_cation_anion.dat",  np.c_[time_axis, timeseries_msdMJ_catani])


quit()

