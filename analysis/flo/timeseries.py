import numpy as np
import MDAnalysis
from newanalysis.diffusion import msdCOM, unfoldTraj, msdMJ, msdMJdecomp, msdMJcrossterms
from newanalysis.correl import correlateParallel
from newanalysis.functions import atomsPerResidue, residueFirstAtom # "apr", "rfa", needed for faster calculation of centerOfMassByResidue, dipoleMomentByResidue
from newanalysis.functions import centerOfMassByResidue, dipoleMomentByResidue
import os, sys
import time
from datetime import timedelta
from pathlib import Path
import re
import argparse
import h5py

sys.path.append("/home/florian/pythonfiles")
from helpers import msd_com, msd_mj

#for use with conda env "mda2.0"
#python msd_charged.py [--insidebash]

############
# argparse
###########
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--insidebash", help="use if script is run inside another script", action="store_true")
parser.add_argument("-s", "--short", help="non verbose option, no printing of current frame calculation", action="store_true")
args = parser.parse_args()

#################################################################################################################################
# Trajectory
#################################################################################################################################
firstfile=1
lastfile=100
skip=1

residues = ["IM1H", "OAC", "IM1", "HOAC"]

if args.insidebash:
    # automatically tries to find correct psf and dcd files.
    # needed structure: dir1/dir2, dir1/traj/polarizable_nvt
    # file must be executed from within dir2
    # psf needs to be in dir1

    #path = pathlib.Path(__file__).parent.absolute #path for dir of script being run
    path = Path().absolute() #current working dir
    base = path.parent
    try:
        psf_generator = base.glob("*.psf")
        psf = list(psf_generator)[0]
        print("PSF file:", psf)
    except:
        print("Error, no psf found")
    try:
        dcd_base = base / "traj/"
        f = [f for f in os.listdir(dcd_base) if re.match("[a-zA-Z0-9]+_[a-zA-Z0-9]+_[0-9]+_[a-zA-Z0-9]+_[a-zA-Z0-9]+_[0-9]+_nvt_1.dcd", f)]
        file_base = f[0][:-5]
        print("DCD base:", dcd_base/file_base)
        dcd = ["%s%d.dcd" % (dcd_base/file_base,i) for i in range(firstfile,lastfile+1)]
    except:
        print("Error, no dcd found")

else:    
    # check manually for right paths to psf and dcd files!
    base="/site/raid2/florian/pt/tests/npt_nvt_sd200/"
    psf=base+"im1h_oac_200_im1_hoac_300.psf"
    dcd_base = "traj/im1h_oac_200_im1_hoac_300_nvt_"
    dcd = ["%s%d.dcd" % (base+dcd_base,i) for i in range(firstfile,lastfile+1)]
    print("PSF file:", psf)
    print("DCD base:", base+dcd_base)

print(f"Using dcd files from {firstfile} to {lastfile} with skip {skip}")

u=MDAnalysis.Universe(psf,dcd)
boxl = np.float64(round(u.coord.dimensions[0],4))
dt=round(u.trajectory.dt,4)

#n = u.trajectory.numframes//skip  # used integer division (//) instead of (/) -> so dtype is int, otherwise problems with line 42 com_cat=np.zeros(...)
n = int(u.trajectory.n_frames/skip)
if u.trajectory.n_frames%skip != 0:
    n+=1
print(f"Boxlength is {boxl}")
#################################################################################################################################
# molecular species
#################################################################################################################################
coms = {}
sel = {}
n_res = {}
masses = {}
charges = {}
apr = {}
rfa = {}
coms = {}
md = {}
for residue in residues:
    sel[residue] = u.select_atoms(f"resname {residue}")
    n_res[residue] = sel[residue].n_residues
    masses[residue] = sel[residue].masses
    charges[residue] = sel[residue].charges
    apr[residue] = atomsPerResidue(sel[residue])
    rfa[residue] = residueFirstAtom(sel[residue])
    coms[residue] = np.zeros((n,n_res[residue],3), dtype=np.float64)
    md[residue] = np.zeros((n,3),dtype=np.float64)
    print(f"Number of {residue}: {n_res[residue]}")

md_total = np.zeros((n,3),dtype=np.float64)
#################################################################################################################################
# Running through trajectory
#################################################################################################################################
ctr=0

start=time.time()
print("")

for ts in u.trajectory[::skip]:
    if not args.short:
        print("\033[1AFrame %d of %d" % (ts.frame,len(u.trajectory)), "\tElapsed time: %.2f hours" % ((time.time()-start)/3600))#"    Est. time left: %s" % (str(timedelta(seconds=(time.time()-start)/ts.frame * (u.trajectory.numframes-ts.frame)))[:7]))
    
    for residue in residues:
        coms[residue][ctr] = sel[residue].center_of_mass(compound="residues")
        coor = np.ascontiguousarray(sel[residue].positions, dtype='double') #array shape needed for C
        md[residue][ctr] += np.sum(dipoleMomentByResidue(sel[residue], coor=coor,charges=charges[residue],masses=masses[residue],com=coms[residue][ctr], apr=apr[residue], rfa=rfa[residue],axis=0), axis=0)
        md_total[ctr]     += md[residue][ctr]

    ctr+=1

hdf5_file = h5py.File(f"msd_{firstfile}-{lastfile}_skip{skip}.hdf5", "w")
hdf5_file.attrs["n"] = n
hdf5_file.attrs["dt"] = dt
hdf5_file.attrs["boxl"] = boxl
hdf5_file.attrs["skip"] = skip
hdf5_file.attrs["firstfile"] = firstfile
hdf5_file.attrs["lastfile"] = lastfile

#hdf5_file.create_group("folded_trajs")
for residue in residues:
    hdf5_file.create_dataset(f"folded_{residue.lower()}", data=coms[residue])

#################################################################################################################################
# Post-processing
#################################################################################################################################
print("unfolding coordinates ..")
for residue in residues:
    unfoldTraj(coms[residue],boxl)

#hdf5_file.create_group("unfolded_trajs")
for residue in residues:
    hdf5_file.create_dataset(f"unfolded_{residue.lower()}", data=coms[residue])

hdf5_file.close()

time_axis = np.arange(0, n*skip*dt, skip*dt) # eigentlich ab dt bis lastvalue included

timestart = 0
timeend = int(n*skip*dt)

print("calculating msd ..")
msdfolder = "./msd"
Path(msdfolder).mkdir(parents=True, exist_ok=True)

for residue in residues:
    #msd_old = msdCOM(coms[residue])
    #np.savetxt(f"{msdfolder}/old_msd_{residue.lower()}_{firstfile}-{lastfile}_{skip}.dat", np.c_[time_axis, msd_old])
    msd = msd_com(coms[residue])
    #np.savetxt(f"{msdfolder}/msd_{residue.lower()}_{firstfile}-{lastfile}_{skip}.dat", np.c_[time_axis, msd])
    np.savetxt(f"{msdfolder}/msd_{residue.lower()}_{timestart}-{timeend}.dat", np.c_[time_axis, msd])

print("calculating msdMJ ..")
msdmjfolder = "./msdmj"
Path(msdmjfolder).mkdir(parents=True, exist_ok=True)

#msdmj_old = msdMJ(coms["IM1H"],coms["OAC"])
msdmj = msd_mj((coms["IM1H"], 1), (coms["OAC"], -1))

timeseries_msdMJ_cat = msdMJdecomp(coms["IM1H"], coms["OAC"], coms["IM1H"],  1)
timeseries_msdMJ_ani = msdMJdecomp(coms["IM1H"], coms["OAC"], coms["OAC"],  -1)

timeseries_msdMJ_catcat = msdMJcrossterms(coms["IM1H"], coms["IM1H"], 1,  1)
timeseries_msdMJ_aniani = msdMJcrossterms(coms["OAC"],  coms["OAC"], -1, -1)
timeseries_msdMJ_catani = msdMJcrossterms(coms["IM1H"], coms["OAC"],  1, -1)

#np.savetxt(f"{msdmjfolder}/old_msdMJ_{firstfile}-{lastfile}_{skip}_total.dat",        np.c_[time_axis, msdmj_old]) #Change these file names to save to different files
np.savetxt(f"{msdmjfolder}/msdMJ_{timestart}-{timeend}_total.dat",        np.c_[time_axis, msdmj]) 
np.savetxt(f"{msdmjfolder}/msdMJ_{timestart}-{timeend}_total_cation.dat", np.c_[time_axis, timeseries_msdMJ_cat])
np.savetxt(f"{msdmjfolder}/msdMJ_{timestart}-{timeend}_total_anion.dat",  np.c_[time_axis, timeseries_msdMJ_ani])
np.savetxt(f"{msdmjfolder}/msdMJ_{timestart}-{timeend}_cation_cation.dat", np.c_[time_axis, timeseries_msdMJ_catcat])
np.savetxt(f"{msdmjfolder}/msdMJ_{timestart}-{timeend}_anion_anion.dat",   np.c_[time_axis, timeseries_msdMJ_aniani])
np.savetxt(f"{msdmjfolder}/msdMJ_{timestart}-{timeend}_cation_anion.dat",  np.c_[time_axis, timeseries_msdMJ_catani])

print("calculating <MD(0)MD(t)>")
mdmdfolder = "./mdmd"
Path(mdmdfolder).mkdir(parents=True, exist_ok=True)

# total
md_total = np.ascontiguousarray(md_total.T)
md0mdt = np.zeros(n)
correlateParallel(md_total,md_total,md0mdt,ltc=1)
np.savetxt(f"{mdmdfolder}/md0mdt_all.dat", np.c_[time_axis, md0mdt])

for residue in residues:
    md[residue] = np.ascontiguousarray(md[residue].T)
    md0mdt = np.zeros(n)
    correlateParallel(md[residue],md_total,md0mdt,ltc=1)
    np.savetxt(f"{mdmdfolder}/md0mdt_{residue.lower()}.dat", np.c_[time_axis, md0mdt])


quit()

#np.save(f"timeseries/timeseries_{resn_an}_{firstfile}-{lastfile}_{skip}", com_an, allow_pickle=False)
#np.save(f"timeseries/timeseries_{resn_cat}_{firstfile}-{lastfile}_{skip}", com_cat, allow_pickle=False)
##to open timeseries for later analysis:
#with open(f"timeseries_{resn_an}_{firstfile}-{lastfile}_{skip}.npy",  "rb") as f:
#    com_an2 = np.load(f, allow_pickle=False)
#with open(f"timeseries_{resn_cat}_{firstfile}-{lastfile}_{skip}.npy", "rb") as f:
#    com_cat2 = np.load(f, allow_pickle=False)
#np.save(f"timeseries/unfolded_timeseries_{resn_an}_{firstfile}-{lastfile}_{skip}", com_an, allow_pickle=False)
#np.save(f"timeseries/unfolded_timeseries_{resn_cat}_{firstfile}-{lastfile}_{skip}", com_cat, allow_pickle=False)
