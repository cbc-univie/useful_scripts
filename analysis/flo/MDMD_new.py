from __future__ import print_function

import argparse
import os
import re
import time
from pathlib import Path

import MDAnalysis
import numpy as np

# from MDAnalysis.newanalysis.correl import correlateParallel
from newanalysis.correl import correlateParallel
from newanalysis.functions import (  # "apr", "rfa", needed for faster calculation of centerOfMassByResidue, dipoleMomentByResidue
    atomsPerResidue,
    centerOfMassByResidue,
    dipoleMomentByResidue,
    residueFirstAtom,
)

# for use with conda env "mda_py3"
# python msd_charged.py [--insidebash]

############
# argparse
###########
parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--insidebash",
    help="use if script is run inside another script",
    action="store_true",
)
parser.add_argument(
    "-s",
    "--short",
    help="non verbose option, no printing of current frame calculation",
    action="store_true",
)
args = parser.parse_args()

#################################################################################################################################
# Trajectory
#################################################################################################################################
firstfile = 1
lastfile = 1
skip = 1

if args.insidebash:
    # automatically tries to find correct psf and dcd files.
    # needed structure: dir1/dir2, dir1/traj/polarizable_nvt
    # file must be executed from within dir2
    # psf needs to be in dir1

    # path = pathlib.Path(__file__).parent.absolute #path for dir of script being run
    path = Path().absolute()  # current working dir
    base = path.parent
    try:
        psf_generator = base.glob("*.psf")
        psf = list(psf_generator)[0]
        print("PSF file:", psf)
    except:
        print("Error, no psf found")
    try:
        dcd_base = base / "traj/"
        f = [f for f in os.listdir(dcd_base) if re.match("nvt_1.dcd", f)]
        # f = [
        #     f
        #     for f in os.listdir(dcd_base)
        #     if re.match(
        #         "[a-zA-Z0-9]+_[a-zA-Z0-9]+_[0-9]+_[a-zA-Z0-9]+_[a-zA-Z0-9]+_[0-9]+_nvt_1.dcd",
        #         f,
        #     )
        # ]
        file_base = f[0][:-5]
        print("DCD base:", dcd_base / file_base)
        dcd = [
            "%s%d.dcd" % (dcd_base / file_base, i)
            for i in range(firstfile, lastfile + 1)
        ]
    except:
        print("Error, no dcd found")

else:
    # check manually for right paths to psf and dcd files!

    # base="/site/raid11a/student3/Florian_Joerg/python_ph_1/"
    base = "/site/raid2/florian/pt/tests/npt_nvt_sd200/"

    psf = base + "im1h_oac_200_im1_hoac_300.psf"
    dcd_base = "traj/im1h_oac_200_im1_hoac_300_nvt_"
    dcd = ["%s%d.dcd" % (base + dcd_base, i) for i in range(firstfile, lastfile + 1)]
    print("PSF file:", psf)
    print("DCD base:", base + dcd_base)

print("Using dcd files from", firstfile, "to", lastfile)
print("Skip is", skip)

u = MDAnalysis.Universe(psf, dcd)
boxl = np.float64(round(u.coord.dimensions[0], 4))
dt = round(u.trajectory.dt, 4)

# n = u.trajectory.numframes//skip  # used integer division (//) instead of (/) -> so dtype is int, otherwise problems with line 42 com_cat=np.zeros(...)
n = int(u.trajectory.n_frames / skip)
if u.trajectory.n_frames % skip != 0:
    n += 1

print("boxl", boxl)

###########################################################################################################################
# molecular species
###########################################################################################################################
sel_cat1 = u.select_atoms("resname IM1H")
ncat1 = sel_cat1.n_residues
mass_cat1 = sel_cat1.masses
charge_cat1 = sel_cat1.charges
apr_cat1 = atomsPerResidue(sel_cat1)
rfa_cat1 = residueFirstAtom(sel_cat1)
com_cat1 = np.zeros((n, ncat1, 3), dtype=np.float64)
mdcat1 = np.zeros((n, 3), dtype=np.float64)
print("Number IM1H   = ", ncat1)

sel_an1 = u.select_atoms("resname OAC")
nan1 = sel_an1.n_residues
mass_an1 = sel_an1.masses
charge_an1 = sel_an1.charges
apr_an1 = atomsPerResidue(sel_an1)
rfa_an1 = residueFirstAtom(sel_an1)
com_an1 = np.zeros((n, nan1, 3), dtype=np.float64)
mdan1 = np.zeros((n, 3), dtype=np.float64)
print("Number OAC   = ", nan1)

sel_cat2 = u.select_atoms("resname IM1")
ncat2 = sel_cat2.n_residues
mass_cat2 = sel_cat2.masses
charge_cat2 = sel_cat2.charges
apr_cat2 = atomsPerResidue(sel_cat2)
rfa_cat2 = residueFirstAtom(sel_cat2)
com_cat2 = np.zeros((n, ncat2, 3), dtype=np.float64)
mdcat2 = np.zeros((n, 3), dtype=np.float64)
print("Number IM1   = ", ncat2)

sel_an2 = u.select_atoms("resname HOAC")
nan2 = sel_an2.n_residues
mass_an2 = sel_an2.masses
charge_an2 = sel_an2.charges
apr_an2 = atomsPerResidue(sel_an2)
rfa_an2 = residueFirstAtom(sel_an2)
com_an2 = np.zeros((n, nan2, 3), dtype=np.float64)
mdan2 = np.zeros((n, 3), dtype=np.float64)
print("Number OAC   = ", nan2)

###########################################################################################################################
# Analysis
###########################################################################################################################
md = np.zeros((n, 3), dtype=np.float64)

ctr = 0
start = time.time()
print("")

for ts in u.trajectory[::skip]:
    print(
        "\033[1AFrame %d of %d" % (ts.frame, u.trajectory.n_frames),
        "\tElapsed time: %.2f hours" % ((time.time() - start) / 3600),
    )

    # efficiently calculate center-of-mass coordinates
    coor_cat1 = np.ascontiguousarray(
        sel_cat1.positions, dtype="double"
    )  # array shape needed for C
    com_cat1 = centerOfMassByResidue(
        sel_cat1, coor=coor_cat1, masses=mass_cat1, apr=apr_cat1, rfa=rfa_cat1
    )
    coor_an1 = np.ascontiguousarray(
        sel_an1.positions, dtype="double"
    )  # array shape needed for C
    com_an1 = centerOfMassByResidue(
        sel_an1, coor=coor_an1, masses=mass_an1, apr=apr_an1, rfa=rfa_an1
    )
    coor_cat2 = np.ascontiguousarray(
        sel_cat2.positions, dtype="double"
    )  # array shape needed for C
    com_cat2 = centerOfMassByResidue(
        sel_cat2, coor=coor_cat2, masses=mass_cat2, apr=apr_cat2, rfa=rfa_cat2
    )
    coor_an2 = np.ascontiguousarray(
        sel_an2.positions, dtype="double"
    )  # array shape needed for C
    com_an2 = centerOfMassByResidue(
        sel_an2, coor=coor_an2, masses=mass_an2, apr=apr_an2, rfa=rfa_an2
    )

    mdcat1[ctr] += np.sum(
        dipoleMomentByResidue(
            sel_cat1,
            coor=coor_cat1,
            charges=charge_cat1,
            masses=mass_cat1,
            com=com_cat1,
            apr=apr_cat1,
            rfa=rfa_cat1,
        ),
        axis=0,
    )
    mdan1[ctr] += np.sum(
        dipoleMomentByResidue(
            sel_an1,
            coor=coor_an1,
            charges=charge_an1,
            masses=mass_an1,
            com=com_an1,
            apr=apr_an1,
            rfa=rfa_an1,
            axis=0,
        ),
        axis=0,
    )
    mdcat2[ctr] += np.sum(
        dipoleMomentByResidue(
            sel_cat2,
            coor=coor_cat2,
            charges=charge_cat2,
            masses=mass_cat2,
            com=com_cat2,
            apr=apr_cat2,
            rfa=rfa_cat2,
        ),
        axis=0,
    )
    mdan2[ctr] += np.sum(
        dipoleMomentByResidue(
            sel_an2,
            coor=coor_an2,
            charges=charge_an2,
            masses=mass_an2,
            com=com_an2,
            apr=apr_an2,
            rfa=rfa_an2,
            axis=0,
        ),
        axis=0,
    )

    # mdcat1[ctr] += np.sum(sel_cat1.dipoleMomentByResidue(coor=coor_cat1,charges=charge_cat1,masses=mass_cat1,com=com_cat1),axis=0)
    # mdan1[ctr]  += np.sum(sel_an1.dipoleMomentByResidue(coor=coor_an1,charges=charge_an1,masses=mass_an1,com=com_an1),axis=0)
    md[ctr] += mdcat1[ctr] + mdan1[ctr] + mdcat2[ctr] + mdan2[ctr]
    # print(mdcat1[ctr])

    ctr += 1

print("calculating <MD(0)MD(t)>")
print(md)
print(mdcat1)
print(np.shape(md))
outfolder = "./mdmd_orig"
Path(outfolder).mkdir(parents=True, exist_ok=True)

# total
md = np.ascontiguousarray(md.T)
md0mdt = np.zeros(n)
correlateParallel(md, md, md0mdt, ltc=1)

f1 = open(f"{outfolder}/md0mdt_all.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdcat1 = np.ascontiguousarray(mdcat1.T)
md0mdt = np.zeros(n)
correlateParallel(mdcat1, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_im1h.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdan1 = np.ascontiguousarray(mdan1.T)
md0mdt = np.zeros(n)
correlateParallel(mdan1, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_oac.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdcat2 = np.ascontiguousarray(mdcat2.T)
md0mdt = np.zeros(n)
correlateParallel(mdcat2, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_im1.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdan2 = np.ascontiguousarray(mdan2.T)
md0mdt = np.zeros(n)
correlateParallel(mdan2, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_hoac.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()
