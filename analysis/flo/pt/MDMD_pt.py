from __future__ import print_function

import argparse
import os
import re
import sys
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
from newanalysis.helpers import dipByResidue

sys.path.append("/home/florian/pythonfiles")
from cond_helpers import Charges

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
lastfile = 100
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
        file_base = f[0][:-5]
        print("DCD base:", dcd_base / file_base)
        dcd = [
            "%s%d.dcd" % (dcd_base / file_base, i)
            for i in range(firstfile, lastfile + 1)
        ]
    except:
        print("Error, no dcd found")
    try:
        charges_base = base / "out/"
        file_base = "charge_changes"
        print("Charge base:", charges_base / file_base)
        charge_files = [
            f"{charges_base}/{file_base}_{i}.out"
            for i in range(firstfile, lastfile + 1)
        ]
    except FileNotFoundError as err:
        print("Error, charges file not found")
        print(err)

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

print("Frames to process:", n)
print("boxl", boxl)

charges = Charges()
charges.from_file(charge_files)
charges.info()
charge_steps = charges.data.keys()
assert len(charge_steps) == len(
    u.trajectory
), f"{len(charge_steps)=}, {len(u.trajectory)=}"

points = int(charges.time_between_updates / (dt * skip))
print(f"{points=}, {charges.steps_between_updates=}")

###########################################################################################################################
# molecular species
###########################################################################################################################
residues = ["IM1H", "OAC", "IM1", "HOAC"]

coms = {}
sel = {}
n_res = {}
n_at = {}
apr = {}
idx_firstlast = {}
idx = 0
for residue in residues:
    sel[residue] = u.select_atoms(f"resname {residue}")
    n_res[residue] = sel[residue].n_residues
    n_at[residue] = sel[residue].n_atoms
    apr[residue] = sel[residue].n_atoms // n_res[residue]
    idx_firstlast[residue] = [idx]
    idx += n_res[residue]
    idx_firstlast[residue].append(idx)
    # masses[residue] = sel[residue].masses
    # charges[residue] = sel[residue].charges
    # apr[residue] = atomsPerResidue(sel[residue])
    # rfa[residue] = residueFirstAtom(sel[residue])
    # coms[residue] = np.zeros((n,n_res[residue],3), dtype=np.float64)
    # md[residue] = np.zeros((n,3),dtype=np.float64)
    print(
        f"Number of {residue}: {n_res[residue]} ranging from {idx_firstlast[residue]}"
    )

sel_im1h = u.select_atoms("resname IM1H and resid 1")
nim1h = sel_im1h.n_residues
mass_im1h = sel_im1h.masses
charge_im1h = sel_im1h.charges
apr_im1h = atomsPerResidue(sel_im1h)
rfa_im1h = residueFirstAtom(sel_im1h)
mdim1h = np.zeros((n, 3), dtype=np.float64)
print("Number IM1H   = ", nim1h)

sel_oac = u.select_atoms("resname OAC and resid 1")
noac = sel_oac.n_residues
mass_oac = sel_oac.masses
charge_oac = sel_oac.charges
apr_oac = atomsPerResidue(sel_oac)
rfa_oac = residueFirstAtom(sel_oac)
mdoac = np.zeros((n, 3), dtype=np.float64)
print("Number OAC   = ", noac)

sel_im1 = u.select_atoms("resname IM1 and resid 1")
nim1 = sel_im1.n_residues
mass_im1 = sel_im1.masses
charge_im1 = sel_im1.charges
apr_im1 = atomsPerResidue(sel_im1)
rfa_im1 = residueFirstAtom(sel_im1)
mdim1 = np.zeros((n, 3), dtype=np.float64)
print("Number IM1   = ", nim1)

sel_hoac = u.select_atoms("resname HOAC and resid 1")
nhoac = sel_hoac.n_residues
mass_hoac = sel_hoac.masses
charge_hoac = sel_hoac.charges
apr_hoac = atomsPerResidue(sel_hoac)
rfa_hoac = residueFirstAtom(sel_hoac)
mdhoac = np.zeros((n, 3), dtype=np.float64)
print("Number HOAC   = ", nhoac)

sel_im = u.select_atoms("resname IM1H or resname IM1")
sel_ac = u.select_atoms("resname OAC or resname HOAC")
print(f"{sel_im.n_atoms=}")
n_im = sel_im.n_residues
n_ac = sel_ac.n_residues


def prepare_dipbyresidue_arguments(
    coor, com, charge, mass, charges, residue_charge, apr
):
    """
    corr: np.array
        coordinates of either im or ac selection
    com: np.array
        center-of-mass coordinates of either im or ac selection
    charge:
        charge of the atoms of a single residue (im1h, oac, im1, hoac)
    mass:
        mass of the atoms of a single residue (im1h, oac, im1, hoac)
    charges:
        charges array with the total charges of either the im or ac selection
    residue_charge: int
        total charge of the residue (im1h: 1, oac: -1, im1,hoac: 0)
    apr: int
        atoms per residue -  number of atoms for either im1h, oac, im1, hoac - attention: im1h/im1 or hoac/oac must have same number of atoms!
    """
    com_short = np.reshape(com[charges == residue_charge], (-1, 3))
    n_res = int(len(com_short))
    charges_list = np.tile(charge, n_res)
    masses_list = np.tile(mass, n_res)
    charges_atoms = np.repeat(charges, apr)
    charges_atoms3D = np.broadcast_to(  # expand the charges to 3D
        charges_atoms, (3, len(charges_atoms))
    ).T  # only view, saves memory, changes in charges_atoms changes charges! for copy use np.tile()
    coor_short = np.reshape(coor[charges_atoms3D == residue_charge], (-1, 3))
    apr_list = np.repeat(apr, n_res)
    rfa_list = np.insert(np.cumsum(apr_list)[:-1], 0, 0).astype(np.int32)
    return coor_short, charges_list, masses_list, n_res, apr_list, rfa_list, com_short


###########################################################################################################################
# Analysis
###########################################################################################################################
md = np.zeros((n, 3), dtype=np.float64)

ctr = 0
start = time.time()
print("")

for ts, step in zip(u.trajectory[::skip], list(charge_steps)[::skip]):
    print(
        "\033[1AFrame %d of %d" % (ts.frame, u.trajectory.n_frames),
        "\tElapsed time: %.2f hours" % ((time.time() - start) / 3600),
    )
    coor_im1h = np.ascontiguousarray(sel_im1h.positions, dtype="double")
    com_im1h = sel_im1h.center_of_mass(compound="residues")

    coor_im = np.ascontiguousarray(
        sel_im.positions, dtype="double"
    )  # array shape needed for C
    coor_ac = np.ascontiguousarray(
        sel_ac.positions, dtype="double"
    )  # array shape needed for C

    com_im = sel_im.center_of_mass(compound="residues")
    com_ac = sel_ac.center_of_mass(compound="residues")

    charges_tmp = np.asarray(charges.data[step])
    charges_im = charges_tmp[
        np.r_[
            idx_firstlast["IM1H"][0] : idx_firstlast["IM1H"][1],
            idx_firstlast["IM1"][0] : idx_firstlast["IM1"][1],
        ]
    ]
    charges_ac = charges_tmp[
        np.r_[
            idx_firstlast["OAC"][0] : idx_firstlast["OAC"][1],
            idx_firstlast["HOAC"][0] : idx_firstlast["HOAC"][1],
        ]
    ]
    assert sum(charges_im) + sum(charges_ac) == 0

    arguments_im1h = prepare_dipbyresidue_arguments(
        coor_im, com_im, charge_im1h, mass_im1h, charges_im, 1, apr_im1h[0]
    )
    arguments_im1 = prepare_dipbyresidue_arguments(
        coor_im, com_im, charge_im1, mass_im1, charges_im, 0, apr_im1[0]
    )
    arguments_oac = prepare_dipbyresidue_arguments(
        coor_ac, com_ac, charge_oac, mass_oac, charges_ac, -1, apr_oac[0]
    )
    arguments_hoac = prepare_dipbyresidue_arguments(
        coor_ac, com_ac, charge_hoac, mass_hoac, charges_ac, 0, apr_hoac[0]
    )

    mdim1h[ctr] += np.sum(dipByResidue(*arguments_im1h), axis=0)
    mdim1[ctr] += np.sum(dipByResidue(*arguments_im1), axis=0)
    mdoac[ctr] += np.sum(dipByResidue(*arguments_oac), axis=0)
    mdhoac[ctr] += np.sum(dipByResidue(*arguments_hoac), axis=0)

    md[ctr] += mdim1h[ctr] + mdim1[ctr] + mdoac[ctr] + mdhoac[ctr]
    ctr += 1

print("calculating <MD(0)MD(t)>")

outfolder = "./mdmd"
Path(outfolder).mkdir(parents=True, exist_ok=True)

# total
md = np.ascontiguousarray(md.T)
md0mdt = np.zeros(n)
correlateParallel(md, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_all.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdim1h = np.ascontiguousarray(mdim1h.T)
md0mdt = np.zeros(n)
correlateParallel(mdim1h, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_im1h.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdim1 = np.ascontiguousarray(mdim1.T)
md0mdt = np.zeros(n)
correlateParallel(mdim1, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_im1.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdoac = np.ascontiguousarray(mdoac.T)
md0mdt = np.zeros(n)
correlateParallel(mdoac, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_oac.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()

mdhoac = np.ascontiguousarray(mdhoac.T)
md0mdt = np.zeros(n)
correlateParallel(mdhoac, md, md0mdt, ltc=1)
f1 = open(f"{outfolder}/md0mdt_hoac.dat", "w")
for i in range(len(md0mdt)):
    f1.write("%5.5f\t%5.5f\n" % (i * skip * dt, md0mdt[i]))
f1.close()
