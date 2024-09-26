import os
import time
from pathlib import Path
import re
import argparse
from matplotlib import pyplot as plt
import numpy as np


def read_files(path: str, basename: str, extension: str, firstfile: int, lastfile: int) -> str:
    """ Reads files and returns them as long concatenated string """
    files = [f"{path}/{basename}{i}.{extension}" for i in range(firstfile, lastfile+1)]
    long_string = ""
    for file in files:
        with open(file, "r") as f:
            for line in f.readlines():
                long_string += line
    return long_string

def read_charge_files(path: str, basename: str, extension: str, firstfile: int, lastfile: int) -> str:
    """ Reads files and returns them as long concatenated string """
    files = [f"{path}/{basename}{i}.{extension}" for i in range(firstfile, lastfile+1)]
    long_string = ""
    for file_nr,file in enumerate(files):
        if file_nr >= 1:
            start = 2
        else:
            start = 0
        with open(file, "r") as f:
            for line in f.readlines()[start:]:
                #if start != 0:
                #    number: str = line.split()[0]
                #    line = line.replace(number, str(int(number)+last_number))
                long_string += line
        #last_number = int(long_string.split('\n')[-2].split()[0])
    return long_string

def n_updates(data: str) -> int:
    updates = 0
    for line in data.split("\n"):
        if "len(candidate_pairs)" in line:
            updates += int(line.split("=")[1])
    return updates

def transfers(data: str) -> dict:
    """ Get the proton transfers for each species. 
    Returns dict[frozenset((str,str)), list[tuple[int,int]]] of reaction before it happend
    SO entry in IM1H OAC means that the proton was transfered from IM1H to OAC and got IM1 HOAC
    """
    transfers = {frozenset(("IM1H", "OAC")): [], frozenset(("IM1","HOAC")): [], frozenset(("IM1","IM1H")): [], frozenset(("OAC","HOAC")): []}
    for line in data.split("\n"):
        if "pair accepted" in line:
            #Format:
            #res1:residx1:charge1-res2:residx2:charge2
            # i.e. HOAC:191:0-IM1:124:0
            line = line.split(":")
            line_middle = line[2].rsplit("-",1) # To not split a charge of -1 in the middle term :-1-IM1
            res1, idx1, res2, idx2 = line[0], int(line[1]), line_middle[1], int(line[3])
            transfers[frozenset((res1, res2))].append((idx2, idx1))
    return transfers

def update_trials(data: str) -> int:
    updates = 0
    for line in data.split("\n"):
        if "Update trial" in line:
            updates += 1
    return updates

def get_opposite_res(residue: str):
    """ stupid function to get name after protonation/deprotonation"""
    if residue == "IM1H":
        return "IM1"
    if residue == "IM1":
        return "IM1H"
    if residue == "OAC":
        return "HOAC"
    if residue == "HOAC":
        return "OAC"


def changes(data: str) -> dict:
    "how the evolution with time is"
    transfers = {"IM1H": [150], "OAC": [150], "IM1": [350], "HOAC": [350]}
    for line in data.split("\n"):
        if "pair accepted" in line:
            #Format:
            #res1:residx1:charge1-res2:residx2:charge2
            # i.e. HOAC:191:0-IM1:124:0
            line = line.split(":")
            line_middle = line[2].rsplit("-",1) # To not split a charge of -1 in the middle term :-1-IM1
            res1, idx1, res2, idx2 = line[0], int(line[1]), line_middle[1], int(line[3])
            prev_number1 = transfers[res1][-1]
            prev_number1_opposite = transfers[get_opposite_res(res1)][-1]
            prev_number2 = transfers[res2][-1]
            prev_number2_opposite = transfers[get_opposite_res(res2)][-1]
            transfers[res1].append(prev_number1-1)
            transfers[get_opposite_res(res1)].append(prev_number1_opposite+1)
            transfers[res2].append(prev_number2-1)
            transfers[get_opposite_res(res2)].append(prev_number2_opposite+1)

    return transfers

def current_state(data: str) -> dict:
    "how many of each species right now"
    transfers = {"IM1H": 150, "OAC": 150, "IM1": 350, "HOAC": 350}
    for line in data.split("\n"):
        if "pair accepted" in line:
            #Format:
            #res1:residx1:charge1-res2:residx2:charge2
            # i.e. HOAC:191:0-IM1:124:0
            line = line.split(":")
            line_middle = line[2].rsplit("-",1) # To not split a charge of -1 in the middle term :-1-IM1
            res1, idx1, res2, idx2 = line[0], int(line[1]), line_middle[1], int(line[3])
            transfers[res1] -= 1
            transfers[res2] -= 1
            transfers[get_opposite_res(res1)] += 1
            transfers[get_opposite_res(res2)] += 1

    return transfers

def concentration_timeseries(charge_changes: str) -> list:
    update_steps = []
    im1h, oac, im1, hoac = [],[],[],[]
    for line in charge_changes.split("\n")[2:-2]:
        #print(line)
        line = line.split("\t")
        update_steps.append(int(line[0]))
        anion = line[1].count("-1")
        neutral = line[1].count("0") / 2
        im1h.append(anion)
        oac.append(anion)
        im1.append(neutral)
        hoac.append(neutral)
    return np.asarray(update_steps), [im1h,oac,im1,hoac]






path = "../out/"
basename = "nvt_"
extension = "out"

firstfile=1
lastfile=59
#skip=50

file_str = read_files(path, basename, extension, firstfile, lastfile)
n_updates = n_updates(file_str)
update_trials = update_trials(file_str)
transfers = transfers(file_str)
print(f"{update_trials=}")
print(f"{n_updates=}")
print("IM1H+OAC:", len(transfers[frozenset(("IM1H","OAC"))]))
print("IM1+HOAC:", len(transfers[frozenset(("IM1", "HOAC"))]))
print("IM1H+IM1:", len(transfers[frozenset(("IM1H", "IM1"))]))
print("OAC+HOAC:", len(transfers[frozenset(("OAC", "HOAC"))]))
final_numbers = current_state(file_str)
assert final_numbers["IM1H"] == final_numbers["OAC"]
assert final_numbers["IM1"] == final_numbers["HOAC"]
print("IM1H: ", final_numbers["IM1H"], "OAC: ", final_numbers["OAC"], "IM1: ", final_numbers["IM1"], "HOAC: ", final_numbers["HOAC"])
print("IM1H/OAC : IM1/HOAC",final_numbers["IM1H"]/(final_numbers["IM1H"]+final_numbers["IM1"]), ":", final_numbers["IM1"]/(final_numbers["IM1H"]+final_numbers["IM1"]))
#print(transfers)

charge_changes = read_charge_files(path, "charge_changes_", "out", firstfile, lastfile)
time, n_molecules = concentration_timeseries(charge_changes)
np.savetxt("im1h_series.dat", np.c_[time, np.asarray(n_molecules[0])])
np.savetxt("im1_series.dat", np.c_[time, n_molecules[2]])


quit()
changes = changes(file_str)
plt.plot(changes["IM1H"], label="IM1H")
plt.plot(changes["OAC"], label="OAC")
plt.plot(changes["IM1"], label="IM1")
plt.plot(changes["HOAC"], label="HOAC")
plt.legend()
plt.show()
