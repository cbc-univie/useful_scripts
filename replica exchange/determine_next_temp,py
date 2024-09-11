import numpy as np
import re
from itertools import permutations, product
import time
import sys
import os
from pathlib import Path

#paths and stats
path = "/site/raid6/marion/il/dil2"
max_tempdiff = 200 #maximal temperature difference between two runs that can be switched same unit as in npt script

#list the directories which should be changed, here they are called run1-4
runs = ["run1","run2","run3","run4"]

log = "./templog"
Path(log).mkdir(parents=True, exist_ok=True)

#counter
cnt = sys.argv[1]

#open files as strings
inp_f = {}

for run in runs:
    with open(f"{path}/{run}/npt.py","r") as f:
        inp_f[run] = f.read()

#create dict to match runs to temperatures
temperatures = {}

for run in runs:
    match = re.search("temperature = ", inp_f[run])
    if match:
        temperatures[run] = float(match.string[match.end():][0:6])

#create possible switching pairs
combs = list(product(temperatures, repeat=2))

possible_switches=[]

for comb in combs:
    if np.abs(temperatures[comb[0]] - temperatures[comb[1]]) <= max_tempdiff:
        possible_switches.append(comb)

#determine which are switched

switches = []

while len(possible_switches) != 0:
    switch = np.random.randint(0,len(possible_switches))
    switch_tuple = possible_switches[switch]
    switches.append(switch_tuple)
    with open(f"{path}/{switch_tuple[0]}/npt.py", "r") as f1:
        lines = f1.readlines()
    with open(f"{path}/{switch_tuple[0]}/npt.py", "w") as f1:
        for line in lines:
            f1.write(re.sub(f'{temperatures[switch_tuple[0]]}', f'{temperatures[switch_tuple[1]]}', line))
    with open(f"{path}/{switch_tuple[1]}/npt.py", "r") as f2:
        lines = f2.readlines()
    with open(f"{path}/{switch_tuple[1]}/npt.py", "w") as f2:
        for line in lines:
            f2.write(re.sub(f'{temperatures[switch_tuple[1]]}', f'{temperatures[switch_tuple[0]]}', line))
    possible_switches = [sw for sw in possible_switches if sw[0] not in switch_tuple and sw[1] not in switch_tuple]

#create templog
with open(f"{path}/templog/templog_{cnt}.dat", "w") as f3:
    f3.write(f"temperature has been changed: {switches}")
f3.close
