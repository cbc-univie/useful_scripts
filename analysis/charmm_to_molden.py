import sys

# Write a molden normal mode file from charmm vibran output.
# This script assumes you have printed your coordinates via
# 'coor prin' prior to vibran, and that all normal modes are present.

try:
    charmm_file = sys.argv[1]
except IndexError:
    print("Usage: python charmm_to_molden.py charmm_output")
    sys.exit()

ang2bohr = 1.8897259886

freq = []
X = []
Y = []
Z = []
coords_x = []
coords_y = []
coords_z = []
atoms = []

vibmod   = "  VIBRATION MODE"
eigvec   = "   EIGENVECTOR:"
coordmod = "          COORDINATE FILE MODULE"

try:
    input = open(charmm_file, 'r')
except IOError:
    print("Cannot read file " + charmm_file)
    sys.exit(-1)
x = input.readline()
while len(x) != 0:
    if x[:len(coordmod)] == coordmod:
        l = x.split()
        while l[1] != "EXT":
            x = input.readline()
            l = x.split()
        x = input.readline()
        l = x.split()
        while len(l) != 0:
            atoms.append(l[3][0])
            coords_x.append(float(l[4]))
            coords_y.append(float(l[5]))
            coords_z.append(float(l[6]))
            x = input.readline()
            l = x.split()

    if x[:len(vibmod)] == vibmod:
        #l = string.split(x)
        l = x.split()
        if len(l[3]) > 12:
            if isinstance(l[3][11:], float):
                freq.append(float(l[3][11:]))
            else:
                freq.append(0)
        else:
            freq.append(float(l[4]))
        x = input.readline()
        while x[:len(eigvec)] != eigvec:
            x = input.readline()
        # get next line
        l = input.readline().split()
        while len(l) == 7:
            X.append(float(l[4]))
            Y.append(float(l[5]))
            Z.append(float(l[6]))
            l = input.readline().split()
    x = input.readline()
input.close()

f = open(charmm_file+".molden", "w")
f.write("[Molden Format]\n")
f.write("\n")
f.write("[FREQ]\n")

for i in range(len(freq)-6):
    f.write(f"{freq[i+6]: 23.10f}\n")
f.write("\n")
f.write("[FR-COORD]\n")

for i in range(len(atoms)):
    coord_string=(
    f"{atoms[i]}{coords_x[i]*ang2bohr: 22.10f}"
    f"{coords_y[i]*ang2bohr: 20.10f}{coords_z[i]*ang2bohr: 20.10f}\n"
    )
    f.write(coord_string)

f.write("\n")
f.write("[FR-NORM-COORD]\n")
for i in range(3*len(atoms)-6):
    f.write(f"vibration {i+1}\n") 
    for j in range(len(atoms)):
        mode_string=(
        f"{X[(i+6)*len(atoms)+j]: 23.10f}{Y[(i+6)*len(atoms)+j]: 20.10f}"
        f"{Z[(i+6)*len(atoms)+j]: 20.10f}\n"
        )
        f.write(mode_string)
f.close()