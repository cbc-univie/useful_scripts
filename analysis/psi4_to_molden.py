import sys
import string

# Write a molden normal mode file from psi4 output.
# This script assumes you have printed your normal modes
# once at the end of your file.
# NOTE: Psi4 can do this natively, this script is for the forgetful

try:
    psi4_file = sys.argv[1]
except IndexError:
    print("Usage: python psi4_to_molden.py psi4_output")
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

freq_str   = "  Freq [cm^-1]"
coordmod   = "    Final optimized geometry"

try:
    input = open(psi4_file, 'r')
except IOError:
    print("Cannot read file " + psi4_file)
    sys.exit(-1)
x = input.readline()
while len(x) != 0:
    if x[:len(coordmod)] == coordmod:
        empty_lines = 0
        while empty_lines < 2:
            x = input.readline()
            if x == '\n':
                empty_lines += 1
        x = input.readline()
        l = x.split()
        while len(l) != 0:
            atoms.append(l[0])
            coords_x.append(float(l[1]))
            coords_y.append(float(l[2]))
            coords_z.append(float(l[3]))
            x = input.readline()
            l = x.split()

    if x[:len(freq_str)] == freq_str:
        l = x.split()
        for i in range(3):
            freq.append(float(l[i + 2]))
        l = input.readline().split()
        while l[0][0] not in string.digits:
            l = input.readline().split()
        bx1 = []
        by1 = []
        bz1 = []
        bx2 = []
        by2 = []
        bz2 = []
        bx3 = []
        by3 = []
        bz3 = []
        while len(l) == 11:
            bx1.append(float(l[2]))
            by1.append(float(l[3]))
            bz1.append(float(l[4]))
            bx2.append(float(l[5]))
            by2.append(float(l[6]))
            bz2.append(float(l[7]))
            bx3.append(float(l[8]))
            by3.append(float(l[9]))
            bz3.append(float(l[10]))
            l = input.readline().split()
        X = X + bx1 + bx2 + bx3
        Y = Y + by1 + by2 + by3
        Z = Z + bz1 + bz2 + bz3
    x = input.readline()

input.close()

f = open(psi4_file+".molden", "w")
f.write("[Molden Format]\n")
f.write("\n")
f.write("[FREQ]\n")

for i in range(len(freq)):
    f.write(f"{freq[i]: 23.10f}\n")
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
        f"{X[i*len(atoms)+j]: 23.10f}{Y[i*len(atoms)+j]: 20.10f}"
        f"{Z[i*len(atoms)+j]: 20.10f}\n"
        )
        f.write(mode_string)
f.close()