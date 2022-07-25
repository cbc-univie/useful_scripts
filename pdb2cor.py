import parmed
import sys

infile = sys.argv[1]
basename = infile.rsplit(".")[0]

coords = parmed.load_file(infile)
coords.save(f"{basename}.crd")
