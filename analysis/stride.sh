#!/bin/bash
#SBATCH -p lgpu

module load parallel
parallel /home/andras/bin/vmd -dispdev text -eofexit -e {} ::: stride_cluster*.tcl
