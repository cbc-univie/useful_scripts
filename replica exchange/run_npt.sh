#!/bin/bash

#SBATCH -p gpu
#SBATCH --gres=gpu

source /home/marta/miniconda3/bin/activate openmm7.6

inp_file="npt.py"

python3 ${inp_file} $1
