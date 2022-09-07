#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu

cd /home/andras/gpu-burn
/home/andras/stress/usr/bin/stress --cpu `nproc` --vm `nproc` --vm-bytes 1GB --io `nproc` --timeout 86400 &
./gpu_burn 86400
