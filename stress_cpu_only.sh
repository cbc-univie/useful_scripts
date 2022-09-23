#!/bin/bash
#SBATCH -p lcpu

/home/andras/stress/usr/bin/stress --cpu `nproc` --vm `nproc` --vm-bytes 1GB --io `nproc` --timeout 86400            
