#!/bin/bash

source /home/marion/miniconda3/bin/activate openmm

inp_file="npt.py"
last=10
path="/site/raid6/marion/il/dil2"

for run in run1 run2 run3 run4
do
	cd ${run}
	sbatch -J ${run}_$1 -o out/npt_$1.out run_npt.sh $1
	cd ../
done

for run in run1 run2 run3 run4
do
	while [ ! -f ${path}/${run}/traj/npt_$1.rst ] 
	do
        	sleep 0.1
	done
done

python3 determine_next_temp.py $1

ret_val=$?

if [ $ret_val -eq 0 -a $1 -lt $last ] ; then
	let next=$1+1
	sh startall.sh $next
elif [ ! $ret_val -eq 0 ] ; then
	echo "Error in run $1 .." >> out/error.log
fi
