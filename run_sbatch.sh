#!/bin/sh

for (( i=1; i <= 40; ++i ))
do
    mPMT=`echo "start=0;dx=100;start+$i*dx" | bc -l` 
    #echo "sbatch -J fit_${i} sbatch_script.sh ${mPMT}"
    #sbatch -J fit_${i} sbatch_script.sh ${mPMT}
    #sbatch -J fit_${i} sbatch_script_eff.sh ${mPMT}
    sbatch -J fit_${i} sbatch_script_source.sh ${mPMT}
    #sbatch -J fit_${i} sbatch_script.sh ${i}
done