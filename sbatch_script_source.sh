#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -p compute
#SBATCH -o /scratch/kmtsui/slurm-%j.out
#SBATCH -e /scratch/kmtsui/slurm-%j.err

#run the application:
source /hepstore/kmtsui/hyperk/setup.sh
cd /bundle/data/T2K/users/kmtsui/LI/fitter
suffix="fullring"

for (( i=1; i <= 10; ++i ))
do
    cat config/config_source_template.toml |\
        sed "s|__LEDERROR__|${i}|g" |\
        sed "s|__NPMT__|${1}|g" > config/config_template_${1}_${i}.toml
    optical_fit -c config/config_template_${1}_${i}.toml -o systematic_study/source_4/${i}pc_mPMT${1}_corrected.root -n 8
    rm config/config_template_${1}_${i}.toml
done
