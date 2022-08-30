#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -p short
#SBATCH -o /scratch/kmtsui/slurm-%j.out
#SBATCH -e /scratch/kmtsui/slurm-%j.err

#run the application:
source /hepstore/kmtsui/hyperk/setup.sh
cd /bundle/data/T2K/users/kmtsui/LI/fitter
suffix="noring2_pol"
cat config/config_template.toml |\
        sed "s|__PMTMASK__|${1}|g" > config/config_template_${1}_${suffix}.toml
optical_fit -c config/config_template_${1}_${suffix}.toml -o mPMT_number_study/diffuser4_400nm_nominal_mPMT_${1}_${suffix}.root -n 8
rm config/config_template_${1}_${suffix}.toml