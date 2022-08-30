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
cat config/config_template.toml |\
        sed "s|__PMTMASK__|${1}|g" > config/config_template_${1}_${suffix}.toml

# for (( i=1; i <= 100; ++i ))
# do
#     optical_fit -c config/config_template_${1}_${suffix}.toml -o joint_eff_study/diffuser4_400nm_nominal_mPMT_${1}_${suffix}_${i}.root -s ${i} -n 8
# done

optical_fit -c config/config_template_${1}_${suffix}.toml -o systematic_study/source_3/diffuser4_400nm_nominal_mPMT_${1}_${suffix}.root -n 8 -t 100

#optical_fit -c config/config_template_${1}_${suffix}.toml -o mPMT_number_study/diffuser4_400nm_nominal_mPMT_${1}_${suffix}.root -n 8
rm config/config_template_${1}_${suffix}.toml