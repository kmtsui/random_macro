#!/bin/sh

for (( i=0; i < 100; ++i ))
do
    echo "[INFO]: Running ${i}-th collimator file..."
    filename=$(dirac-dms-lfn-replicas /hyperk.org/beta-production/hk-calibration-with-mpmt/wcsim/nickwp-2.0.3-22-g9178f1d/nominal/LI/collimator/400nm/32/out/collimator32_400nm_nominal_${i}.root | awk '/root:\/\// {print $3}')
    WCSIM_TreeConvert -f ${filename} -o out_collimator32_nominal_${i}.root -l 400 -d
done
