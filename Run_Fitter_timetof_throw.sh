#!/bin/sh

for (( i=1; i <= 100; ++i ))
do
    echo "[INFO]: Running ${i}-th fitter..."
    optical_fit -c config/config_mPMT.toml -o TN/time_throw/fitoutput_diffuser4_400nm_nominal_${i}.root -n 30 -s ${i}
done
