#!/bin/sh

#for (( i=26; i <= ${NUMRUNS}; ++i ))
for (( i=0; i <= 20; ++i ))
do
    timetof=`echo "start=-948;dx=0.5;start+$i*dx" | bc -l` 
    echo ${timetof}
    timetof_string=$(printf "%1.1f" $timetof)

    cat config/config_BnLPMT_template.toml |\
        sed "s|__TIMETOF__|${timetof_string}|g" > config/config_BnLPMT_template_${timetof_string}_run.toml
    echo "[INFO]: Running ${i}-th fitter..."
    optical_fit -c config/config_BnLPMT_template_${timetof_string}_run.toml -o TN/time_cut/diffuser4_400nm_nominal_BnLPMT_timetof_${timetof_string}.root -n 30
    rm config/config_BnLPMT_template_${timetof_string}_run.toml
done
