for (( i=1; i <= 1; ++i ))
do
    #timetof=`echo "start=-945;dx=0.1;start+$i*dx" | bc -l` 
    #nmpmts=`echo "start=0;dx=100;start+$i*dx" | bc -l` 
    intensityErr=`echo "start=0;dx=100;start+$i*dx" | bc -l` 
    echo ${intensityErr}
    #timetof_string=$(printf "%1.1f" $timetof)
    dirname=$(printf "TN_update/ensemble/noRing2/%imPMTs/" $intensityErr)
    echo ${dirname}
    mkdir -p ${dirname}

    cat config/config_template.toml |\
        sed "s|__NMPMTS__|${intensityErr}|g" > config/config_mPMT_template_${intensityErr}_run.toml
    echo "[INFO]: Running ${i}-th fitter..."
    optical_fit -c config/config_mPMT_template_${intensityErr}_run.toml -o ${dirname}/diffuser4_400nm.root -n 8 -t 100
    rm config/config_mPMT_template_${intensityErr}_run.toml
done
