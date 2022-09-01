#!/bin/bash
# ./make_run.sh #fits #filesToDivide
TOTAL=$1
NFILE=$2
FILECOUNTER=0
STARTCOUNTER=0
ENDCOUNTER=0
COUNTERZEROS=0


while [ $FILECOUNTER -lt $NFILE ]; do
FILECOUNTERS=`printf $FILECOUNTER`
touch run_$FILECOUNTERS.sh
(
  let STARTCOUNTER=TOTAL/NFILE*FILECOUNTER
  let ENDCOUNTER=STARTCOUNTER+TOTAL/NFILE
  while [ $STARTCOUNTER -lt $ENDCOUNTER ]; do
      COUNTERZEROS=`printf $STARTCOUNTER`
      printf 'nohup genWeightsThreeTrack.exe -p ../run2w_prod6D_%02d.root -o run2w_prod6D_%02d.root -r MC >out_2w_%02d.txt & \n' $COUNTERZEROS $COUNTERZEROS $COUNTERZEROS 
      let STARTCOUNTER=STARTCOUNTER+1
  done  
)  > run_$FILECOUNTERS.sh
chmod 744 run_$FILECOUNTERS.sh
    let FILECOUNTER=FILECOUNTER+1
done

# for f in ToyMC_scripts/RunToyMC_*.sh; do
#     n=`echo $f | sed -e s/[^0-9]//g`
#     mv $f `printf ToyMC_scripts/RunToyMC_%05d.sh $n`
# done

