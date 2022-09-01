#!/bin/sh

NUMRUNS=1

densitySwitch_Static=2

MA=0.95
b_proton_pinull=3.0
b_neutron_piplus=1.5

which_resonanceModel=0

dirname="ND280_GiBUU_CH_numuCC"

numTimeSteps=150

mediumSwitch_Delta=".true."

singlePiModel=1

EQS_Type=5

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      -name)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      dirname="$2"
      echo "[INFO]: Making directory ${2}"
      shift # past argument
      ;;

      -n|--num-cc-jobs)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NUMRUNS=${2}
      echo "[INFO]: Running ${2}x stats..."
      shift # past argument
      ;;

      -e)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      EQS_Type=${2}
      echo "[INFO]: EQS_Type = ${2}"
      shift # past argument
      ;;

      -d|--density-switch)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      densitySwitch_Static=${2}
      echo "[INFO]: densitySwitch_Static = ${2}"
      shift # past argument
      ;;

      -r|--resoance-model)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      which_resonanceModel=${2}
      echo "[INFO]: which_resonanceModel = ${2}"
      shift # past argument
      ;;

      -b|--bnl-switch)

      MA=1.3
      b_proton_pinull=6.0
      b_neutron_piplus=3.0
      echo "[INFO]: Using BNL fit values"
      ;;

      -nofsi)

      numTimeSteps=0
      echo "[INFO]: No FSI in C"
      ;;

      -B)

      mediumSwitch_Delta=".false."
      echo "[INFO]: Do not use collisional broadening of baryonic resonances"
      ;;

      -N)

      singlePiModel=0
      echo "[INFO]: Using pi production model according to Nieves et al (hep-ph/0701149) "
      ;;

      *)
              # unknown option
      echo "Unknown option $1"
      exit 1
      ;;
  esac
  shift # past argument or value
done

CARDDIR=${GIBUUTOOLSROOT}/../../jobcards

#if ! hash MakeBUUInputSoftlink &> /dev/null; then
#  echo "[ERROR]: Cannot find the command \"GiBUU.x\", have you built GiBUU?"
#  exit 1
#fi

#if ! hash GiBUU.x &> /dev/null; then
#  echo "[ERROR]: Cannot find the command \"GiBUU.x\", have you built GiBUU?"
#  exit 1
#fi

#if ! hash GiBUUToStdHep &> /dev/null; then
#  echo "[ERROR]: Cannot find the command \"GiBUUToStdHep\", have you built the GiBUU tools?"
#  exit 1
#fi

OLDPWD=$(pwd)

mkdir -p ${dirname}

cd ${dirname}

#MakeBUUInputSoftlink
ln -s /hepstore/kmtsui/T2K/GiBUU/buuinput BUUInput

cat ${CARDDIR}/ND280_C_numuCC.job |\
 sed "s|__BUUINPUT__|../BUUInput|g" |\
 sed "s|__NUMRUNS__|${NUMRUNS}|g" |\
 sed "s|__DWS__|${densitySwitch_Static}|g" |\
 sed "s|__MA__|${MA}|g" |\
 sed "s|__BPP__|${b_proton_pinull}|g" |\
 sed "s|__BNP__|${b_neutron_piplus}|g" |\
 sed "s|__WRS__|${which_resonanceModel}|g" |\
 sed "s|__NTS__|${numTimeSteps}|g" |\
 sed "s|__MSD__|${mediumSwitch_Delta}|g" |\
 sed "s|__SPM__|${singlePiModel}|g" |\
 sed "s|__EQS__|${EQS_Type}|g" |\
 sed "s|__SEED__|${RANDOM}|g" > ND280_C_numuCC.job

cat ${CARDDIR}/ND280_H_numuCC.job |\
 sed "s|__BUUINPUT__|../BUUInput|g" |\
 sed "s|__NUMRUNS__|${NUMRUNS}|g" |\
 sed "s|__DWS__|${densitySwitch_Static}|g" |\
 sed "s|__MA__|${MA}|g" |\
 sed "s|__BPP__|${b_proton_pinull}|g" |\
 sed "s|__BNP__|${b_neutron_piplus}|g" |\
 sed "s|__WRS__|${which_resonanceModel}|g" |\
 sed "s|__SPM__|${singlePiModel}|g" |\
 sed "s|__SEED__|${RANDOM}|g" > ND280_H_numuCC.job

mkdir -p C_numuCC
cd C_numuCC
echo "[INFO]: Running C events..."
if ! GiBUU.x < ../ND280_C_numuCC.job  > gibuu.C.run; then
  echo "[ERROR]: Failed to run C events:"
  tail -20 gibuu.C.run
  cd ${OLDPWD}
  exit 1
fi
cd ..

mkdir -p H_numuCC
cd H_numuCC

echo "[INFO]: Running H events..."
if ! GiBUU.x < ../ND280_H_numuCC.job  > gibuu.H.run; then
  echo "[ERROR]: Failed to run H events:"
  tail -20 gibuu.H.run
  cd ${OLDPWD}
  exit 1
fi
cd ..

echo "[INFO]: Translating to stdhep..."
if ! GiBUUToStdHep -NP -u 14 -z 6 -a 12 -W 12 -f C_numuCC/FinalEvents.dat \
  -z 1 -a 1 -W 1 -f H_numuCC/FinalEvents.dat \
  -F numu_flux,BUUInput/neutrino/T2K_ND280_250kA-numu.dat -R i13 \
  -o ND280_CH_numuCC.stdhep.root; then
  echo "[ERROR]: Failed to translate to stdhep format."
  cd ${OLDPWD}
  exit 1
fi

echo "[INFO]: Successfully GiBUU stdhep events."

cd ${OLDPWD}
