#! /bin/bash
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
EXE=PileupConstructer.exe

fitter=${1}
threshold=${2}
caloNum=$((${3}%24+1))

echo -e "submitting on ${threshold}, calorimeter ${caloNum}"

InputFile=/home/chencheng/data/Run4U/PUTest/${fitter}_NewBinning_Cuts${threshold}MeV.root
OutputFile=/home/chencheng/data/Run4U/PUTest/${fitter}_NewBinning_Cuts${threshold}MeV_PUCalos/hists_PU_calo${caloNum}.root

mkdir -p `dirname ${OutputFile}`

./${EXE} -i ${InputFile} -o ${OutputFile} -n ${caloNum} -e ${threshold}
echo -e "Finished on ${threshold} calorimeter ${caloNum}"