#! /bin/bash
# source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

STATUS=write
EXE=PileupHandler2.exe

inputFile=/home/chencheng/data/Run4U/oldfitter/hists.root
# inputFile=/home/chencheng/data/Run4U/oldfitter/hists_PU_All.root
outputFile=/home/chencheng/data/Run4U/oldfitter/hists_PU_triple_func_All_test.root
# mkdir -p ${outputDir}
g++ `root-config --cflags --glibs` -o ${EXE}.exe PileupHandler.cpp &&\
./${EXE}.exe -i ${inputFile} -o ${outputFile} -s ${STATUS}
echo -e "Finished"