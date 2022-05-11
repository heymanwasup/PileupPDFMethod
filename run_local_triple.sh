#! /bin/bash
g++ `root-config --cflags --glibs` -o PileupConstructer.exe PileupConstructer.cpp \
&& ./run_condor.sh Newfitter All 1