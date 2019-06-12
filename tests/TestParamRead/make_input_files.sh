#!/bin/bash

N_divs=11
latticeType=D3Q15
dynamics=2
partitionMethod=metis
nParts=4
nOMP=1
bPreProcess=1
bRestart=0
bTimeAvg=1

MAT_FILE=empty_chan.mat

numTs=1001
tsRepFreq=500
warmupTs=0
plotFreq=500
re=125
dt=0.004
cs=5
bSS=0

#create a test directory and copy over necessary files.

mkdir tmp
cp ../../python/*.py tmp/.
cp make_empty_channel.py tmp/.
cp ../../build/tLBM tmp/.

cd tmp


./make_empty_channel.py

./tLBM_partition.py $MAT_FILE $latticeType $partitionMethod $nParts

./tLBM_write_params.py $MAT_FILE $latticeType $dynamics \
           $partitionMethod $nParts $numTs $tsRepFreq $warmupTs \
           $plotFreq $re $dt $cs $bRestart $bTimeAvg $bSS


