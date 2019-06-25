#!/bin/bash

# $1 - CrayFlag



N_divs=15
latticeType=D3Q15
dynamics=2
partitionMethod=1D
nParts=5
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
cp make_empty_channel_2.py tmp/.
cp ../../../tLBM_build/tLBMexec tmp/.

cd tmp

CRAY_FLAG=$1 # [0 if not on cray | 1 if on cray ]

if [ "$CRAY_FLAG" == "1" ]
then

aprun -n 1 ./make_empty_channel_2.py $N_divs

aprun -n 1 ./tLBM_partition.py $MAT_FILE $latticeType $partitionMethod $nParts

aprun -n 1 ./tLBM_write_params.py $MAT_FILE $latticeType $dynamics \
           $partitionMethod $nParts $numTs $tsRepFreq $warmupTs \
           $plotFreq $re $dt $cs $bRestart $bTimeAvg $bSS

export OMP_NUM_THREADS=$nOMP

aprun -n $nParts -d $nOMP ./tLBMexec

else 
mpirun -np 1 ./make_empty_channel_2.py

mpirun -np 1 ./tLBM_partition.py $MAT_FILE $latticeType $partitionMethod $nParts

mpirun -np 1 ./tLBM_write_params.py $MAT_FILE $latticeType $dynamics \
           $partitionMethod $nParts $numTs $tsRepFreq $warmupTs \
           $plotFreq $re $dt $cs $bRestart $bTimeAvg $bSS

export OMP_NUM_THREADS=$nOMP

#mpirun -np -n $nParts ./tLBMexec

fi

# test complete; clean up

#cd ..
#rm -rf tmp
