#!/bin/bash

# $1 - CrayFlag



N_divs=5
latticeType=D3Q15
dynamics=1
partitionMethod=1D
nParts=1
nOMP=1
bPreProcess=1
bRestart=0
bTimeAvg=0

MAT_FILE=wall_mounted_brick.mat

numTs=1001
tsRepFreq=100
warmupTs=0
plotFreq=100
re=10
dt=0.02
cs=1
bSS=0

CRAY_FLAG=$1 # [0 if not on cray | 1 if on cray ]

if [ "$CRAY_FLAG" == "1" ]
then

aprun -n 1 ./wmb_geom.py $N_divs

aprun -n 1 ./tLBM_partition.py $MAT_FILE $latticeType $partitionMethod $nParts

aprun -n 1 ./tLBM_write_params.py $MAT_FILE $latticeType $dynamics \
           $partitionMethod $nParts $numTs $tsRepFreq $warmupTs \
           $plotFreq $re $dt $cs $bRestart $bTimeAvg $bSS

export OMP_NUM_THREADS=$nOMP

aprun -n $nParts -d $nOMP ./tLBMexec

aprun -n 1 ./processTLBM.py

else 
mpirun -np 1 ./wmb_geom.py $N_divs

mpirun -np 1 ./tLBM_partition.py $MAT_FILE $latticeType $partitionMethod $nParts

mpirun -np 1 ./tLBM_write_params.py $MAT_FILE $latticeType $dynamics \
           $partitionMethod $nParts $numTs $tsRepFreq $warmupTs \
           $plotFreq $re $dt $cs $bRestart $bTimeAvg $bSS

export OMP_NUM_THREADS=$nOMP

mpirun -np $nParts ./tLBMexec

mpirun -np 1 ./processTLBM.py

fi


