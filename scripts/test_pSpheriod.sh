#!/bin/bash

# $1 - CrayFlag
# $2 - nParts


N_divs=33
latticeType=D3Q27
dynamics=1
partitionMethod=metis
nParts=$2
nOMP=1
bPreProcess=1
bRestart=0
bTimeAvg=0

MAT_FILE=prolate_spheroid.mat

numTs=50001
tsRepFreq=1000
warmupTs=0
plotFreq=10000
re=100
dt=0.003
cs=15
bSS=0

CRAY_FLAG=$1 # [0 if not on cray | 1 if on cray ]

if [ "$CRAY_FLAG" == "1" ]
then

if [ "$bPreProcess" == "1" ]

aprun -n 1 ./psph_geom.py $N_divs

aprun -n 1 ./tLBM_partition.py $MAT_FILE $latticeType $partitionMethod $nParts

aprun -n 1 ./tLBM_write_params.py $MAT_FILE $latticeType $dynamics \
           $partitionMethod $nParts $numTs $tsRepFreq $warmupTs \
           $plotFreq $re $dt $cs $bRestart $bTimeAvg $bSS
fi

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


