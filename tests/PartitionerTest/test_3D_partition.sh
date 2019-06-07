#!/bin/bash

# make a temporary test folder

mkdir tmp

# required source files

depFiles="FluidChannel.py "
depFiles+="partition_compare.py partition_suggestion.py "
depFiles+="pyPartition.py "
depFiles+="vtkHelper.py "
depFiles+="tLBM_partition.py "

# copy necessary source files from the src folder to the tmp folder
cd ../../python

for i in $depFiles; do
  cp $i ../tests/PartitionerTest/tmp/.
done

# move to the temporary test directory
cd ../tests/PartitionerTest/tmp

cp ../test_geom.py .

# make a test geometry
./test_geom.py  #<-- outputs test_geom.mat

# partition the geometry
geom="test_geom.mat"
disc="D3Q27"
style="3D"
parts=8

echo "testing 3D partitioning scripts"

./tLBM_partition.py $geom $disc $style $parts 

exitStat=$?
echo "exitStat = $exitStat"

if [ $exitStat == 0 ]
then
  echo "test successful!"
else
  echo "test failed!"
fi


# clean up 
cd ..
rm -rf tmp

