#!/bin/bash

DIRNAME=$(date +"TIME_COMPARISON_TEST_%F_%H:%M")
mkdir -p "$DIRNAME"

s=40
v=10
l=.1
echo "calculating S="$s
echo "calculating Vf="$v
echo "calculating Vc="$l

/usr/bin/time -pva --output=RUN_TIME.txt  ./Richtmyer $s $v $l 1 
/usr/bin/time -pva --output=RUN_TIME.txt  ./RichtmyerHDF5 $s $v $l 1 



mv -- *.dat "./$DIRNAME"
mv -- *.h5 "./$DIRNAME"
mv -- RUN_TIME.txt "./$DIRNAME"


