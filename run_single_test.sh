#!/bin/bash

DIRNAME=$(date +"DataFiles_SINGLE_TEST_%F_%H:%M")
mkdir -p "$DIRNAME"
mkdir -p "$DIRNAME/extrema"

s=20
v=10
l=.1
echo "calculating S="$s
echo "calculating Vf="$v
echo "calculating Vc="$l
./Richtmyer $s $v $l 0 | tee -a Fluid.log
#./AnalysisELEC $(wc -l electro*) $s | tee -a Antenna.log
./TimeSeries $(wc -l slice*) $s | tee -a TimeSeries.log

mv -- *.dat "./$DIRNAME"
mv -- *.log "./$DIRNAME"

