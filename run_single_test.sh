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

FILENAME=$(find slice*.dat)
WORDCOUNT=$(wc -l slice*.dat)
LINENUMBER=${WORDCOUNT% *}

./TimeSeries $LINENUMBER "$FILENAME" $s | tee -a TimeSeries.log

mv -- Extrema* "./$DIRNAME/extrema"
mv -- *.dat "./$DIRNAME"
mv -- *.log "./$DIRNAME"

