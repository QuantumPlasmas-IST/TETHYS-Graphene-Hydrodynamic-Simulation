#!/bin/bash

DIRNAME=$(date +"DataFiles_BATCH_%F_%H:%M")
mkdir -p "$DIRNAME"

v=10
l=10
for s in 10 20 40 60
do 
	echo "calculating S="$s
	echo "calculating Vf="$v
	echo "calculating Vc="$l
	./Richtmyer $s $v $l 0 | tee -a Richt_full_output.log
	mv -- *.dat "./$DIRNAME"
done 

mv -- *.log "./$DIRNAME"
