#!/bin/bash

DIRNAME=$(date +"DataFiles_CFL_TRIAL_BATCH_%F_%H:%M")
mkdir -p "$DIRNAME"

v=10.0
l=0.0
vis=0.0
#for s in 10 20 30 40 50 60 70 80
for s in 60 70 80
do 
	echo "calculating S="$s
	echo "calculating Vf="$v
	echo "calculating Vc="$l
	./TETHYS_2D $s $v $l $vis 0 
	mv -- *.dat "./$DIRNAME"
done 

mv -- *.log "./$DIRNAME"
