#!/bin/bash

DIRNAME=$(date +"DataFiles_CFL_TRIAL_BATCH_%F_%H:%M")
mkdir -p "$DIRNAME"
mkdir -p "$DIRNAME/electronics"

v=5
l=0.0
vis=0.0
#for s in 1 2 4 6 8 10 12 14 16 18 20 25 30 40 50 60
for s in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80
do 
	echo "calculating S="$s
	echo "calculating Vf="$v
	echo "calculating Vc="$l
	./TETHYS_1D $s $v $l $vis 0 
	mv -- *electro_*.dat "./$DIRNAME/electronics"
	mv -- *.dat "./$DIRNAME"
done 

mv -- *.log "./$DIRNAME"
