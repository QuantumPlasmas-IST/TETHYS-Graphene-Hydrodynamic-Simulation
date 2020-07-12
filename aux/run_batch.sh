#!/bin/bash

DIRNAME=$(date +"DataFiles_VISCOSITY_TRIAL_BATCH_FTCS_with_Goduov_splitting_%F_%H:%M")
mkdir -p "$DIRNAME"
#mkdir -p "$DIRNAME/extrema"
#mkdir -p "$DIRNAME/electronics"

v=10
l=0.0
vis=0.0
#for s in 1 2 4 6 8 10 12 14 16 18 20 25 30 40 50 60
for s in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80
do 
	for vis in 0.0 0.02 0.05 0.1 0.12 0.15
	do 
		echo "calculating S="$s
		echo "calculating Vf="$v
		echo "calculating Vc="$l
		./TETHYS_2D $s $v $l $vis 1 
#		FILENAME1=$(find preview*.dat)
#		WORDCOUNT1=$(wc -l preview*.dat)
#		LINENUMBER1=${WORDCOUNT1% *}
#		./TimeSeries "$LINENUMBER1" "$FILENAME1" "$s" 	
#		FILENAME2=$(find electro_*.dat)
#		WORDCOUNT2=$(wc -l electro_*.dat)
#		LINENUMBER2=${WORDCOUNT2% *}
#		./ElectronicAnalysis "$LINENUMBER2" "$FILENAME2" "$s" 
#		mv -- electro_S=*vF=*vis=*l=*.dat "./$DIRNAME/electronics"
#		mv -- *electro_*.dat "./$DIRNAME/electronics"
#		mv -- extrema_* "./$DIRNAME/extrema"
		mv -- *.dat "./$DIRNAME"
		mv -- *.h5 "./$DIRNAME"
	done
done 

mv -- *.log "./$DIRNAME"
