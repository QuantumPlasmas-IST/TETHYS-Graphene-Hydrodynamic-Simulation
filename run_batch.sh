#!/bin/bash

DIRNAME=$(date +"DataFiles_BATCH_%F_%H:%M")
mkdir -p "$DIRNAME"
mkdir -p "$DIRNAME/extrema"
mkdir -p "$DIRNAME/electronics"

v=10
l=0.0
vis=0.0
for s in 1 2 4 6 8 10 12 14 16 18 20 25 30 40 50 60
do 
	echo "calculating S="$s
	echo "calculating Vf="$v
	echo "calculating Vc="$l
	./RichtmyerHDF5 $s $v $l $vis 0 
	FILENAME1=$(find preview*.dat)
	WORDCOUNT1=$(wc -l preview*.dat)
	LINENUMBER1=${WORDCOUNT1% *}
	./TimeSeries "$LINENUMBER1" "$FILENAME1" "$s" 
	FILENAME2=$(find electro_*.dat)
	WORDCOUNT2=$(wc -l electro_*.dat)
	LINENUMBER2=${WORDCOUNT2% *}
	./ElectronicAnalysis "$LINENUMBER2" "$FILENAME2" "$s" 
	mv -- ElectronicProperties* "./$DIRNAME/electronics"
	mv -- electro_*.dat "./$DIRNAME/electronics"
	mv -- Extrema* "./$DIRNAME/extrema"
	mv -- *.dat "./$DIRNAME"
done 

mv -- *.log "./$DIRNAME"
