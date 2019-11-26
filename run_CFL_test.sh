#!/bin/bash

DIRNAME1=$(date +"DataFiles_CFL_%F_%H:%M")
mkdir -p "$DIRNAME1"

l=0
for v in 5 10 15 15 20 25 30
do
	DIRNAME2=$(date +"DataFiles_Batch_Vf=$v")
	mkdir "./$DIRNAME1/$DIRNAME2"
	for ((s = $v; s<=30; s+=5))
		do 
		DIRNAME3=$(date +"DataFiles_Single_Vf=$v S=$s")
		mkdir "./$DIRNAME1/$DIRNAME2/$DIRNAME3"
		echo "calculating S="$s
		echo "calculating Vf="$v
		echo "calculating Vc="$l
		./Richtmyer $s $v $l 0 | tee -a Richt_full_output.log
		echo > gnuplot.in
		echo "set xlabel \"time\"" >> gnuplot.in
		echo "set ylabel \"density\"" >> gnuplot.in
		echo "set title \"den-t Vf=$v S=$s\"" >> gnuplot.in
		echo "set term png; set output \"density.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:2 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "set ylabel \"velocity\"" >> gnuplot.in
		echo "set title \"vel-t Vf=$v S=$s\"" >> gnuplot.in
		echo "set output \"velocity.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:3 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "quit" >> gnuplot.in
		gnuplot gnuplot.in
		mv -- *.dat "./$DIRNAME1/$DIRNAME2/$DIRNAME3"
		mv -- *.png "./$DIRNAME1/$DIRNAME2/$DIRNAME3"
	done 
done

mv -- *.log "./$DIRNAME1"