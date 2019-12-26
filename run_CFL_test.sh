#!/bin/bash

DIRNAME1=$(date +"DataFiles_CFL_%F_%H:%M")
DIRNAMEIMG=$(date +"Image_Vf=$v")
mkdir -p "$DIRNAME1"
mkdir "./$DIRNAME1/$DIRNAMEIMG"

l=0
for((v=3; v<=8; v+=1))
do
	DIRNAME2=$(date +"DataFiles_Batch_Vf=$v")
	mkdir "./$DIRNAME1/$DIRNAME2"
	for ((s = v; s<= v + 4; s+=2))
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
		echo "set term png; set output \"density_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:2 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "set ylabel \"velocity\"" >> gnuplot.in
		echo "set title \"vel-t Vf=$v S=$s\"" >> gnuplot.in
		echo "set output \"velocity_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:3 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "quit" >> gnuplot.in
		gnuplot gnuplot.in
		mv -- *.dat "./$DIRNAME1/$DIRNAME2/$DIRNAME3"
		mv -- *.png "./$DIRNAME1/$DIRNAMEIMG"
		done 
	done
	
for((v=9; v<=15; v+=2))
do
	DIRNAME2=$(date +"DataFiles_Batch_Vf=$v")
	mkdir "./$DIRNAME1/$DIRNAME2"
	for ((s = v; s<= v + 2 ; s+=2))
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
		echo "set term png; set output \"density_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:2 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "set ylabel \"velocity\"" >> gnuplot.in
		echo "set title \"vel-t Vf=$v S=$s\"" >> gnuplot.in
		echo "set output \"velocity_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:3 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "quit" >> gnuplot.in
		gnuplot gnuplot.in
		mv -- *.dat "./$DIRNAME1/$DIRNAME2/$DIRNAME3"
		mv -- *.png "./$DIRNAME1/$DIRNAMEIMG"
		done 
	done

for((v=15; v<=30; v+=2))
do
	DIRNAME2=$(date +"DataFiles_Batch_Vf=$v")
	mkdir "./$DIRNAME1/$DIRNAME2"
	for ((s = v; s<= v 2 4; s+=2))
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
		echo "set term png; set output \"density_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:2 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "set ylabel \"velocity\"" >> gnuplot.in
		echo "set title \"vel-t Vf=$v S=$s\"" >> gnuplot.in
		echo "set output \"velocity_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:3 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "quit" >> gnuplot.in
		gnuplot gnuplot.in
		mv -- *.dat "./$DIRNAME1/$DIRNAME2/$DIRNAME3"
		mv -- *.png "./$DIRNAME1/$DIRNAMEIMG"
		done 
	done


for((v=30; v<=50; v+=5))
do
	DIRNAME2=$(date +"DataFiles_Batch_Vf=$v")
	mkdir "./$DIRNAME1/$DIRNAME2"
	for ((s = v; s<= v + 5; s+=5))
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
		echo "set term png; set output \"density_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:2 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "set ylabel \"velocity\"" >> gnuplot.in
		echo "set title \"vel-t Vf=$v S=$s\"" >> gnuplot.in
		echo "set output \"velocity_$s.00Vf=$v.00l=$l.00.png\"; plot \"slice_S=$s.00Vf=$v.00l=$l.00.dat\" u 1:3 with linespoints ps 0 lw 0.1" >> gnuplot.in
		echo "quit" >> gnuplot.in
		gnuplot gnuplot.in
		mv -- *.dat "./$DIRNAME1/$DIRNAME2/$DIRNAME3"
		mv -- *.png "./$DIRNAME1/$DIRNAMEIMG"
		done 
	done

mv -- *.log "./$DIRNAME1"