#!/bin/bash

DIR_ROOT=$(date +"DataFiles_CFL_%F_%H:%M")
DIR_IMG="$DIR_ROOT/Images"

mkdir -p "$DIR_ROOT"
mkdir -p "$DIR_IMG"

rm *.dat
rm *.log

l=0

V_MIN=3
V_RANGE=50
V_ITER=1
S_ITER=2
S_RANGE=4

v_iterator () {
	V_ITER=1
	if (( $1 >= 8 )); then
		V_ITER=2
	fi
	if (( $1 >= 30 )); then
		V_ITER=5
	fi
}

s_iterator () {
	S_ITER=2
	if (( $1 >= 30 )); then
		S_ITER=5
	fi
}

s_range () {
	S_RANGE=4
	if (( $1 >= 8 )); then
		S_RANGE=2
	fi
	if (( $1 >= 15 )); then
		S_RANGE=4
	fi
	if (( $1 >= 30 )); then
		S_RANGE=5
	fi
}


for((v=V_MIN; v<=V_RANGE; v+=V_ITER))
do
	v_iterator $v
	s_iterator $v
	s_range $v
        DIR_SUB="./$DIR_ROOT/Vf=$v"
	mkdir -p "$DIR_SUB"
	echo "calculating Vf="$v
	for ((s = v; s<= v + S_RANGE; s+=S_ITER))
		do 
		echo "	calculating S="$s
		IMG_DEN_NAME="S=$s.Vf=$v.l=$l.DEN.png"
		IMG_VEL_NAME="S=$s.Vf=$v.l=$l.VEL.png"
		
		echo "	Running ......"; ./Richtmyer $s $v $l 0  &> /dev/null

		FILE_NAME="slice_S=$s.00vF=$v.00l=$l.00.dat"

		echo > gnuplot.in
		echo "set term png
		set o '$IMG_DEN_NAME'
		set title 'den-t Vf=$v S=$s'
		plot '$FILE_NAME' u 1:2 w l noti
		set o '$IMG_VEL_NAME'
		set title 'vel-t Vf=$v S=$s'
		plot '$FILE_NAME' u 1:3 w l noti
		set o; set term pop" >> gnuplot.in
		gnuplot gnuplot.in
		mv -- *.dat "$DIR_SUB"
		mv -- *.png "$DIR_IMG"
	done 
done
	
mv -- *.log "$DIR_ROOT"
