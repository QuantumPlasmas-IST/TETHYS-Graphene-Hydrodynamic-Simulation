#!/bin/bash

DIRNAME=$(date +"DataFiles_SINGLE_TEST_%F_%H:%M")
mkdir -p "$DIRNAME"

s=20
v=10
l=10
echo "calculating S="$s
echo "calculating Vf="$v
echo "calculating Vc="$l
./Richtmyer $s $v $l 0 | tee -a Fluid.log

mv -- *.dat "./$DIRNAME"
mv -- *.log "./$DIRNAME"
