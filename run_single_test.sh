#!/bin/bash

DIRNAME=$(date +"DataFiles_SINGLE_TEST_%F_%H:%M")
mkdir -p "$DIRNAME"

s=20
v=10
echo "calculating S="$s
echo "calculating Vf="$v
./Richtmyer $s $v 0 | tee -a Fluid.log

mv -- *.dat ./$DIRNAME
mv -- *.log ./$DIRNAME
