

DIRNAME=$(date +"DataFiles_SINGLE_TEST_%F_%H:%M")
mkdir -p $DIRNAME

s=20

echo "calculating S="$s
./Richtmyer $s 0 | tee -a Fluid.log

mv *.dat ./$DIRNAME
mv *.log ./$DIRNAME
