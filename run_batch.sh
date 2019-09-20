

DIRNAME=$(date +"DataFiles_BATCH_%F_%H:%M")
mkdir -p $DIRNAME


for s in {10,20,40,60}
do 
	echo "calculating S="$s
	./Richtmyer $s 0 | tee -a Richt_full_output.log
	mv *.dat ./$DIRNAME
done 

mv *.log ./$DIRNAME
