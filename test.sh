#! /bin/bash
# Run your client for 100 times
# 
# 
# Usage:
#	./test.sh
# Notes:
#	This shell script should be in the same directory as your
#	client implementation and your client should have portnumber 
#	and the number of strings in theArray as command line paramters;
#


# Parameters
Sizes=(5300 13000 22500)
Processes=(1 2 4 8)
Duplicates=50

clear

echo "Specify the Program name and press [Enter]"
read program
echo "Compiling, Running, and Performing initial Testing..." 
mpicc -g -Wall -o $program $program.c Lab4_IO.c -lm
mpirun -n 4 ./$program
#gcc -o $program $program.c Lab4_IO.c -lm
#./$program
./serialtester

for SIZE in ${Sizes[@]}; do
	./datagen -b $SIZE
	for Process in ${Processes[@]}; do
		echo "Starting Size $SIZE for $Process Process"
		ATTEMPT=0
		while [[ $ATTEMPT -ne $Duplicates ]]; do
			let ATTEMPT+=1
			mpirun -n $Process ./$program
			#./$program
			sleep .1
		done
	done
done
