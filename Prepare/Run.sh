#!/bin/bash

ROOTScript="Splitter.C"

declare -a FileArray=("single_electron_5-0.8.root" "input_file_2")  # raw SF_e
declare -a RunOverFile=(1 0)
InputDirectory="Input"

for ((i=0;i<${#FileArray[@]};++i))
do
	if [ ${RunOverFile[i]} -eq 1 ]
	then
	#Getting the variables at the current iterator

	File="${FileArray[i]}"

	root -l -b -q $ROOTScript'("'$File'")'
	fi
done

echo ""
echo "All done"
