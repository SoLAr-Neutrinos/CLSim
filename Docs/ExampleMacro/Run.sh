#!/bin/bash

#Define .C script that we want ROOT to interprete later on
ROOTScript="macro.C"

declare -a FileArray=("out")
InputDirectory="."

declare -a RunOverFile=( 1 )

#Defining input and output directory with respect to current directory
#Looping over available  files
for ((i=0;i<${#FileArray[@]};++i))
do
	if [ ${RunOverFile[i]} -eq 1 ]
	then
	#Getting the variables at the current iterator

	File="${FileArray[i]}"
	dEdx="${dEdxArray[i]}"
	Coverage="${Coverage[i]}"


	echo "Currently working on: "$File

	#Define output directories for PNGs and ROOT file
	OutputDirectory="Output"
	PNGOutputFolderName=$OutputDirectory"/PNGs/"$File
	ROOTOutputFolderName=$OutputDirectory"/ROOT"

	#Create output directories (the -p option checks if the output directory already exists and does nothing in case they do)
	mkdir -p $PNGOutputFolderName
	mkdir -p $ROOTOutputFolderName

	echo ": "$File

	root -l -b -q $ROOTScript'("'$File'","'$InputDirectory'", "'$OutputDirectory'")'
	fi
done

echo ""
echo "All done"
