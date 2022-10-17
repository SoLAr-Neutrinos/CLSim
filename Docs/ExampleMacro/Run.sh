#!/bin/bash

#Define .C script that we want ROOT to interprete later on
ROOTScript="macro.C"

# Root files without the .root extension
declare -a InputFileArray=("root_file_from_simulation_1" "root_file_from_simulation_2" )

# Where are the root files stored
InputDirectory="./Input"

# Over which files do we want to loop
declare -a RunOverFile=( 1 0 )

#Defining input and output directory with respect to current directory
for ((i=0;i<${#InputFileArray[@]};++i))
do
	if [ ${RunOverFile[i]} -eq 1 ]
	then
	#Getting the variables at the current iterator

	InputFile="${InputFileArray[i]}"


	echo "Currently working on: "$InputFile

	#Define output directories for PNGs and ROOT file
	OutputDirectory="Output"
	PNGOutputFolderName=$OutputDirectory"/PNGs/"$InputFile
	ROOTOutputFolderName=$OutputDirectory"/ROOT"

	#Create output directories (the -p option checks if the output directory already exists and does nothing in case they do)
	mkdir -p $PNGOutputFolderName
	mkdir -p $ROOTOutputFolderName

	root -l -b -q $ROOTScript'("'$InputFile'", "'$InputDirectory'", "'$OutputDirectory'")'
	fi
done

echo ""
echo "All done"
