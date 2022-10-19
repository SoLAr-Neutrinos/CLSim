#!/bin/bash

#Define .C script that we want ROOT to interprete later on
ROOTScript="macro.C"

# Root files without the .root extension
declare -a SigInputFileArray=("single_electron_10_CL" "single_electron_5_CL" )
declare -a BackInputFileArray=("Ar39_background_CL" "Ar39_background_CL" )

# Where are the root files stored
InputDirectory="./Input"

# Over which files do we want to loop
declare -a RunOverFile=( 1 1 )

#Defining input and output directory with respect to current directory
for ((i=0;i<${#SigInputFileArray[@]};++i))
do
	if [ ${RunOverFile[i]} -eq 1 ]
	then
	#Getting the variables at the current iterator

	InputFile="${InputFileArray[i]}"
	Sig=${SigInputFileArray[i]}
	Back=${BackInputFileArray[i]}

	echo "Currently working on: "$InputFile

	#Define output directories for PNGs and ROOT file
	OutputDirectory="Output"
	PNGOutputFolderName=$OutputDirectory"/PNGs/"$Sig

	rm -rf $PNGOutputFolderName

	ROOTOutputFolderName=$OutputDirectory"/ROOT"

	#Create output directories (the -p option checks if the output directory already exists and does nothing in case they do)
	mkdir -p $PNGOutputFolderName
	mkdir -p $ROOTOutputFolderName

	root -l -b -q $ROOTScript'("'$Sig'", "'$Back'","'$InputDirectory'", "'$OutputDirectory'")'
	fi
done

echo ""
echo "All done"
