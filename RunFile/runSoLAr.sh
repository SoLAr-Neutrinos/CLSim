for f in ../Input/SoLAr111/*
do
	echo $f
	Outfile=`basename $f`
	Outfile=$(echo $Outfile | sed 's/.root/_denseC_L.root/')
	Outfile=$Outfile

	Infofile=$(echo $f | sed 's/\.\.\///')
	echo $Infofile

	stud=`basename $f`

	run=run_$stud
	echo "cd /pc2014-data3/tdieminger && source setup_lAr.sh && cd /pc2014-data3/tdieminger/SoLAr_Env/CLSim && ./analyze_light $Infofile PlacementFiles/26161/SiPM 0.6 $Outfile --charge PlacementFiles/26161/densePixels --pixSize 1 --number 10000" > $run.sh

	chmod +x $run.sh
	echo "Submitting $run"
	sleep 1
	 while [[ $( screen -ls | tail -2 | head -1 | awk '{print $1}') -gt 40 ]]
	 do
	 	echo "Waiting for a screen to free up - 60 seconds"
	 	sleep 60
	 done
	screen -S $run -m -d bash $run.sh
done
