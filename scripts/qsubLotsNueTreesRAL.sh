#! /bin/bash

OUTPUTD="/stage/minos-data6/vahle/NueTrees/"


if [ $1 == mc ]; then
    echo $1
    OUTPUTDIR=${OUTPUTD}/near_mc
    echo $2
    if [ $2 == "L010185" ]; then
	echo "LE10 185"
#	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'n130[1-2]1*L010185.sntp.R1_18_2.root')
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'n13011*L010185.sntp.R1_18_2.root')
    elif [ $2 == "L010170" ]; then
	echo "LE10 170"
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'n130[1-2]1*L010170.sntp.R1_18_2.root')
    elif [ $2 == "L010200" ]; then
	echo "LE10 200"
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'n130[1-2]1*L010200.sntp.R1_18_2.root')
    elif [ $2 == "L100200" ]; then
	echo "ME 200"
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'n1301101*L100200.sntp.R1_18_2.root')
    elif [ $2 == "L250200" ]; then
	echo "HE 200"
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'n13011*L250200.sntp.R1_18_2.root')
    elif [ $2 == "L010000" ]; then
	echo "Horn off"
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'n13011*L010000.sntp.R1_18_2.root')
    elif [ $2 == "farbeam" ]; then
	echo "far"
	OUTPUTDIR=${OUTPUTD}/far_mc
	echo $OUTPUTDIR
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'f210*L010185.sntp.R1_18_2.root')
    elif [ $2 == "farnue" ]; then
	echo "far"
	OUTPUTDIR=${OUTPUTD}/far_mc
	echo $OUTPUTDIR
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'f211*L010185.sntp.R1_18_2.root')
    fi
    
elif [ $1 == "data" ]; then
    echo $1
    echo $2
    OUTPUTDIR=${OUTPUTD}/near_data
    echo $OUTPUTDIR
    # all 1e20 nd plus a bit
    if [ $2 == "all" ]; then
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'N000*.spill.sntp.R1_18_2.0.root')
    elif [ $2 == "L010000" ]; then
	filelist=$(find /stage/minos-data6/vahle/hornoff-sntp-data -name 'N000*.spill.sntp.R1_18_2.0.root')
    elif [ $2 == "some" ]; then
	filelist=$(find /stage/minos-data3/dcm_catalogue -name 'N0000[7-8]*.spill.sntp.R1_18_2.0.root')
    fi
elif [ $1 == "test" ]; then
    filelist=/stage/minos-data6/nearMC_R1_18.2/sntp/n1301100[3-7]_0000_L010185.sntp.R1_18_2.root
fi

for f in $filelist; do
#    echo $f
    $HOME/MinosSoft/test/NueAna/scripts/qsubOneNueTreeRAL.sh $f $OUTPUTDIR
done