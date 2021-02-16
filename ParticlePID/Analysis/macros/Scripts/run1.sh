#!/bin/sh

RELEASE=/minos/app/scavan/Systematics

args=(${2//^/ })
saveDir=${args[2]}
type=${args[0]}
sys=${args[1]}


tmpDir=${_CONDOR_SCRATCH_DIR}
if [ -z $tmpDir ]; then
	tmpDir=$PWD
fi

export MINOS_SETUP_DIR=/afs/fnal.gov/files/code/e875/general/minossoft/setup
echo SETTING UP UPS

unset SETUP_UPS SETUPS_DIR
. /afs/fnal.gov/files/code/e875/general/ups/etc/setups.sh


setup_minos()
{
. $MINOS_SETUP_DIR/setup_minossoft_FNALU.sh $*
}

setup_minos -r S09-09-18-R2-00    #S09-07-10-R2-00

echo "release $RELEASE"
echo "type $type"
echo "sys $sys"
echo "tmpdir $tmpDir"

cd $RELEASE
srt_setup -a


cd $tmpDir


#requires name of script to be same as name of directory
runit()
{
	echo $1 $2 $3

	rm -f log
	rm -f out.root
	rm -f runtime
	rm -f filelist
	#rm -f memuse

	date
	date > runtime

	#need to currently manually kill the job!
	#$RELEASE/NueAna/ParticlePID/Analysis/macros/Scripts/memcheck.sh >memuse&

	loon -b -q $RELEASE/NueAna/macros/LoadLibs.C  "$RELEASE/NueAna/ParticlePID/Analysis/macros/Scripts/$1/$1.C(\"$2\",$3)" >&log
	date
	date >>runtime

	mkdir -p $saveDir/$3/$1
	
	mv log $saveDir/$3/$1
	mv out.root $saveDir/$3/$1
	mv runtime $saveDir/$3/$1
	mv filelist $saveDir/$3/$1
	
	#mv memuse $saveDir/$3/$1

	chmod -R 777 $saveDir/$3/$1
	ls

	

	exit
}






######input file locations
files_near_mc_mrcc_hornON_normal=/minos/data/users/scavan/Dogwood1_Near_ANN14_newMRCC/out/minos/data2/nue_group_tmp/MRCC/MC/L010185N
files_near_mc_standard_hornON_normal=
files_near_data_mrcc_hornON_normal=
files_near_data_standard_hornON_normal=
files_near_mc_mrcc_hornON_ParticlePID=
files_near_mc_standard_hornON_ParticlePID=
files_near_data_mrcc_hornON_ParticlePID=
files_near_data_standard_hornON_ParticlePID= 
files_near_mc_mrcc_hornOFF_normal=
files_near_mc_standard_hornOFF_normal=
files_near_mc_standard_hornOFF_ParticlePID= 
files_near_data_mrcc_hornOFF_normal=
files_near_data_standard_hornOFF_normal= 
files_near_data_standard_hornOFF_ParticlePID= 
files_far_mc_mrcc_normal= 
files_far_mc_mrcc_ParticlePID=
files_far_mc_standard_normal=
files_far_mc_standard_ParticlePID=


#######
if [ "$type" = "1" ]; then
runit near_mc_mrcc_hornON_normal $files_near_mc_mrcc_hornON_normal $sys
fi

if [ "$type" = "2" ]; then
runit near_mc_standard_hornON_normal $files_near_mc_standard_hornON_normal $sys
fi

if [ "$type" = "3" ]; then
runit near_data_mrcc_hornON_normal $files_near_data_mrcc_hornON_normal $sys
fi

if [ "$type" = "4" ]; then
runit near_data_standard_hornON_normal $files_near_data_standard_hornON_normal $sys
fi

if [ "$type" = "5" ]; then
runit near_mc_mrcc_hornON_ParticlePID $files_near_mc_mrcc_hornON_ParticlePID $sys
fi

if [ "$type" = "6" ]; then
runit near_mc_standard_hornON_ParticlePID $files_near_mc_standard_hornON_ParticlePID $sys
fi

if [ "$type" = "7" ]; then
runit near_data_mrcc_hornON_ParticlePID $files_near_data_mrcc_hornON_ParticlePID $sys
fi

if [ "$type" = "8" ]; then
runit near_data_standard_hornON_ParticlePID $files_near_data_standard_hornON_ParticlePID $sys
fi


#######near mc horn off
if [ "$type" = "9" ]; then
runit near_mc_mrcc_hornOFF_normal $files_near_mc_mrcc_hornOFF_normal $sys
fi

if [ "$type" = "10" ]; then
runit near_mc_standard_hornOFF_normal $files_near_mc_standard_hornOFF_normal $sys
fi

if [ "$type" = "11" ]; then
runit near_mc_standard_hornOFF_ParticlePID $files_near_mc_standard_hornOFF_ParticlePID $sys
fi

#######near data horn off
if [ "$type" = "12" ]; then
runit near_data_mrcc_hornOFF_normal $files_near_data_mrcc_hornOFF_normal $sys
fi

if [ "$type" = "13" ]; then
runit near_data_standard_hornOFF_normal $files_near_data_standard_hornOFF_normal $sys
fi

if [ "$type" = "14" ]; then
runit near_data_standard_hornOFF_ParticlePID $files_near_data_standard_hornOFF_ParticlePID $sys
fi



######far
if [ "$type" = "15" ]; then
runit far_mc_mrcc_normal $files_far_mc_mrcc_normal $sys
fi

if [ "$type" = "16" ]; then
runit far_mc_mrcc_ParticlePID $files_far_mc_mrcc_ParticlePID $sys
fi

if [ "$type" = "17" ]; then
runit far_mc_standard_normal $files_far_mc_standard_normal $sys
fi

if [ "$type" = "18" ]; then
runit far_mc_standard_ParticlePID $files_far_mc_standard_ParticlePID $sys
fi


###############






