#!/bin/sh

cd /data/minos/NueAnaPID
#setup_minos
source /minosbase/minossoft/setup/setup_minossoft_macfink.sh R2.0.1
srt_setup -a
cd -


#requires name of script to be same as name of directory
runit()
{
	echo $1
	cd $1
	rm -f log
	rm -f out.root
	loon -b -q /data/minos/NueAnaPID/NueAna/macros/LoadLibs.C  $1.C >&log
	date
	cd ..
}

date
date > runtime


#######

runit near_mc_mrcc_hornON_normal
runit near_mc_standard_hornON_normal
runit near_data_mrcc_hornON_normal
runit near_data_standard_hornON_normal


runit near_mc_mrcc_hornON_ParticlePID
runit near_mc_standard_hornON_ParticlePID
runit near_data_mrcc_hornON_ParticlePID
runit near_data_standard_hornON_ParticlePID



#######near mc horn off

runit near_mc_mrcc_hornOFF_normal
runit near_mc_standard_hornOFF_normal
runit near_mc_standard_hornOFF_ParticlePID


#######near data horn off

runit near_data_mrcc_hornOFF_normal
runit near_data_standard_hornOFF_normal
runit near_data_standard_hornOFF_ParticlePID


######far
runit far_mc_mrcc_normal
runit far_mc_mrcc_ParticlePID
runit far_mc_standard_normal
runit far_mc_standard_ParticlePID



###############

root -b -q GatherDirectories.C


date
date >>runtime
