#! /bin/bash

# wrap the loon executable in a script for PBS submission
 
# when the job runs, the user is logged
# need to redo the minossoft setup here in order to pickup
# the srt_setup() function
echo "setting up minos soft $rel release"
source /stage/minos-data5/software/minossoft/setup/setup_minossoft_csf.sh ${rel}
 
pushd $HOME/MinosSoft/test/
echo "after pushd, path is " 
pwd
#set
srt_setup -a
popd
 
echo "after popd, path is " 
pwd
cd ${WORKDIR}
echo "cd to workdir: now in " 
pwd
export NUE_OUT_DIR="./"
echo "set NUE_OUT_DIR=$NUE_OUT_DIR"
echo "running job"
loon -bq $HOME/MinosSoft/test/NueAna/macros/MakeAnaNueTree.C ${file}
echo "all done with loon.  Contents of working dir are:"
ls -alhtr

echo "Copying output to ${finaldir}"
for f in AnaNue*; do
    mv $f ${finaldir}
done

echo "At end of job, we've got:"
ls -alhtr
echo "Goodbye and goodluck"