#! /bin/bash

#OUTPUTDIR="/stage/minos-data6/vahle/NueTrees/"
RELEASE="new"

qsub -M${USER}@hep.ucl.ac.uk -m a -j oe -v "file=$1, finaldir=$2, rel=$RELEASE" $SRT_PRIVATE_CONTEXT/NueAna/scripts/runNueTreesRAL.sh