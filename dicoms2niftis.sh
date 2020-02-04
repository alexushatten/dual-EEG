#!/bin/bash
#1st argument: subject ID
#2nd argument: data type (anat, dwi, func, fmap)
#3rd argument: input directory
#4th argument: output directory
#filepath of script with filename and extension

#DIR="$( dirname "${0}" )" # Get the directory where this script is stored

DIR=$4
SUBJ_DIR=$1
DATA_TYPE_DIR=$2
INPUT_DIR=$3
ID=$1

BIDS_DIR="${DIR}/${SUBJ_DIR}/${DATA_TYPE_DIR}"

mkdir -p $BIDS_DIR

#Run command for conversion of DICOM to NIfTI using dcm2niix

echo ============================================= >> $DIR/log.txt
echo dcm2niix -b y -z y -v y -o $BIDS_DIR -f "${ID}_%d_%t_${DATA_TYPE_DIR}" $INPUT_DIR >> $DIR/log.txt
dcm2niix -b y -z y -v y -o $BIDS_DIR -f "${ID}_%d_%t_${DATA_TYPE_DIR}" $INPUT_DIR >> $DIR/log.txt
cat $DIR/log.txt