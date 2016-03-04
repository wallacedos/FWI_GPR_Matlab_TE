#!/bin/bash

# Useage: qsub serial.sh
# Output: <JOB_NAME>.o<JOB_ID>
# lines begin with "#$" are to set qsub parameters.
# lines begin with "#" except "#!" and "#$" are comments.

#$ -S /bin/bash
#$ -cwd
#$ -m beas
#$ -j y
#$ -V

######### set job's name
#$ -N JOB_NAME

WORKPATH=`pwd`
echo "Current PATH = ${WORKPATH}"
echo "Begin computing, please have a rest ..."
#########  execute PROGRAM_NAME
cd ${WORKPATH}
matlab -nodisplay -r RTM03_181
