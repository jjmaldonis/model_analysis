#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -N vor_rmc
#$ -q all.q
#$ -V
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#
/home/jjmaldonis/OdieCode/vor/hrmc/vor_hrmc $JOB_ID
