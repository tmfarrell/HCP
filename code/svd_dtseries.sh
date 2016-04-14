#!/bin/bash
# 
# Runs svd_dtseries.py for 1 'partition' of HCP subjects (for list of 
# partitions, see $base_dir/data/id_paritions.txt). Saves output and 
# error to $base_dir/data/out/QNL-HCP-svd_dtseries.o<ja_task_id>.SGE_TASK_ID
# and $base_dir/data/err/QNL-HCP-svd_dtseries.e<ja_task_id>.SGE_TASK_ID, respectively. 
#
# Tim Farrell, tmf@bu.edu
# QNL, BU
# 20160219

## for qsub
#$ -N QNL-HCP-svd_dtseries
#$ -o /projectnb/bohland/HCP/data/out/
#$ -e /projectnb/bohland/HCP/data/err/
#$ -m ae 
#$ -M tfarrell01@gmail.com 
# Note: change -M for email notifications.  

## prelims 
set -e
code_dir=/projectnb/bohland/HCP/code
if [ $1 == '-h' ]; then 
    echo -e "\nUsage:\t\t$ ./svd_dtseries.sh [--allow-recomputing] [--save-plots]"
    echo -e "\nParams:" 
    echo -e "--save-plots\tSaves all the preprocessing plots for all the runs." 
    echo -e "\t\t\tFalse by default." 
    echo -e "--allow-recomputing\tComputes all svds, even if subject-run svds already in $HCP/data/svds." 
    echo -e "\t\t\tFalse by default. Note: since by default svds already in $HCP/data/svds are not recomputed," 
    echo -e "\t\t\tthis will influence how many subjects are actually computed." 
    exit
fi 

## get optional args
args=""
if [ ! -z $1 ]; then
    args=$1
fi
if [ ! -z $2 ]; then
    args="$args $2"
fi

## main 
# compute svd for dtseries of subjects in partition SGE_TASK_ID
time python $code_dir/svd_dtseries.py $SGE_TASK_ID $args

# eof
