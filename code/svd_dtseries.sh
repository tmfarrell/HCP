#!/bin/bash
# 
# Runs svd_dtseries.py for partition $SGE_TASK_ID in $HCP/data/id_partitions.txt. 
# Saves output and error to $HCP/data/out/QNL-HCP-svd_dtseries.o<ja_task_id>.$SGE_TASK_ID
# and $HCP/data/err/QNL-HCP-svd_dtseries.e<ja_task_id>.$SGE_TASK_ID, respectively. 
#
# Tim Farrell, tmf@bu.edu
# QNL, BU
# 20160219

## args for qsub
#$ -N QNL-HCP-svd_dtseries
#$ -o /projectnb/bohland/HCP/data/out/
#$ -e /projectnb/bohland/HCP/data/err/
#$ -m ae 
#$ -M tfarrell01@gmail.com 
# See http://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/
# for more on qsub args.   

## prelims 
set -e
code_dir=/projectnb/bohland/HCP/code
if [ "$1" == "-h" ]; then 
    echo -e "\nUsage:\t\t$ ./svd_dtseries.sh [args]"
    echo -e "\nParams:" 
    echo -e "args\t\tArgs for svd_dtseries.py (e.g. --allow-recomputing --save-plots)." 
    exit
fi

## main 
# compute svd for dtseries of subjects in partition SGE_TASK_ID
# passing in optional args (${@:1}) 
time python $code_dir/svd_dtseries.py "${@:1}" $SGE_TASK_ID

# plot results 
# TODO: need to add cmdline args to plot_results.py
#       so that it only plots those svds from this partition
#python plot_results.py $SGE_TASK_ID

# eof
