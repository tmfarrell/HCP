#!/bin/bash
# 
# (a) Splits first $1 HCP subjects (taken from /projectnb/connectomedb/Q6 on scc1) 
#     into $2 partitions; and 
# (b) Parallelizes `svd_dtseries.sh` over $2 partitions using qsub. 
# 
# Tim Farrell, tmf@bu.edu 
# QNL, BU 
# 20160401

# prelims
set -e
if [ $1 == '-h' ]; then 
    echo -e "\nUsage:\n\t$ svd_dtseries_parallel.sh <#subjects> <#partitions> [--save-plots] [--allow-recomputing]"
    echo -e "\nParams:" 
    echo -e "#subjects\t\tNumber of subjects to run." 
    echo -e "#partitions\t\tNumber of partitions to run them over." 
    echo -e "--save-plots\t\tSaves all the preprocessing plots for all the runs. False by default." 
    echo -e "--allow-recomputing\tComputes all svds, even if subject-run svds already in $HCP/data/svds. False by default." 
    echo -e "\t\t\tNote: since by default svds already in $HCP/data/svds are not recomputed," 
    echo -e "\t\t\tthis will influence which subjects are and aren't computed (if you're computing less than all)."
    exit
fi 

# get #subjects we want to run
if [ ! -z $1 ]; then 
    n_subjects=$1
else
    echo "Need to pass #subjects (to run) as 1st arg."
    exit
fi

# get #partitions we want them to run them in
if [ ! -z $2 ]; then 
    n_partitions=$2
else 
    echo "Need to pass #partitions (to run those subjects over) as 2nd arg." 
    exit 
fi

# get optional args
args=""
if [ ! -z $3 ]; then
    args=$3
fi 
if [ ! -z $4 ]; then 
    args="$args $4" 
fi 

 
# init new file to store partitions
base_dir=/projectnb/bohland/HCP/data 
partitions_file=$base_dir/id_partitions.txt
if [ -f $partitions_file ]; then 
    rm -f $partitions_file
fi     
touch $partitions_file

# get #(subjects/partition) and array of subjects  
subjects_per_partition=$(echo "$n_subjects / $n_partitions" | bc) 
subjects=($(ls /projectnb/connectomedb/Q6/ | awk ' { print $1 } '))

# add #(subjects/partition) to each line of partitions file 
for p in $(seq 0 $(($n_partitions - 1))); do
    for s in $(seq 1 $subjects_per_partition); do
	i=$(( ($p * $subjects_per_partition) + $s ))
	if [ "1" == $s ]; then 
	    l=${subjects[$i]}
	else
	    l="$l\t${subjects[$i]}"
	fi 
    done
    echo -e $l >> $partitions_file
done

# clear error and output dirs
if [ "$(ls -A $base_dir/err/)" ]; then 
    rm -f $base_dir/err/*
fi 
if [ "$(ls -A $base_dir/out/)" ]; then 
    rm -f $base_dir/out/*
fi 

# parallelize svd_dtseries.sh over #partitions 
qsub -V -t 1-$n_partitions svd_dtseries.sh $args
