#!/bin/bash
# 
# (a) Splits first $1 HCP subjects (taken from /projectnb/connectomedb/Q6 on scc1) 
#     into $2 partitions; and 
# (b) Parallelizes `svd_dtseries.sh` over those partitions using qsub. 
# 
# Run `svd_dtseries_parallel.sh -h` for usage instructions. 
# 
# Tim Farrell, tmf@bu.edu 
# QNL, BU 
# 20160401

# prelims
set -e
if [ "$1" == "-h" ]; then 
    echo -e "\nUsage:\n\t$ svd_dtseries_parallel.sh <#subjects> <#partitions> [args]"
    echo -e "\nParams:" 
    echo -e "#subjects\tNumber of subjects to run." 
    echo -e "#partitions\tNumber of partitions to run them over." 
    echo -e "args\t\tArgs for svd_dtseries.py" 
    echo -e ""
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
# passing optional args (${@:3}) 
qsub -V -t 1-$n_partitions svd_dtseries.sh "${@:3}" 