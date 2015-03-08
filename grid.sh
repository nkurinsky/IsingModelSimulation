#!/bin/bash

#sizes="10 30 100"
sizes="50"
dims="2 3 4 5"
qs="2 3 5 10"

script="./run_T_grid.sh"

for dim in $dims
do
    for q in $qs
    do
	for size in $sizes
	do
	    echo $script $dim $size $q
	    bsub -q long $script $dim $size $q
	done
    done
done
