#!/bin/bash

if [ $# -gt 2 ]
then

    dim=$1
    size=$2
    q=$3
    
    tmin=0.5
    tmax=4.0
    tstep=0.1
    
    outfilebase="D-"$dim"_s-"$size"_q-"$q"_T-"
    
    newdir="D-"$dim"_s-"$size"_q-"$q
    if [ -d $newdir ]
    then
	rm -rf $newdir
    fi
    mkdir $newdir
    cd $newdir
    
    outfileroot=".out"
    resfile="../"$newdir"_results.txt"
    cmd="PottsMC $dim $size $q"
    for i in $(seq $tmin $tstep $tmax)
    do
	outfile=$outfilebase$i$outfileroot
	scriptfile=$outfilebase$i$scriptfileroot
	echo $cmd $i $outfile
	$cmd $i $outfile
	tail -n 1 $outfile >> $resfile
    done
else
    echo "Not enough arguments"
    echo "Calling Sequence: "$0" dim size q"
fi    
