#!/bin/bash

for i in $(seq 2 1 5)
do
    resfile=results/${i}cD_results.txt

    if [ -f $resfile ]
    then
	rm -v $resfile
    fi

    dirs=$(ls -p | grep "/" | grep "D" | grep $i | grep "cD")
    for dir in $dirs
    do
	echo $dir | cut -d '/' -f 1
	files=$(ls $dir/ | grep "results")
	for file in $files
	do
	    size=$(echo $file | cut -d "_" -f 2 | cut -d '-' -f 2)
	    if [ $size -eq 10 ]
	    then
		cat $dir/$file | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000}' >> $resfile
	    else
		cat $dir/$file | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16" "$17" "$18" "$19" "$20" "$21" "$22" "$23" "$24" "$25" "$26}' >> $resfile
	    fi
	done
    done
done
