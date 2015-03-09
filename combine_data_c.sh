#!/bin/bash

dirs=$(ls -p | grep "/" | grep "2cD")
for dir in $dirs
do
    files=$(ls $dir/ | grep "results")
    for file in $files
    do
	size=$(echo $file | cut -d "_" -f 2 | cut -d '-' -f 2)
	if [ $size -eq 10 ]
	then
	    cat $dir/$file | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000" "0.000000}'
	else
	    cat $dir/$file | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16" "$17" "$18" "$19" "$20" "$21" "$22" "$23" "$24" "$25" "$26}'
	fi
    done
done
