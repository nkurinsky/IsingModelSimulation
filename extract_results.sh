#!/bin/bash

dtypes="D- cD-"

pushd mcRuns

for i in $(seq 2 1 5)
do
    files=$(ls | grep "D-$i" | grep -v "cD" | grep "results")
    cfiles=$(ls | grep "D-$i" | grep "cD" | grep "results")
    cp -v $files ../${i}D_data/
    cp -v $cfiles ../${i}cD_data/
done

popd