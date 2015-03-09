#!/bin/bash

echo "Extracting"
./extract_results.sh
echo "Combining"
./combine_data.sh
./combine_data_c.sh
