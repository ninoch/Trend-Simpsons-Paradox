#!/bin/bash

###############################
#       BINS PARAMETERS       #
###############################
num_of_bins=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['num_of_bins']")
least_num_of_datapoints=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['least_num_of_datapoints_in_each_bin']")
target_variable_column=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['target_variable_column']")
csv_file_name=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['csv_file_name']")

############################
#       FINDING BINS       #
############################

# make the model directory if it doesn't exist already
mkdir -p temporary_files
mkdir -p temporary_files/$csv_file_name

cd scripts
make
./cppexe -infile:../input/$csv_file_name".csv" -outfolder:../temporary_files/$csv_file_name -ycol:$target_variable_column -nbins:$num_of_bins -mindatapoints:$least_num_of_datapoints

#####################################
#       REGRESSION PARAMETERS       #
#####################################

mkdir -p ../output
mkdir -p ../store_results

python trend_simpson.py


printf "\nDONE\n"

