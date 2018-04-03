#!/bin/bash

###############################
#       BINS PARAMETERS       #
###############################

num_of_bins=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['num_of_bins']")
least_num_of_datapoints=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['least_num_of_datapoints_in_each_bin']")
target_variable_column=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['target_variable_column']")
csv_file_name=$(cat "input_info.json" | python -c "import sys, json; print json.load(sys.stdin)['csv_file_name']")


##################################
#       MAKING OUT FOLDERS       #
##################################

# make the model directory if it doesn't exist already

printf "\n *** BASH: Making folders for outputs\n"

mkdir -p temporary_files
mkdir -p temporary_files/$csv_file_name
mkdir -p output
mkdir -p store_results

############################
#       FINDING BINS       #
############################

printf "\n *** BASH: Run the binning algorithm and store results in temporary_files/ \n"

cd scripts
make
./cppexe -infile:../input/$csv_file_name".csv" -outfolder:../temporary_files/$csv_file_name -ycol:$target_variable_column -nbins:$num_of_bins -mindatapoints:$least_num_of_datapoints

######################################
#       FINDING SIMPSONS PAIRS       #
######################################

printf "\n *** BASH: Run Trend Simpson's Paradox algorithm to find pairs and store results in output/ \n"

python trend_simpson.py 2> /dev/null 

printf "\n *** BASH: DONE\n"

