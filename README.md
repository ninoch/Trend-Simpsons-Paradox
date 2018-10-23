# What does it do? 

We describe a data-driven discovery method that leverages Simpson's paradox to uncover interesting patterns in behavioral data. Our method systematically disaggregates data to identify subgroups within a population whose behavior deviates significantly from the rest of the population. 
Given an outcome of interest and a set of covariates, the method follows three steps. First, it disaggregates data into subgroups, by conditioning on a particular covariate, so as minimize the variation of the outcome within the subgroups. Next, it models the outcome as a linear function of another covariate, both in the subgroups and in the aggregate data. Finally, it compares trends to identify disaggregations that produce subgroups with different behaviors from the aggregate.

More details about the algorithm is available in our paper on ICWSM 2018: "[Using Simpson’s Paradox to Discover Interesting Patterns in Behavioral Data](https://arxiv.org/abs/1805.03094)". 
Please cite our paper if you use this code in your own work:

	@article{alipourfard2018using,
	  title={Using Simpson's Paradox to Discover Interesting Patterns in Behavioral Data},
	  author={Alipourfard, Nazanin and Fennell, Peter G and Lerman, Kristina},
	  booktitle={Twelfth International AAAI Conference on Web and Social Media (ICWSM)},
	  year={2018}
	}


# Input and Running: 
For running the algorithm you need to follow three steps:
#### Step 1: 
Add .csv data file to input/ directory.

	- Put your data file (.csv) into input/ directory
	- Your data file must consists of the header row for name of the variables. Each row is a datapoint.

#### Step 2: 
Update input_info.json file.

	- num_of_bins: Maximum number of bins for disaggregating data
	- least_num_of_datapoints_in_each_bin: Minimum number of datapoints in each bin
	- target_variable: name of the target variable in input .csv file. The Y variable **must** be binary in this version. 
	- target_variable_column: Column number of target variable in .csv file (zero-based).
	- level_of_significance: Level of significance for chi-square deviance test
	- csv_file_name: Name of the file you put to the input/ directory in step 1
	- ignore_columns: An array of name of the variables for not including them in the algorithm variables. You should list all the columns with string or float values.
	- log_scales: A list of variable name, which shows you prefer log scale for axis for that variable in output/ plots.

Before running the algorithm, please make sure the json format of the input_info.json file is [valid](https://jsonformatter.curiousconcept.com/).

#### Step 3: 
Run the run.sh script.

# Output: 
After running, three directories will be added. 

	- temporary_files/
	- output/
	- store_results/

The information about the bins, will be available in temporary_files/ directory. 

	- bins.csv
	  	- L rows (L = number of features/covariates)
	  	- Row l is the bin splits for feature L
	  	- the bins then are given by the intervals
			- [(l,1) (l,2)]
			- ((l,2) (l,3)]
			- ((l,3) (l,4)]
			…
			- ((l,nl-1) (l,nl)]
	    where (l,i) is the (i+1)’st entry of row l (with 0 indexing, element (l,0) is the name of the feature.
	- N_tree.csv: Number of datapoints in each bin
	- R2improvements.csv: Total R^2 improvement for each variable. 
	- ybar_tree.csv: Average value of Y in each bin. 


Plots of Trend Simpson's Pairs will be available in output/ directory. For each pair, there is a PDF file. First plot is logistic fit to aggregated data. Second one is logistic for each of the bins. The third plots are histogram and heatplot of the Paradox and Conditioning variables. 


The details of the algorithm will be printed in the terminal. You can use them as log. Beside that, all the informations about the logistic fits, simpson's pairs and deviance values are availabe in store_results/ directory as python pickle objects. You can use load functions in trend_simpsons.py file for loading them. 

	- aggregated_vars_params.obj: Logistic fit information (Intercept, coefficients, errors) for aggregated data for different variables.
	- disaggregated_vars_params.obj: For every pair of variables, information of logistic fit to disaggregated data. 
	- simpsons_pairs.obj: Trend Simpson's Paradox pairs (before finalizing them by chi-squared test).
	- loglikelihoods.obj: Loglikelihoods for full / null, aggregated / disaggregated models, for using to compute deviance of models and chi-squared test. 


#### Notes:
Please keep in mind that since we are using c++, the end-of-line is really important. Since [the end-of-line depends on OS](https://www.loginradius.com/engineering/eol-end-of-line-or-newline-characters/), if your data is generated in different OS, and you are runnig this code on another OS, consider using "[cleaning_end_of_lines.py](https://github.com/ninoch/Trend-Simpsons-Paradox/blob/master/scripts/cleaning_end_of_lines.py)" script to clean up the end-of-lines. 

