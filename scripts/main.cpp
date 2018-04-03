#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <string>
#include <algorithm>

// Have multiple conditions

using namespace std;

const int DL_thresh = 3;    // threshold for loss function

#include "main.h"

int main(int argc, char* argv[])
{
    
    cout << "\n* BINNING ALGORITHM *\n" << endl;
    
    clock_t start = clock();
    
    // ------------------------------
    // 1: Read in inputs
    
    
    ifstream datafile;      // input .csv file
    string outfolder;       // folder where to store the output data
    int N = 0, N_declared=0;    // Number of rows of .csv file
    int y_col=0;            // Position of ycolumn (default 1st column)
    double lambda = 0.000000001; // default lambda parameter
    int nbins = 10000, mindatapoints = 0; // number of bins, minimumn number of datapoints per bin
    string lambda_val, nbins_val, mindatapoints_val;
    bool lambda_stopping = true, nbins_stopping = false, mindatapoints_stopping = false; // is this stopping method to be used for the binning.
    
    for(int i=0; i<argc; i++)
    {
        if(strncmp(argv[i],"-infile:",8)==0){
            datafile.open(string(argv[i]).erase(0,8).c_str());
            cout << " Reading data from " << string(argv[i]).erase(0,8) << endl;
        }
        if(strncmp(argv[i],"-outfolder:",11)==0){
            outfolder = string(argv[i]).erase(0,11);
            cout << " Writing data to folder " << outfolder << endl;
        }
        if(strncmp(argv[i],"-ycol:",6)==0)
            y_col = atol(string(argv[i]).erase(0,6).c_str());
        if(strncmp(argv[i],"-n_rows:",8)==0)
        {
            N_declared = 1;
            N = atol(string(argv[i]).erase(0,8).c_str());
        }
        if(strncmp(argv[i],"-lambda:",8)==0)
        {
            lambda = atof(string(argv[i]).erase(0,8).c_str());
            lambda_val = to_string(lambda);
            lambda_stopping = true;
            cout << " - lambda = " << lambda << endl;
        }
        if(strncmp(argv[i],"-nbins:",7)==0)
        {
            nbins = atoi(string(argv[i]).erase(0,7).c_str());
            nbins_val = to_string(nbins);
            nbins_stopping = true;
            cout << " - Number of bins = " << nbins << endl;
        }
        if(strncmp(argv[i],"-mindatapoints:",15)==0)
        {
            mindatapoints = atoi(string(argv[i]).erase(0,15).c_str());
            mindatapoints_val = to_string(mindatapoints);
            mindatapoints_stopping = true;
            cout << " - Minimum number of datapoints per bin = " << mindatapoints << endl;
        }
    }
    if(!datafile.is_open())
    {
        cout << "ERROR: datafile not found" << endl;
        return 1;
    }
    
    // ---------------------------------------------------
    // 2: Read in dependent and independent variable names
    
    cout << "\nFEATURES:" << endl;
    vector<string> features;
    get_column_names(datafile, features);
    int N_columns = features.size();
    for(auto it =features.begin(); it!=features.end(); ++it)
        if(*it != features[y_col])
            cout << *it << endl;
    cout << endl;
    
    cout << "TARGET VARIABLE:\n" << features[y_col] << endl << endl;
    
    cout << "Reading in Y data..." << endl;
    double y, ybar=0, SST=0, R2=0;
    if(N_declared == 0) // if the number of rows has not been passed as as input
        while(!datafile.eof())
        {
            N++; // N is earlier initialized to 0
            read_single_column(y_col, datafile, y);
            ybar += y;
            SST += y*y;
        }
    else // if the number of rows has been declared
        for(int i=0; i<N; i++)
        {
            read_single_column(y_col, datafile, y);
            ybar += y;
            SST += y*y;
        }
    
    ybar /= N;
    SST -= N*ybar*ybar; // SST = \sum y_i^2 - N*ybar^2
    
    cout << " - " << N << " rows of data\n - SampleMean(Y) = " << ybar << "\n - SampleSTD(Y) = " << sqrt(1.0*SST/(N-1)) << "\n - SST(Y) = " << SST<< endl;
    
    // ------------------------------
    // 3: Open outfiles & add headers
    
    ofstream bins_file(string(outfolder).append("/bins.csv").c_str());
    if(!bins_file.is_open())
    {
        cout << "ERROR: Outfolder files are not open." << endl;
        return 1;
    }
    ofstream N_file(string(outfolder).append("/N_tree.csv").c_str());
    if(!N_file.is_open())
    {
        cout << "ERROR: Outfolder files are not open." << endl;
        return 1;
    }
    ofstream ybar_file(string(outfolder).append("/ybar_tree.csv").c_str());
    if(!ybar_file.is_open())
    {
        cout << "ERROR: Outfolder files are not open." << endl;
        return 1;
    }
    ofstream R2improvements_file(string(outfolder).append("/R2improvements.csv").c_str());
    if(!R2improvements_file.is_open())
    {
        cout << "ERROR: Outfolder files are not open." << endl;
        return 1;
    }

    // ----------------------------
    // 4: Declare variables
    
    // Groupings
    int G=1; // number of groups
    vector<unsigned int> groups(N,0); // group number of each row of data
    
    // Cumulative
    vector<map<double, unsigned int> > cumcount(G,map<double, unsigned int>());
    vector<map<double, double> > cumsum(G,map<double, double>());;
    
    // Feature splitting variables
    double R2_improvement, best_R2_improvement, R2_improvement_s;
    vector<double> R2_improvements, DL;
    set<double> range;
    vector<double> splits;
    map<double, double> deltaR2_p, split_p;
    double pl, pu, pl1, pu1, pl2, pu2, split;
    int mindatapoints_feature; // minimum number of datapoints for a feature as the algorithm progresses
    
    // --------------------------
    //  5: BINNING
    // --------------------------
    
    cout << "\nCOMMENCING BINNING" << endl << endl;
    
    for(int j=0; j<N_columns; j++) // for each feature
    {
        
        // 0: check if column is good
        if(j == y_col) // if this it the y_columns
            continue;
        /*if(bad_feature[j]==1) // if its not a wanted feature
        {
            cout << " - skipping bad feature " << features[j] << endl << endl;
            continue;
        }*/
        
        cout << " - processing " << features[j] << endl;
        cout << "   - reading data..." << endl;
        
        // i: read in the feature column from file to create cumcount and cumsum vectors
        create_Xj_info(j, cumcount, cumsum, range, datafile, groups, G, y_col, N);
        // error check on the ybars (deleted as comparing floating point numbers (not good)
        
        cout << "   - range: [" << *(range.begin()) << "," << *(std::prev(range.end())) << "]" << endl;
        cout << "   - splits:" << endl;
        cout << "     - (interval, splits, R2)" << endl;
        
        // ii: initialize the  R2 data and the deltaR2_p map
        R2_improvement = 0;
        R2_improvements.clear();
        deltaR2_p.clear();
        split_p.clear();
        splits.clear();
        splits.push_back(*(range.begin())); // first and last values
        splits.push_back(*(std::prev(range.end())));
        DL.clear();
        
        // iii: find the best first split
        pl = *(range.begin());
        pu = *(std::prev(range.end())); // pl and pu are lower and upper bounds of the partition
        calculate_best_split(cumcount, cumsum, pl, pu, range, G, split, R2_improvement_s);
        best_R2_improvement = R2_improvement_s/SST;
        R2_improvement = best_R2_improvement;
        update_R2_DL(best_R2_improvement, R2_improvements, R2, DL, lambda);
        
        // iv: create the two subpartitions from the best split, [Xl,split] and (split,Xu)
        splits.push_back(split);
        pl2 = *(std::next(range.find(split))); // next element after split
        pu2 = pu;
        pl1 = pl;
        pu1 = split;
        mindatapoints_feature = min(cumcount[0][pu1], cumcount[0][pu2] - cumcount[0][pu1]);
            
        
        // iii: calculate the rest of the splits until the stopping condition has been reached
        for(int nbins_feature = 2; ;nbins_feature ++)
        {
            
            // 0: check to see if stopping condition is reached
            if(lambda_stopping == true) // lambda optimization
                if(lambda_stopping_condition(DL, best_R2_improvement, R2_improvements, splits, R2_improvement))
                {
                    cout << "     - Loss function has reached minimum (neglect last " << DL_thresh << " splits); breaking" << endl;
                    break;
                }
            if(nbins_stopping == true) // number of bins optimization
                if(nbins_feature > nbins)
                {
                    splits.pop_back(); // remove the last split as it has resulted in a partition with more than nbins
                    cout << "     - nbins = " << nbins << "; breaking" << endl;
                    break;
                }
            if(mindatapoints_stopping == true) // min number of datapoints optimization
            {
                if(mindatapoints_feature > min(cumcount[0][pu1], cumcount[0][pu2] - cumcount[0][pu1]))
                    mindatapoints_feature = min(cumcount[0][pu1], cumcount[0][pu2] - cumcount[0][pu1]);
                if(mindatapoints_feature < mindatapoints)
                {
                    splits.pop_back(); // remove the last split as it has resulted in a partition with less than mindatapoints
                    cout << "     - mindatapoints: " << mindatapoints_feature << " datatpoints in next split (< " << mindatapoints << "); breaking" << endl;
                    break;
                }
            }
            // output the information
            cout << "     - [" << pl << "," << pu << "] split into [" << pl << "," << split << "] and [" << pl2 << "," << pu << "], ";
            cout << "R2 = " << R2_improvement << endl; //<< ", DL = " << DL.back() << endl;

            // a: change cumcount and cumsum for the split (only have the change the second partition
            recalculate_partitioned_cums(cumcount, cumsum, pl1, pl2, pu2, G);
            
            // b: calculate the best splits and R2 improvements for the splitted partition
            calculate_best_split(cumcount, cumsum, pl1, pu1, range, G, split, R2_improvement_s);
            deltaR2_p[pu1] = R2_improvement_s/SST;
            split_p[pu1] = split;
            
            calculate_best_split(cumcount, cumsum, pl2, pu2, range, G, split, R2_improvement_s);
            deltaR2_p[pu2] = R2_improvement_s/SST;
            split_p[pu2] = split;
            
            // c: calculate the best partition and split overall (start from second element, first is placeholder);
            best_R2_improvement = deltaR2_p.begin()->second;
            split = split_p.begin()->second;
            pl = *(range.begin());
            pu = deltaR2_p.begin()->first;
            for(auto it = std::next(deltaR2_p.begin()); it!=deltaR2_p.end(); ++it)
            {
                if(it->second > best_R2_improvement)
                {
                    best_R2_improvement = it->second;
                    split = split_p[it->first];
                    pl = *std::next(range.find(std::prev(it)->first)); //  the bottom of the this interval is the next element after the end of the last interval
                    pu = it->first;
                }
            }
            
            if(best_R2_improvement <= 0) // if there are no R2 improvements
            {
                cout << "     - No further possible improvements in R2; breaking" << endl;
                break;
            }
            
            // d: create the two subpartitions from the best split
            splits.push_back(split);
            pl2 = *(std::next(range.find(split))); // next element after split
            pu2 = pu;
            pl1 = pl;
            pu1 = split;
            
            // e: update the R2 improvement and DL
            update_R2_DL(best_R2_improvement, R2_improvements, R2, DL, lambda);
            R2_improvement += best_R2_improvement;
            
        } // splits calculated
        
        /*
        cout << "     - [" << pl << "," << pu << "] split into [" << pl << "," << split << "] and [" << pl2 << "," << pu << "], R2 = " << R2_improvement << endl; // ", DL = " << DL.back() << endl;
        */
        cout << "   - Chosen partitions: ";
        std::sort(splits.begin(), splits.end());
        cout << "[" << *splits.begin() << "," << *std::next(splits.begin()) << "], ";
        for(auto it = std::next(std::next(splits.begin())); it!=splits.end(); ++it)
            cout << "(" << *std::prev(it) << "," << *it << "], ";
        cout << endl;
        cout << "   - R2 improvement: " << R2_improvement << endl << endl;
        
        // sort the splits
        std::sort(splits.begin(), splits.end());

        // f: Save all the info
        
        //  - i) save the splits
        bins_file << features[j] << ',';
        for(auto it = splits.begin(); it!= std::prev(splits.end()); ++it)
            bins_file << *it << ',';
        bins_file << *(std::prev(splits.end())) << endl;
        
        
        //  - ii) save the Number of elements per partition
        N_file << features[j] << ',';
        for(auto it = std::next(splits.begin()); it!= std::prev(splits.end()); ++it)
            N_file << cumcount[0][*it] << ',';
        N_file << cumcount[0][*(std::prev(splits.end()))] << endl;
        
        //  - iii) save the Number of elements per partition
        ybar_file << features[j] << ',';
        for(auto it = std::next(splits.begin()); it!= std::prev(splits.end()); ++it)
            ybar_file << cumsum[0][*it]/(1.0*cumcount[0][*it]) << ',';
        ybar_file << cumsum[0][*(std::prev(splits.end()))]/(1.0*cumcount[0][*(std::prev(splits.end()))]) << endl;
        
        //  - iv) save the R2 improvements
        R2improvements_file << features[j] << ',' << R2_improvement << endl;

        
    } // end of loop over all features

    clock_t end = clock();
    cout << "\nBINNING COMPLETE\nTotal procedure time = " << 1.0*(end - start)/CLOCKS_PER_SEC << " seconds." <<  endl;
    
    return 0;
    
}

