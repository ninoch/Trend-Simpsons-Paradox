#ifndef main_h
#define main_h

void read_single_column(int col, ifstream & datafile, double & val)
{
    // skip the the first col-1 entries
    for(int i=0; i<col; i++)
        datafile.ignore(numeric_limits<streamsize>::max(),',');
    
    // take the next entry
    datafile >> val;
    
    // skip the start of the next line
    datafile.ignore(numeric_limits<streamsize>::max(),'\n');
    
}


void read_x_y_from_column(ifstream & datafile, int xcol, double & x, int ycol, double & y)
{
    if(xcol < ycol)
    {
        // skip the the first col-1 entries
        for(int i=0; i<xcol; i++)
            datafile.ignore(numeric_limits<streamsize>::max(),',');
        
        // take the next entry and ignore the following comma
        datafile >> x;
        datafile.ignore(1,',');
        
        // skip the next remaining entries to y
        for(int i=xcol+1; i<ycol; i++)
            datafile.ignore(numeric_limits<streamsize>::max(),',');
        
        // take the next entry
        datafile >> y;
        
        // skip the start of the next line
        datafile.ignore(numeric_limits<streamsize>::max(),'\n');
    }
    else
    {
        // skip the the first col-1 entries
        for(int i=0; i<ycol; i++)
            datafile.ignore(numeric_limits<streamsize>::max(),',');
        
        // take the next entry and ignore the following comma
        datafile >> y;
        datafile.ignore(1,',');
        
        // skip the next remaining entries to y (+1 as we are at the next column after ycol)
        for(int i=ycol+1; i<xcol; i++)
            datafile.ignore(numeric_limits<streamsize>::max(),',');
        
        // take the next entry
        datafile >> x;
        
        // skip the start of the next line
        datafile.ignore(numeric_limits<streamsize>::max(),'\n');
    }
    
}

void set_ifstream_to_row1(ifstream & datafile)
{
    // go to beginning of file
    datafile.clear();
    datafile.seekg(0, ios::beg);
    
    // skip first line
    datafile.ignore(numeric_limits<streamsize>::max(),'\n');
}

void create_Xj_info(int & j, vector<map<double,unsigned int> > & cumcount, vector<map<double,double> > & cumsum, set<double> & range, ifstream & datafile, vector<unsigned int> & groups, int & G, int y_col, int N)
{
    // initialize the cumsum and cumcount vectors as empty
    cumcount = vector<map<double,unsigned int> > (G, map<double,unsigned int>());
    cumsum = vector<map<double,double> > (G, map<double,double>());
    range.clear();
    
    // go to start of data
    set_ifstream_to_row1(datafile);
    
    // read in the data
    double Xj;
    double y;
    
    for(int i=0; i<N; i++)
    {
        // read the (y,x) value
        read_x_y_from_column(datafile, j, Xj, y_col, y);
        
        // check if the Xj value already exists in the group
        if(cumcount[groups[i]].find(Xj) == cumcount[groups[i]].end())
        {
            // value doesn;t exist
            cumcount[groups[i]][Xj] = 1;
            cumsum[groups[i]][Xj] = y;
        }
        else
        {
            // value does exist
            cumcount[groups[i]][Xj] ++;
            cumsum[groups[i]][Xj] += y;
        }
    }
    
    
    // All data read in - now make it cumulative for each group and create range of X vector
    for(int g=0; g<G; g++)
    {
        // check to see if the group is greater than one (it is zero, nothing goes, if it is one then cum is done alread)
        if(cumcount[g].size() == 0)
            continue;
        
        // iterate through the groups
        range.insert(cumcount[g].begin()->first); // add the first key to the range as it wont be added after
        for(auto it = std::next(cumcount[g].begin()); it!=cumcount[g].end(); it++)
        {
            if(it->second == 0)
                cout << "ERROR: Have a 0 cumcount element" << endl;
            it->second += (std::prev(it))->second;
            cumsum[g][it->first] += cumsum[g][(std::prev(it))->first];
            range.insert(it->first); // add the key to the range
        }
        
    }
    
}


void calculate_best_split(vector<map<double,unsigned int> > &cumcount, vector<map<double,double> > &cumsum, double Xl, double Xu, set<double> &range, int G, double & best_split, double & best_R2_improvement)
{// calculate the best split s in [Xl,Xu) to split the data and the corresponding R2 improvement
    
    
    // i) Declare R2 improvement for each element x in [Xl,Xu) and initialize to 0
    map<double,double> R2_improvements;
    for(auto it = range.find(Xl); it!=range.find(Xu); ++it)
        R2_improvements[*it] = 0;
    //cout << "Range size = " << R2_improvements.size() << endl;
    
    // i') If there is only one element, cannot split so return
    if(R2_improvements.size() == 0)
        best_R2_improvement = 0;
    else
    {
        
        // ii) loop over every group, for each group finding the contribution of each x in [Xl,Xu) to the R2
        unsigned int countterm_s, countend;
        double ysum_s, ysumend, subtraction_term;
        map<double,unsigned int>::iterator it_end;
        for(int g=0; g<G; g++)
        {
            if(cumcount[g].size() == 0)
                continue;
            
            // 0) if there are no elements in the group:
            // i.e., if all the elements are at \geq Xu or all elements are \leq Xl
            if((cumcount[g].begin()->first >= Xu)||(std::prev(cumcount[g].end())->first <= Xl))
                continue;
            else
            {
                // a) take the larges x <= Xu as this will frequently be needed for the subraction
                it_end = std::prev(cumcount[g].end());
                
                // cout << 1 << " " << flush;
                // group end is greater than Xl, now find largest element in range of group that is not larger than Xu;
                while(it_end->first > Xu)
                    it_end --;
                //cout << 1.1 << " " << flush;
                countend = it_end->second;
                //cout << 2 << " " << flush;
                if(countend == 0) // if there are no datapoints in the group
                    continue;
                ysumend = cumsum[g][it_end->first];
                subtraction_term = ysumend*ysumend/countend;
                //cout << 3 << " " << flush;
            }
            
            // b) iterate through the elements and give contribution to R2 (dont count the last element as it cannot be split there (last element is in group at least))
            for( auto it = cumcount[g].begin(); it != std::prev(cumcount[g].end()); ++it)
            {
                // only interested in Xj values in the range [Xl,Xu)
                if(it->first < Xl)
                    continue;
                if(it->first >= Xu)
                    break;
                
                //s = it->first;
                countterm_s = it->second;
                if(countterm_s == 0)
                    cout << "WARNINIG: countterm_s=0; " << it->first << " in group " << g << " (lims [" << cumcount[g].begin()->first << "," << std::prev(cumcount[g].end())->first << ")) in [" << Xl << "," << Xu << ")" << endl;
                ysum_s = cumsum[g][it->first];
                
                R2_improvements[it->first] += (ysum_s*ysum_s/countterm_s) + (ysumend-ysum_s)*(ysumend-ysum_s)/(countend-countterm_s) - subtraction_term;
                
            }
            
        }
        //cout << 4 << " " << endl;
        
        // Have calculated R2_improvements for each split, now return the best
        best_R2_improvement = 0;
        for(auto it = R2_improvements.begin(); it != R2_improvements.end(); it++)
        {
            // cout << it->first << "\t" << it->second << endl;
            if(it->second > best_R2_improvement)
            {
                best_R2_improvement = it->second;
                best_split = it->first;
            }
            
            
        }
    }
    
}

void update_R2_DL(double best_R2_improvement, vector<double> & R2_improvements, double R2, vector<double> & DL, double lambda)
{
    // update the R2 for the current variable
    R2_improvements.push_back(best_R2_improvement);
    
    // update the change in loss function for the current variable
    DL.push_back(-(best_R2_improvement/(1-R2)) + lambda);
    
}


void get_column_names(ifstream & datafile, vector<string> & columns)
{
    string s;
    char c;
    
    datafile.get(c);
    // read columns until get to end of the line
    while(true)
    {
        
        // read characters into s until get to a ,
        while((c != ',')&&(c != '\n'))
        {
            s.push_back(c);
            datafile.get(c);
        }
        
        // add string to columns and clear s (and advance to next symbol)
        columns.push_back(s);
        s.clear();
        
        // check to see if reached end of line
        if(c == '\n')
            break;
        else
            datafile.get(c);
        
    }
    
}

bool lambda_stopping_condition(vector<double> & DL, double best_R2_improvement, vector<double> & R2_improvements, vector<double> & splits, double & total_R2_improvement)
{
    int n_consec_increasing = 0;
    for(int i=DL.size()-1; i>=0; i--)
        if(DL[i] > 0)
            n_consec_increasing ++;
        else
            break;
    
    if(n_consec_increasing >= DL_thresh)
    {
        // restructure R2_improvements and splits_vectors
        R2_improvements.resize(R2_improvements.size() - DL_thresh);
        total_R2_improvement = 0;
        for(int i=0; i<R2_improvements.size(); i++)
            total_R2_improvement += R2_improvements[i];
        splits.resize(splits.size() - DL_thresh);
        return true;
    }
    else
    {
        // if there is no R2 improvement then quit (after only chosing elements for which LF was decreasing
        if(best_R2_improvement == 0)
        {
            R2_improvements.resize(R2_improvements.size() - n_consec_increasing);
            total_R2_improvement = 0;
            for(int i=0; i<R2_improvements.size(); i++)
                total_R2_improvement += R2_improvements[i];
            splits.resize(splits.size() - n_consec_increasing);
            return true;
        }
        else
            return false;
    }
    
}

void recalculate_partitioned_cums(vector<map<double,unsigned int> > &cumcount, vector<map<double,double> > &cumsum, double & pl1, double & pl2, double & pu2, int G)
{
    
    unsigned int cumcount1, total_cumcount1 = 0, total_cumcount2 = 0;
    double cumsum1, total_cumsum1 = 0, total_cumsum2 = 0;
    int i=0;
    // for each group, subtract from each element in the second partitions (i.e., [pl2,pu2)) the cumcount and cumsum from the end of the previous partition (<pl2)
    for(int g=0; g<G; g++)
    {
        if(cumcount[g].size() == 0)
            continue;
        
        // The first Xj element of group g has to be below the start of the current group
        if(cumcount[g].begin()->first >= pl2)
            continue;
        
        for(auto it = std::next(cumcount[g].begin()); it!=cumcount[g].end(); ++it) // (start at second element as first is redundant (no previous elements))
        {
            // iterate through cumcount[g] until arrive at pl2. If that doesn't happen (i.e., if end is first), grand
            if(it->first >= pl2)
            {
                // previous Xj value should also be in the previous interval
                if(std::prev(it)->first >= pl1)
                {
                    cumcount1 = std::prev(it)->second;
                    total_cumcount1 += cumcount1;
                    if(cumcount1 == 0)
                        cout << "WARNING 2: cumcount1=0"<<endl;
                    cumsum1 = cumsum[g][std::prev(it)->first];
                    total_cumsum1 += cumsum1;
                    //cout << "Last element of previous interval = " << std::prev(it)->first << ", cumcount1 = " << cumcount1 << ", last element of interval = " << pu2 <<  endl;
                    i=0;
                    while((it->first <= pu2)&&(it != cumcount[g].end()))
                    {
                        //cout << it->first << " - Before: " << it->second;
                        it->second -= cumcount1;
                        //cout << "\tAfter: " << it->second << endl;;
                        if(it->second == 0)
                        {
                            cout << "FIRST WARNING: cumcount[g]=0 at position "<< i << " in interval" << endl;
                            cout << "Prev Xj value is " << std::prev(it)->first << ", current value is " << it->first << endl;
                        }
                        cumsum[g][it->first] -= cumsum1;
                        it++;
                        //cout << "IT distance: " << std::distance(cumcount[g].end(), it) << endl;
                        //if(it == cumcount[g].end())
                        //break;
                    }
                    
                }
                
                // all done, so break to next group
                break;
                
            }
            
        }
        
        total_cumcount2 += std::prev(cumcount[g].end())->second;
        total_cumsum2 += std::prev(cumsum[g].end())->second;
        
    }
    
    // cout << "\nN1 = " << total_cumcount1 << ", N2 = " << total_cumcount2 << ", sum = " << total_cumcount1+total_cumcount2 << endl;
    // cout << "ybar1 = " << total_cumsum1/total_cumcount1 << ", ybar2 = " << total_cumsum2/total_cumcount2 << endl;
    
    
}

#endif
