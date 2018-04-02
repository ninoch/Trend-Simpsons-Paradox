import numpy as np
import pandas as pd 
from scipy import stats
import statsmodels.api as sm
import bisect
import copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import pickle
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import chi2
import matplotlib.colors as clrs
import json

df = pd.DataFrame()
possible_values_df = dict()
target_variable = None
markers = ['o','v','+','s','p','x','h','d','o','v','+','s','p','x','h','d']

######################################################
##                 INPUT FUNCTIONS                  ##
######################################################


def read_the_csv_file(file_name, tar_var, drop_variables = None, filter_variables = None):
    global df
    df = pd.read_csv(file_name)
    if drop_variables:
        df = df.drop(drop_variables, axis=1)
    if filter_variables: 
        df = df.loc[filter_variables]
    name_of_variables = df.columns.values
    
    var = list(name_of_variables)
    target_variable = tar_var
    var.remove(target_variable)

    return var


def read_the_bining_file(file_name):
    global possible_values_df
    last_var = ""
    for line in open(file_name):
        csv_row = line.split(',')
        last_var = csv_row[0]
        possible_values_df[csv_row[0]] = csv_row[1:]
        possible_values_df[csv_row[0]][-1] = possible_values_df[csv_row[0]][-1][:-1]
        for ind in range(len(possible_values_df[csv_row[0]])):
            possible_values_df[csv_row[0]][ind] = float(possible_values_df[csv_row[0]][ind])

    possible_values_df[csv_row[0][:-1]] = possible_values_df[csv_row[0]]


def different_vals_of(var):
    if var in possible_values_df:
        arr = possible_values_df[var]
        arr[0] = arr[0] - 1
        return list(arr)
    return None


def read_input_info():
    with open("../input_info.json", 'r') as f:
        input_json_info = json.load(f)
    target_variable = input_json_info["target_variable"]
    level_of_significance = input_json_info["level_of_significance"]
    
    name_of_the_variables = read_the_csv_file("../input/" + input_json_info["csv_file_name"] + ".csv", target_variable, input_json_info["ignore_columns"])
    read_the_bining_file("../temporary_files/" + input_json_info["csv_file_name"] + "/bins.csv")
    log_scales = input_json_info["log_scales"]

    return target_variable, level_of_significance, name_of_the_variables, log_scales

#####################################################
##                 PAIR FUNCTIONS                  ##
#####################################################


def set_paradox_conditioning_pairs(name_of_variables, paradox_vars = None, conditioning_vars = None, my_pairs = None):
    pairs = []

    if not paradox_vars:
        paradox_vars = copy.copy(name_of_variables)
    if not conditioning_vars:
        conditioning_vars = copy.copy(name_of_variables)

    if my_pairs:
        for par_var, cond_var in my_pairs:
            if par_var != cond_var:
                pairs.append((par_var, cond_var))
    else:
        for par_var in paradox_vars:
            for cond_var in conditioning_vars:
                if par_var != cond_var:
                    pairs.append((par_var, cond_var))
    return pairs


########################################################
##                 DRAWING FUNCTIONS                  ##
########################################################


def compute_mean_std(var, lvar, rvar, cond = -1, lcond = -1, rcond = -1, confidence = 0.95):
    a = []
    if cond == -1:
        a = df.loc[(df[var] > lvar) & (df[var] <= rvar)][target_variable].values
    else:
        a = df.loc[(df[cond] > lcond) & (df[cond] <= rcond)].loc[(df[var] > lvar) & (df[var] <= rvar)][target_variable].values
    n = len(a)
    se = stats.sem(a)
    h = se * stats.t._ppf((1+confidence)/2., n-1)
    return np.mean(a), h


def draw(trend_simpsons_pair, aggregated_vars_params, disaggregated_vars_params, log_scales, max_group = 10):
    for var, cond in trend_simpsons_pair:
    	# Making first page
    	print "Making", str(var + '-vs-' + cond + '.pdf'), " file "
    	if disaggregated_vars_params[var + "," + cond]["params"] == []:
    		continue
        pp = PdfPages("../output/" + var + '-vs-' + cond + '.pdf') #TODO
        plt.figure() 
        plt.axis('off')
        text = var + " with conditioning on " + cond 
        plt.text(0.5, 0.5, text, ha='center', va='center')
        pp.savefig(bbox_inches='tight', papertype='a4')
        plt.close() 

        # Drawing group plots
        print "\t Drawing chart aggregated"
        plt.clf()
        coef, inter = aggregated_vars_params[var]["params"][1], aggregated_vars_params[var]["params"][0]
            
        possible_values = different_vals_of(var)
        y_actual, y_err, y_hat = [], [], []
        for indx in range(len(possible_values) - 1):
            m, e = compute_mean_std(var, possible_values[indx], possible_values[indx + 1])
            y_actual.append(m)
            y_err.append(e)
        
        x_hat = np.arange((possible_values[0] + possible_values[1]) / 2, (possible_values[-2] + possible_values[-1]) / 2, float(possible_values[-1] - possible_values[0]) / 100.0)
        y_hat = 1 / (1 + np.exp(-(coef * x_hat + inter)))
        plt.plot(x_hat, y_hat, linewidth=1, linestyle='dashed', color='k',label='logistic fit')
        plt.errorbar(np.array(np.array(possible_values[1:]) + np.array(possible_values[:-1])) / 2, y_actual, yerr=[y_err, y_err], alpha=0.75, color='black', label='data', fmt='o')
        plt.xlabel(var)
        plt.ylabel(target_variable)

        if log_scales[var]:
            plt.xscale('log')

        plt.legend(loc='best')
        plt.title(var)
        pp.savefig(bbox_inches='tight', papertype='a4')
        plt.close()

        print "\t Drawing chart disaggregated"
        disaggregated_vars_params[var + "," + cond]["params"] = np.array(disaggregated_vars_params[var + "," + cond]["params"])
        coefs = disaggregated_vars_params[var + "," + cond]["params"][:, 1]
        inters = disaggregated_vars_params[var + "," + cond]["params"][:, 0]
        
        plt.clf()
        conditioning_groups = different_vals_of(cond)
        coefs_ind = 0
        for ind in range(len(conditioning_groups) - 1):
            X = df.loc[(df[cond] > conditioning_groups[ind]) & (df[cond] <= conditioning_groups[ind + 1])][var].values
            if len(X) < 100:
                continue
            if ind >= max_group:
                break
            X_lables = possible_values[(bisect.bisect_left(possible_values, X.min()) - 1):(bisect.bisect_right(possible_values, X.max()) + 1)]
            y_actual, y_err = [], []
            for indx in range(len(X_lables)-1): 
                m, e = compute_mean_std(var, X_lables[indx], X_lables[indx + 1], cond, conditioning_groups[ind], conditioning_groups[ind + 1])
                y_actual.append(m)
                y_err.append(e)

            colorr = float(float(ind + 1)/float(min(max_group, len(conditioning_groups))))
            plt.errorbar(np.array(np.array(X_lables[1:]) + np.array(X_lables[:-1])) / 2, y_actual, yerr=[y_err, y_err], color=(colorr,0,1.0 - colorr), alpha=0.75, fmt=markers[ind % len(markers)], label= '(' + str(conditioning_groups[ind]) + "-" + str(conditioning_groups[ind + 1]) + ']' )
        	
            if len(X_lables) > 1:
                x_hat = np.arange((X_lables[0] + X_lables[1]) / 2, (X_lables[-2] + X_lables[-1]) / 2, float(X_lables[-1] - X_lables[0]) / 100.0)
                y_hat = 1 / (1 + np.exp(-(coefs[coefs_ind] * x_hat + inters[coefs_ind])))
                plt.plot(x_hat, y_hat, color=(colorr,0,1.0 - colorr), linewidth=1, linestyle='dashed')
            coefs_ind += 1

        plt.xlabel(var)
        if log_scales[var]:
            plt.xscale('log')
        
        plt.ylabel(target_variable)
        try:
        	plt.legend(loc='best')
        except:
        	pass
        plt.title(var + " with conditioning on " + cond)
        pp.savefig(bbox_inches='tight', papertype='a4')
        plt.close()

        # Drawing Pcolormesh Plot 
        mat_mean = np.zeros((len(conditioning_groups) - 1, len(possible_values) - 1))
        mat_freq = np.zeros((len(conditioning_groups) - 1, len(possible_values) - 1))

        
        print "\t Drawing pcolormesh plots"
        for x in range(len(possible_values) - 1):
            for y in range(len(conditioning_groups) - 1):
                lvar, rvar = possible_values[x], possible_values[x + 1]
                lcond, rcond = conditioning_groups[y], conditioning_groups[y + 1]
                target_array = df.loc[(df[cond] > lcond) & (df[cond] <= rcond)].loc[(df[var] > lvar) & (df[var] <= rvar)][target_variable].values
                mat_mean[y][x] = np.mean(target_array)
                mat_freq[y][x] = len(target_array) 

        mat_mean = np.ma.masked_where(np.isnan(mat_mean),mat_mean)

        plt.clf()
        fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True)#, gridspec_kw={'width_ratios':[1,1,1.25]})

        cmap1 = cm.bwr
        cmap1.set_bad('lightgray', 1.)
        im = ax[0].pcolormesh(np.array([possible_values[1:]] * (len(conditioning_groups) - 1)), np.array([conditioning_groups[1:]] * (len(possible_values) - 1)).T, mat_mean,
	                            vmin=mat_mean.min() - 0.01,
	                            vmax=mat_mean.max() + 0.01,
	                            cmap = cmap1)

        fig.colorbar(im, ax=ax[0])
        if log_scales[var]:
            ax[0].set_xscale('log')
        if log_scales[cond]:
            ax[0].set_yscale('log')

        cmap2 = cm.YlGn
        cmap2.set_bad('lightgray', 1.)
        im = ax[1].pcolormesh(np.array([possible_values[1:]] * (len(conditioning_groups) - 1)), np.array([conditioning_groups[1:]] * (len(possible_values) - 1)).T, mat_freq,
	                            norm=clrs.LogNorm(vmin=1, vmax=mat_freq.max()),
	                            cmap = cmap2)

        fig.colorbar(im, ax=ax[1])
        if log_scales[var]:
            ax[0].set_xscale('log')
            ax[1].set_xscale('log')
        if log_scales[cond]:
            ax[1].set_yscale('log')

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.grid(False)
        plt.xlabel(var)
        plt.ylabel(cond)
        pp.savefig(bbox_inches='tight', papertype='a4')
        pp.close()
        plt.close()



########################################################
##                 STORING FUNCTIONS                  ##
########################################################

def store_info(file_name, obj):
    with open("../store_results/" + file_name, 'wb') as handle: #TODO
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)


def store_all_info(trend_simpsons_pairs, aggregated_vars_params, disaggregated_vars_params):
    store_info("simpsons_pairs.obj", trend_simpsons_pairs)
    store_info("aggregated_vars_params.obj", aggregated_vars_params)
    store_info("disaggregated_vars_params.obj", disaggregated_vars_params)

def load_info(file_name):
    with open("../store_results/" + file_name, 'rb') as handle: #TODO
        obj = pickle.load(handle)
    return obj

def load_all_info():
    trend_simpsons_pairs = load_info("simpsons_pairs.obj")
    aggregated_vars_params = load_info("aggregated_vars_params.obj")
    disaggregated_vars_params = load_info("disaggregated_vars_params.obj")
    paradox_ranking = load_info("paradox_ranking.obj")
    return trend_simpsons_pairs, aggregated_vars_params, disaggregated_vars_params, paradox_ranking

################################################################
##                 TREND SIMPSON'S FUNCTIONS                  ##
################################################################

def sign(val, pval = 0):
    if pval >= level_of_significance or val == 0:
        return 0
    elif val < 0:
        return -1
    return 1

def compute_mean_of_error(lr, Y, X):
    delta = lr.predict(X)
    inds = np.array(Y == 0)
    delta[inds] = 1 - delta[inds]
    return 1.0 - float(sum(delta) / float(len(delta)))

def logistic_regression(X, Y):
    if max(Y) == min(Y):
        if max(Y) == 1:
            return [20, 0], [0, 1], 0 # 
        else:
            return [-20, 0], [0, 1], 0 # 

    logistic_model = sm.Logit(Y,X)
    try:
        lr = logistic_model.fit(maxiter=200, disp=False) 
    except: # Singular Matrix
        intercept = -1 * np.log((1 / np.mean(Y)) - 1)
        return [intercept, 0], [0, 1], 0 #  

    #print lr.summary()
    pvals = lr.pvalues
    params = lr.params

    return lr.params, lr.pvalues, compute_mean_of_error(lr, Y, X)


def find_trend_simpsons_pairs(pairs):
    trend_simpsons_pairs = []
    aggregated_vars_params = {}
    disaggregated_vars_params = {}

    for var, cond in pairs:
        agg_sign = 0
        disagg_sign = 0

        print ""
        print "*** Running Logistic Regression on", var, "***"
        if var in aggregated_vars_params:
            agg_params, agg_pvalues, agg_error = aggregated_vars_params[var]["params"], aggregated_vars_params[var]["pvalues"], aggregated_vars_params[var]["error"]
        else:
            X = df[var].values
            X = np.array([X]).T
            X = sm.tools.tools.add_constant(X, has_constant='add')
            Y = df[target_variable].values
            agg_params, agg_pvalues, agg_error = logistic_regression(X, Y)
            aggregated_vars_params[var] = {}
            aggregated_vars_params[var]["params"], aggregated_vars_params[var]["pvalues"], aggregated_vars_params[var]["error"] = agg_params, agg_pvalues, agg_error

        agg_sign = sign(agg_params[1], agg_pvalues[1])

        #print "Coefficient: ", agg_params[1], "(", agg_pvalues[1], "), ", " Intercept: ", agg_params[0], "(", agg_pvalues[0] ,")", " Mean of errors: ", agg_error
        
        print ""
        print "*** Running Logistic Regression on", var, "conditioning", cond, "***"
        disagg_params = []
        disagg_pvalues = []
        disagg_errors = []

        conditioning_groups = different_vals_of(cond)
        for ind in range(len(conditioning_groups) - 1):
            the_df = df.loc[(df[cond] > conditioning_groups[ind]) & (df[cond] <= conditioning_groups[ind + 1])]
            if len(the_df) < 100:
                print ">>>>>>>>>>> BINNING :", "len(df[(df[", cond ,"] > ", conditioning_groups[ind], ") & (df[", cond, "] <= ", conditioning_groups[ind + 1], ")]) = ", len(the_df)
                continue
            X = the_df[var].values
            X = np.array([X]).T
            X = sm.tools.tools.add_constant(X, has_constant='add')
            Y = the_df[target_variable].values
            pars, pval, err = logistic_regression(X, Y)

            disagg_params.append(pars)
            disagg_pvalues.append(pval)
            disagg_errors.append(err)
            
            #print "For group", (ind + 1), " : Coefficient: ", pars[1], "(", pval[1], "), ", " Intercept: ", pars[0], "(", pval[0] ,")", " Mean of errors: ", err

            disagg_sign += sign(pars[1], pval[1])
        
        disagg_sign = sign(float(disagg_sign) / float(len(disagg_params)))
        disaggregated_vars_params[var + "," + cond] = {}
        disaggregated_vars_params[var + "," + cond]["params"] = disagg_params
        disaggregated_vars_params[var + "," + cond]["pvalues"] = disagg_pvalues 
        disaggregated_vars_params[var + "," + cond]["errors"] = disagg_errors

        if agg_sign != 0 and (disagg_sign != agg_sign):
            print ""
            print ">>>>> Trend Simpson's instance for var", var, "with conditioning on", cond, "SIGNS: ", agg_sign, disagg_sign, "<<<<<"
            print ""
            trend_simpsons_pairs.append([var, cond])

    return trend_simpsons_pairs, aggregated_vars_params, disaggregated_vars_params

def show_deviance_ranking(pairs, deviance_ranking):
    print "######################################################"
    print "## Deviance improvement ranking for finalized pairs ##"
    print "######################################################"

    mvar = np.unique([i[0] for i in pairs])
    for var in mvar: 
        for key, dr in deviance_ranking[::-1]:
            if key.startswith(var + ","):
                print key, float("{0:.2f}".format(dr))
        print "------------------------------"


def ranking_deviance(finalized_pairs):
    tmp = load_info("loglikelihoods.obj")
    
    null_agg = tmp['null_agg']
    full_agg = tmp['full_agg']
    null_disagg = tmp['null_disagg']
    full_disagg = tmp['full_disagg']
    rank = dict()

    for var, cond in finalized_pairs:
        key = var + "," + cond
        #n_a = null_agg[var]
        f_a = full_agg[var] # remove
        f_da = full_disagg[key]

        rank[key] = 1 - float(float(f_da) / float(f_a))
        
    return sorted(rank.iteritems(), key=lambda (k,v): (v,k))


def chi_sq_deviance():
    print "################################################################################"
    print "## Applying Chi-squared Deviance Test to all pairs and finding finalized ones ##"
    print "################################################################################"

    tmp = load_info("loglikelihoods.obj")
    
    null_agg = tmp['null_agg']
    full_agg = tmp['full_agg']
    null_disagg = tmp['null_disagg']
    full_disagg = tmp['full_disagg']
    finalized_pairs = []

    for key, n_da in null_disagg.iteritems():
        f_da = full_disagg[key]
        var, cond = key.split(',')[0], key.split(',')[1]
        n_a = null_agg[var]
        f_a = full_agg[var]

        chi_agg = 2 * (f_a - n_a)
        df_agg = 1
        chi_dagg = 2 * (f_da - n_da)
        df_dagg = len(different_vals_of(cond))
        print "Chi-squared test for pair", key, ": agg = ", (1 - chi2.cdf(chi_agg, df_agg)), " disagg = ", (1 - chi2.cdf(chi_dagg, df_dagg))

        if f_a < f_da and (1 - chi2.cdf(chi_agg, df_agg)) < level_of_significance and (1 - chi2.cdf(chi_dagg, df_dagg)) < level_of_significance:
            print "\t Pass"
            finalized_pairs.append((var, cond))
        else:
            print "\t Not pass"
    return finalized_pairs


def deviance_all_pairs(pairs, aggregated_vars_params, disaggregated_vars_params):
    print "##############################################################################################"
    print "## Computing Loglikelihood for full / null, aggregated / disaggregated models for all pairs ##"
    print "##############################################################################################"

    theta_0 = np.mean(df[target_variable].values)

    null_agg = dict()
    full_agg = dict()
    null_disagg = dict()
    full_disagg = dict()

    for var, cond in pairs:
    	print "Computing Deviance for pair [", var, ", ", cond, "]"

    	# Aggregated Deviance 
        if var not in null_agg:
    	    coef, inter = aggregated_vars_params[var]["params"][1], aggregated_vars_params[var]["params"][0]

            null_agg_tmp = 0
            full_agg_tmp = 0
            for val, y_i in zip(df[var].values, df[target_variable].values):
                y_hat_i = np.sum(1 / (1 + np.exp(-(coef * val + inter))))
                if y_i: 
                    null_agg_tmp += y_i * np.log(theta_0)
                    full_agg_tmp += y_i * np.log(y_hat_i) 
                else:
                    null_agg_tmp += (1 - y_i) * np.log(1 - theta_0)
                    full_agg_tmp += (1 - y_i) * np.log(1 - y_hat_i)

            null_agg[var] = null_agg_tmp
            full_agg[var] = full_agg_tmp


        # Disaggregated Deviance
        disaggregated_vars_params[var + "," + cond]["params"] = np.array(disaggregated_vars_params[var + "," + cond]["params"])
        coefdisagg = disaggregated_vars_params[var + "," + cond]["params"][:, 1]
        interdisagg = disaggregated_vars_params[var + "," + cond]["params"][:, 0]

        null_disagg_tmp = 0
        full_disagg_tmp = 0
        disagg_ind = 0

        conditioning_groups = different_vals_of(cond)
        for ind in range(len(conditioning_groups) - 1):
            the_df = df.loc[(df[cond] > conditioning_groups[ind]) & (df[cond] <= conditioning_groups[ind + 1])]
            theta_1 = np.mean(the_df[target_variable].values)
            if len(the_df) < 100:
                continue
            for val, y_i in zip(the_df[var].values, the_df[target_variable].values):
                y_hat_i = np.sum(1 / (1 + np.exp(-(coefdisagg[disagg_ind] * val + interdisagg[disagg_ind]))))
                if y_i:
                    null_disagg_tmp += y_i * np.log(theta_1) 
                    full_disagg_tmp += y_i * np.log(y_hat_i)
                else:
                    null_disagg_tmp += (1 - y_i) * np.log(1 - theta_1)
                    full_disagg_tmp += (1 - y_i) * np.log(1 - y_hat_i)


            disagg_ind += 1

        null_disagg[var + "," + cond] = null_disagg_tmp
        full_disagg[var + "," + cond] = full_disagg_tmp

        print "\t aggregated( null / full ):", null_agg[var], full_agg[var], " disaggregated( null / full ):", null_disagg[var + "," + cond], full_disagg[var + "," + cond]

    store_info("loglikelihoods.obj", {'null_agg': null_agg, 'full_agg': full_agg, 'null_disagg': null_disagg, 'full_disagg': full_disagg})


if __name__ == "__main__":
    # Reading info from input_info.txt
    target_variable, level_of_significance, name_of_the_variables, log_scales = read_input_info()

    # Finding all Trend Simpson's Pairs 
    pairs = set_paradox_conditioning_pairs(name_of_the_variables)
    trend_simpsons_pairs, aggregated_vars_params, disaggregated_vars_params = find_trend_simpsons_pairs(pairs)
    store_all_info(trend_simpsons_pairs, aggregated_vars_params, disaggregated_vars_params)
    # To prevent computing trend simpson's pairs again, you can load results from your previous run
    #trend_simpsons_pairs, aggregated_vars_params, disaggregated_vars_params, paradox_ranking = load_all_info() 
    print "Number of all pairs found as Simpson's: ", len(trend_simpsons_pairs)
    print trend_simpsons_pairs

    # Finding finalized pairs based on Deviace chi-squared measure
    deviance_all_pairs(trend_simpsons_pairs, aggregated_vars_params, disaggregated_vars_params)
    finalized_pairs = chi_sq_deviance()
    print "Number of all finalized pairs: ", len(finalized_pairs)
    print finalized_pairs

    # Ranking for pairs 
    deviance_ranking = ranking_deviance(finalized_pairs)
    show_deviance_ranking(finalized_pairs, deviance_ranking)
    
    # Drawing charts for finalized pairs
    draw(finalized_pairs, aggregated_vars_params, disaggregated_vars_params, log_scales)

