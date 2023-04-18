#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.3.0"
__maintainer__ = "Jose Rodriguez"
__email__ = "andrea.laguillo@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# import modules
import os
import sys
import argparse
import configparser
import glob
import logging
import pandas as pd
import numpy as np
import concurrent.futures
from tqdm import tqdm
from itertools import repeat
pd.options.mode.chained_assignment = None  # default='warn'

#infile = r"C:\Users\Andrea\Desktop\SHIFTS-4\testing\cXcorr_Len_Rank_Results_TargetData_Calibration.txt"

# TODO if empty space in one column, ignore whole row (all modules)

###################
# Local functions #
###################
def concatInfiles(infile):
    '''    
    Concat input files...
    '''
    
    # read input file
    # df = pd.read_csv(infile, sep="\t", float_precision='high', low_memory=False)
    df = pd.read_feather(infile)
    #df['Experiment'] = infile[0]
    df['Filename'] = infile
    # add folder name into column
    #foldername = os.path.dirname(infile[0])
    #df['Experiment'] = foldername
    # add filename column
    #df['Filename'] = os.path.basename(infile)
    # assign type to categorical columns
    #df['Experiment'] = df['Experiment'].astype('category')
    df['Filename'] = df['Filename'].astype('category')
    return df

def generate_histogram(df, bin_width):
    '''
    Group by DeltaMass into bins of the size specified.
    '''
    
    def _decimal_places(x):
        s = str(x)
        if not '.' in s:
            return 0
        return len(s) - s.index('.') - 1
    
    # sort by deltamass
    df.sort_values(by=['cal_dm_mh'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    # make bins
    bins = list(np.arange(int(round(df['cal_dm_mh'][0])),
                          int(round(df['cal_dm_mh'].iloc[-1]))+bin_width,
                          bin_width))
    bins = [round(x, _decimal_places(bin_width)) for x in bins]
    df['bin'] = pd.cut(df['cal_dm_mh'], bins=bins)
    
    # make histogram table
    bins_df = df['bin'].value_counts().to_frame().rename(columns = {'bin':'count'})
    bins_df.insert(0, 'bin', bins_df.index)
    bins_df.insert(1, 'midpoint', bins_df['bin'].apply(lambda x: x.mid))
    bins_df.reset_index(drop=True, inplace=True)
    bins_df.sort_values(by=['bin'], inplace=True)
    bins_df.reset_index(drop=True, inplace=True)
    #bins_df['bin'][9].left #access value
    
    return df, bins_df

def _parallelRegression(list_subset, smoothed, second_derivative):
    '''
    Calculate the linear regression line and return the slope
    (and the intercept, if called with smooth == True).
    '''
    bin_subset = list_subset[0]
    pos = list_subset[1]
    result = linear_regression(bin_subset, smoothed, second_derivative)
    if result is tuple:
        result = linear_regression(bin_subset, smoothed, second_derivative)[0]
    result = [pos, result]
    return result

def linear_regression(bin_subset, smoothed, second_derivative):
    '''
    Calculate the linear regression line and return the slope
    (and the intercept, if called with smooth == True).
    '''
    # ignore special cases at beginning and end and use a
    # linear regression function for the rest
    x_list = bin_subset['midpoint'].tolist()
    if smoothed:
        y_list = bin_subset['smooth_count'].tolist()
    elif not second_derivative:
        y_list = bin_subset['count'].tolist()
    else:
        y_list = bin_subset['slope1'].tolist()
    sum1, sum2 = 0, 0
    for i in range(len(x_list)):
        sum1 += (x_list[i] - np.mean(x_list)) * (y_list[i] - np.mean(y_list))
        sum2 += (x_list[i] - np.mean(x_list)) ** 2
    working_slope = sum1 / sum2
    intercept = np.mean(y_list) - working_slope*np.mean(x_list)
    if smoothed or second_derivative:
        return working_slope
    else:
        return working_slope, intercept

def _parallelSmoothing(list_subset):
    bin_subset = list_subset[0]
    pos = list_subset[1]
    working_slope, intercept = linear_regression(bin_subset, False, False)
    #smoothed = pd.Series([pos, working_slope, intercept],
                         #index=['i', 'working_slope', 'intercept'])
    smoothed = [pos, working_slope, intercept]
    return smoothed

def smoothing(bins_df, spoints):
    '''
    Calculate the slope (first derivative) for each bin. Calculate new smoothed
    value for the midpoint using the linear regression line.
    '''
    bins_df['smooth_count'] = None
    # for i in range(spoints, len(bins_df)-spoints):
    #     #working_bin = bins_df.loc[i]
    #     bin_subset = bins_df[i-spoints:i+spoints+1]
    #     working_slope, intercept = linear_regression(bin_subset, False, False)
    #     bins_df.loc[i, 'smooth_count'] = intercept + (working_slope*bins_df.loc[i, 'midpoint'])
    bin_subsets = []
    for i in range(spoints, len(bins_df)-spoints):
        bin_subset = bins_df[i-spoints:i+spoints+1]
        bin_subsets.append([bin_subset, i])
    logging.info("\tSmoothing...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:   
        intercepts = list(tqdm(executor.map(_parallelSmoothing, bin_subsets, chunksize=1000),
                               total=len(bin_subsets)))
    intercepts = pd.DataFrame(intercepts, columns=['i', 'working_slope', 'intercept'])
    # TODO: add to bins_df.loc
    bins_df = pd.merge(bins_df, intercepts, left_index=True, right_on='i', how='outer')
    bins_df.reset_index(drop=True, inplace=True)
    #bins_df['smooth_count'] = bins_df['intercept'] + (bins_df['working_slope']*bins_df['midpoint'])
    bins_df['smooth_count'] = bins_df.apply(lambda x: None if np.isnan(x['intercept'])
                                                           else x['intercept'] + (x['working_slope']*x['midpoint']), axis = 1)
    bins_df[["smooth_count"]] = bins_df[["smooth_count"]].apply(pd.to_numeric)
    bins_df = bins_df.drop(['i', 'working_slope', 'intercept'], axis=1)
    return bins_df

def first_derivative(bins_df, points, spoints):
    '''
    Calculate the slope (first derivative) for each bin.
    Returns the slope of the linear regression line through data points in
    known_y's and known_x's. The slope is the vertical distance divided by
    the horizontal distance between any two points on the line, which is the
    rate of change along the regression line.
    Known_y's  Bins. An array or cell range of numeric dependent data points.
    Known_x's  Count. The set of independent data points.
    '''
    if spoints > 0: #smoothing
        bins_df = smoothing(bins_df, spoints)
        j = 2
    else: #no smoothing
        j = 1
    logging.info("\tFirst derivative...")
    #bins_df['slope1'] = None
    ### TODO: can we make this faster? #######################################
    # for i in range(points*j, len(bins_df)-points*j):
    #     #working_bin = bins_df.loc[i]
    #     bin_subset = bins_df[i-points:i+points+1]
    #     if spoints > 0:
    #         bins_df.loc[i, 'slope1'] = linear_regression(bin_subset, True, False)
    #     else:
    #         bins_df.loc[i, 'slope1'] = linear_regression(bin_subset, False, False)[0] #slope only
    # ##########################################################################
    bin_subsets = []
    for i in range(points*j, len(bins_df)-points*j):
        bin_subset = bins_df[i-points:i+points+1]
        bin_subsets.append([bin_subset, i])
    
    if spoints > 0:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:   
            firsts = list(tqdm(executor.map(_parallelRegression,
                                            bin_subsets, repeat(True), repeat(False),
                                            chunksize=1000), total=len(bin_subsets)))
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:   
            firsts = list(tqdm(executor.map(_parallelRegression,
                                            bin_subsets, repeat(False), repeat(False),
                                            chunksize=1000), total=len(bin_subsets)))
    firsts = pd.DataFrame(firsts, columns=['i', 'slope1'])
    bins_df = pd.merge(bins_df, firsts, left_index=True, right_on='i', how='outer').replace({np.nan: None})
    bins_df = bins_df.drop('i', axis=1)
    bins_df.reset_index(drop=True, inplace=True)
    #bins_df[['slope1']] = bins_df[['slope1']].apply(pd.to_numeric)
    return bins_df

def second_derivative(bins_df, points, spoints):
    '''
    Calculate the second derivative for each bin.
    '''
    if spoints > 0: #smoothed
        j = 3
    else: #not smoothed
        j = 2
    logging.info("\tSecond derivative...")
    #bins_df['slope2'] = None
    bin_subsets = []
    for i in range(points*j, len(bins_df)-points*j):
        bin_subset = bins_df[i-points:i+points+1]
        bin_subsets.append([bin_subset, i])
        #bins_df.loc[i, 'slope2'] = linear_regression(bin_subset, False, True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:
        seconds = list(tqdm(executor.map(_parallelRegression,
                                        bin_subsets, repeat(False), repeat(True),
                                        chunksize=1000), total=len(bin_subsets)))
    seconds = pd.DataFrame(seconds, columns=['i', 'slope2'])
    bins_df = pd.merge(bins_df, seconds, left_index=True, right_on='i', how='outer').replace({np.nan: None})
    bins_df = bins_df.drop('i', axis=1)
    bins_df.reset_index(drop=True, inplace=True)
    #bins_df[['slope2']] = bins_df[['slope2']].apply(pd.to_numeric)           
    return bins_df
#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    
    #Main variables
    bins = float(config._sections['PeakModeller']['bins'])
    slope_points = int(config._sections['PeakModeller']['slope_points'])
    smooth_points = int(config._sections['PeakModeller']['smooth_points'])
    
    logging.info("Reading input file list...")
    if '*' in args.infile: # wildcard
        infiles = glob.glob(args.infile)
        h_outfile = os.path.join(os.path.dirname(args.infile), 'PeakModeller_DMHistogram.tsv')
        t_outfile = os.path.join(os.path.dirname(args.infile), 'PeakModeller_DMTable.feather')
    else:
        with open(args.infile) as f:
            infiles = f.readlines()
        infiles = [x.strip() for x in infiles] # remove whitespace
        infiles = list(filter(None, infiles)) # remove empty lines
        h_outfile = args.infile[:-4] + '_DMHistogram.tsv'
        t_outfile = args.infile[:-4] + '_DMTable.feather'
    for i in infiles:
        logging.info('\t' + str(os.path.basename(i)))
    
    logging.info("Concat input files...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
        df = executor.map(concatInfiles, infiles)
    df = pd.concat(df)
    df.reset_index(drop=True, inplace=True)

    logging.info("Generating DMHistogram...")
    # make bins
    df, bins_df = generate_histogram(df, bins)
    # calculate derivatives
    #grouped_bins_df = bins_df.groupby(['bin'])
    bins_df = first_derivative(bins_df, #does 1st smoothing pass and 2nd normal pass
                               slope_points//2,
                               smooth_points//2)  
    bins_df = second_derivative(bins_df,
                                slope_points//2,
                                smooth_points//2)
        # check which bins pass
    logging.info("Writing output files...")
    # write DMhistogram
    bins_df.to_csv(h_outfile, index=False, sep='\t', encoding='utf-8')
    # write DMtable (input for PeakSelector)
    df = df.astype(str)
    # df.to_csv(t_outfile, index=False, sep='\t', encoding='utf-8')
    df.to_feather(t_outfile)
    logging.info("Peak Modelling finished")

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Modeller',
        epilog='''
        Example:
            python PeakModeller.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Input file with the peak file(s) to be filtered')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    # TODO: output file path

    parser.add_argument('-b', '--bins', help='Width of the bins')
    parser.add_argument('-p', '--slope_points', help='Number of points (bins) to use for slope calculation')
    parser.add_argument('-s', '--smooth_points', help='Number of points (bins) to use for smoothing')

    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.bins is not None:
        config.set('PeakModeller', 'bins', str(args.bins))
        config.set('Logging', 'create_ini', '1')
    if args.slope_points is not None:
        config.set('PeakModeller', 'slope_points', str(args.slope_points))
        config.set('Logging', 'create_ini', '1')
    if args.smooth_points is not None:
        config.set('PeakModeller', 'smooth_points', str(args.smooth_points))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    if '*' in args.infile: # wildcard
        log_file = os.path.join(os.path.dirname(args.infile), 'PeakModeller_log.txt')
        log_file_debug = os.path.join(os.path.dirname(args.infile), 'PeakModeller_log_debug.txt')
    else:
        log_file = args.infile[:-4] + '_log.txt'
        log_file_debug = args.infile[:-4] + '_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')
