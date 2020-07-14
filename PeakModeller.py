#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "andrea.laguillo@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# import modules
import os
import sys
import argparse
import configparser
import logging
import pandas as pd
import numpy as np
import concurrent.futures
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
    df = pd.read_csv(infile, sep="\t", float_precision='high')
    # add folder name into column
    foldername = os.path.dirname(infile)
    df['Experiment'] = foldername
    # add filename column
    df['Filename'] = os.path.basename(infile)
    # assign type to categorical columns
    df['Experiment'] = df['Experiment'].astype('category')
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

def smoothing(bins_df, spoints):
    '''
    Calculate the slope (first derivative) for each bin. Calculate new smoothed
    value for the midpoint using the linear regression line.
    '''
    bins_df['smooth_count'] = None
    for i in range(spoints, len(bins_df)-spoints):
        #working_bin = bins_df.loc[i]
        bin_subset = bins_df[i-spoints:i+spoints+1]
        working_slope, intercept = linear_regression(bin_subset, False, False)
        bins_df.loc[i, 'smooth_count'] = intercept + (working_slope*bins_df.loc[i, 'midpoint'])
    bins_df[["smooth_count"]] = bins_df[["smooth_count"]].apply(pd.to_numeric)
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
    bins_df['slope1'] = None
    for i in range(points*j, len(bins_df)-points*j):
        #working_bin = bins_df.loc[i]
        bin_subset = bins_df[i-points:i+points+1]
        if spoints > 0:
            bins_df.loc[i, 'slope1'] = linear_regression(bin_subset, True, False)
        else:
            bins_df.loc[i, 'slope1'] = linear_regression(bin_subset, False, False)[0] #slope only
    bins_df[['slope1']] = bins_df[['slope1']].apply(pd.to_numeric)
    return bins_df

def second_derivative(bins_df, points, spoints):
    '''
    Calculate the second derivative for each bin.
    '''
    if spoints > 0: #smoothed
        j = 3
    else: #not smoothed
        j = 2
    bins_df['slope2'] = None
    for i in range(points*j, len(bins_df)-points*j):
        bin_subset = bins_df[i-points:i+points+1]
        bins_df.loc[i, 'slope2'] = linear_regression(bin_subset, False, True)
    bins_df[['slope2']] = bins_df[['slope2']].apply(pd.to_numeric)           
    return bins_df

#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    
    logging.info("Reading input file list...")
    with open(args.infile) as f:
        infiles = f.readlines()
    infiles = [x.strip() for x in infiles] # remove whitespace
    
    logging.info("Concat input files...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
        df = executor.map(concatInfiles, infiles)
    df = pd.concat(df)
    df.reset_index(drop=True, inplace=True)

    logging.info("Generating DMHistogram...")
    # make bins
    df, bins_df = generate_histogram(df, float(config._sections['PeakModeller']['bins']))
    # calculate derivatives
    #grouped_bins_df = bins_df.groupby(['bin'])
    bins_df = first_derivative(bins_df, #does 1st smoothing pass and 2nd normal pass
                               int(config._sections['PeakModeller']['points'])//2,
                               int(config._sections['PeakModeller']['spoints'])//2)  
    bins_df = second_derivative(bins_df,
                                int(config._sections['PeakModeller']['points'])//2,
                                int(config._sections['PeakModeller']['spoints'])//2)
        # check which bins pass
    logging.info("Writing output files...")
    # write DMhistogram
    outfile = args.infile[:-4] + '_DMHistogram.txt'
    bins_df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    # write DMtable (input for PeakSelector)
    outfile = args.infile[:-4] + '_DMTable.txt'
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
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
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/PeakModeller.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Input file with the peak file(s) to be filtered')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')

    parser.add_argument('-b', '--bins', help='Width of the bins')
    parser.add_argument('-p', '--points', help='Number of points (bins) to use for slope calculation')
    parser.add_argument('-s', '--spoints', help='Number of points (bins) to use for smoothing')

    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.bins is not None:
        config.set('PeakModeller', 'bins', str(args.bins))
        config.set('Logging', 'create_ini', '1')
    if args.points is not None:
        config.set('PeakModeller', 'points', str(args.points))
        config.set('Logging', 'create_ini', '1')
    if args.spoints is not None:
        config.set('PeakModeller', 'spoints', str(args.spoints))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/PeakModeller.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_log_debug.txt'
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