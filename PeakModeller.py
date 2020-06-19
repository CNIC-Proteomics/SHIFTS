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
from itertools import repeat
pd.options.mode.chained_assignment = None  # default='warn'

#infile = r"C:\Users\Andrea\Desktop\SHIFTS-4\testing\cXcorr_Len_Rank_Results_TargetData_Calibration.txt"

###################
# Local functions #
###################
def concatInfiles(infile, fwhm_fname):
    '''    
    Concat input files...
    '''
  
    # read input file
    # use high precision with the floats
    df = pd.read_csv(infile, sep="\t", float_precision='high')
    # assign type to categorical columns
    df['filename'] = df['filename'].astype('category')
    df['Label'] = df['Label'].astype('category')
    df['IsotpicJump'] = df['IsotpicJump'].astype('category')
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
    df.sort_values(by=['Cal_Delta_MH'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    # make bins
    bins = list(np.arange(int(round(df['Cal_Delta_MH'][0])),
                          int(round(df['Cal_Delta_MH'].iloc[-1]))+bin_width,
                          bin_width))
    bins = [round(x, _decimal_places(bin_width)) for x in bins]
    df['bin'] = pd.cut(df['Cal_Delta_MH'], bins=bins)
    
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
    # TODO: define special cases at beginning and end and use a
    # linear regression function for the rest
    x_list = bin_subset['midpoint'].tolist()
    if smoothed:
        y_list = bin_subset['smooth_count'].tolist()
    elif not second_derivative:
        y_list = bin_subset['count'].tolist()
    else:
        y_list = bin_subset['Slope1'].tolist()
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

def smoothing(bins_df, points):
    '''
    Calculate the slope (first derivative) for each bin. Calculate new smoothed
    value for the midpoint using the linear regression line.
    '''
    bins_df['smooth_count'] = None
    for i in range(points, len(bins_df)-points): #TODO: handle leftovers at start/end
        #working_bin = bins_df.loc[i]
        bin_subset = bins_df[i-points:i+points+1]
        working_slope, intercept = linear_regression(bin_subset, False, False)
        bins_df.loc[i, 'smooth_count'] = intercept + (working_slope*bins_df.loc[i, 'midpoint'])
    bins_df[["smooth_count"]] = bins_df[["smooth_count"]].apply(pd.to_numeric)
    return bins_df

def first_derivative(bins_df, points):
    '''
    Calculate the slope (first derivative) for each bin.
    Returns the slope of the linear regression line through data points in
    known_y's and known_x's. The slope is the vertical distance divided by
    the horizontal distance between any two points on the line, which is the
    rate of change along the regression line.
    Known_y's  Bins. An array or cell range of numeric dependent data points.
    Known_x's  Count. The set of independent data points.
    '''
    bins_df = smoothing(bins_df, points)
    bins_df['Slope1'] = None
    for i in range(points*2, len(bins_df)-points*2): #TODO: handle leftovers at start/end
        #working_bin = bins_df.loc[i]
        bin_subset = bins_df[i-points:i+points+1]
        bins_df.loc[i, 'Slope1'] = linear_regression(bin_subset, True, False)
    bins_df[["Slope1"]] = bins_df[["Slope1"]].apply(pd.to_numeric)
    return bins_df

def second_derivative(bins_df, points):
    '''
    Calculate the second derivative for each bin.
    '''
    bins_df['Slope2'] = None
    for i in range(points*3, len(bins_df)-points*3): #TODO: handle leftovers at start/end
        bin_subset = bins_df[i-points:i+points+1]
        bins_df.loc[i, 'Slope2'] = linear_regression(bin_subset, False, True)
    bins_df[["Slope2"]] = bins_df[["Slope2"]].apply(pd.to_numeric)           
    return bins_df

def filter_peaks():
    '''
    Find peaks that are above the thresholds for slope and PSMs.
    '''
    return

#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    
    logging.info("read input file list")
    with open(args.infile) as f:
        infiles = f.readlines()
    infiles = [x.strip() for x in infiles] # remove whitespace
    
    logging.info("concat input files") # TODO: do we do this for each file or all together
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
        df = executor.map(concatInfiles, infiles, repeat(args.fwhm_filename))
    df = pd.concat(df)
    df.reset_index(drop=True, inplace=True)

    # make bins
    df, bins_df = generate_histogram(df, args.bins)
    # calculate derivatives
    #grouped_bins_df = bins_df.groupby(['bin'])
    bins_df = first_derivative(bins_df, args.points)  #does 1st smoothing pass and 2nd normal pass
    bins_df = second_derivative(bins_df, args.points)
        # check which bins pass
    # write apex list in txt

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
    
    parser.add_argument('-b', '--bins', help='Width of the bins')
    parser.add_argument('-p', '--points', help='Number of points (bins) to use for slope calculation')
    parser.add_argument('-s', '--slope', help='Threshold for slope of DM peak')
    parser.add_argument('-f', '--frequency', help='Threshold for number of PSMs')
    parser.add_argument('-m', '--mode', required=True, help='0=info, 1=filter')
    #parser.add_argument('-m', '--mode', required=True, help='0=filter by slope, 1=filter by frequency, 2=filter by both')
    # ALWAYS FILTER BY BOTH

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
    if args.slope is not None:
        config.set('PeakSelector', 'slope', str(args.slope))
        config.set('Logging', 'create_ini', '1')
    if args.frequency is not None:
        config.set('PeakSelector', 'frequency', str(args.frequency))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/DMcalibrator.ini', 'w') as newconfig:
            config.write(newconfig)
            
    # TODO: add mode to log!

    # logging debug level. By default, info level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')