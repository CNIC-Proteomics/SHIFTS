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
import logging
import re
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

def readHistogram(infile):
    df_hist = pd.read_csv(infile, sep="\t", float_precision='high')
    df_hist = df_hist.dropna() # Remove rows with missing values (will always have some in beginning and end)
    df_hist.reset_index(drop=True, inplace=True)
    return df_hist

def multipleApex(apex_list, apex_massdiff):
    diffs = np.diff(apex_list)
    new_apex_list = []
    for i in range(len(apex_list)):
        check = []
        if i-1 >= 0: check.append(diffs[i-1]) # not the first one
        if i <= len(diffs)-1: check.append(diffs[i]) # not the last one
        if all(diff <= apex_massdiff for diff in check):
            new_apex_list.append(apex_list[i])
    return new_apex_list

def firstAndLastApex(apex_list):
    new_apex_list = []
    new_apex_list.append(apex_list[0])
    new_apex_list.append(apex_list[-1])
    return new_apex_list

def peakSelector(df_hist, slope, frequency, apex_points, decimal_places):
    
    ### MARK BINS ###
    
    df_hist['previous'] = df_hist['slope1'].shift()
    df_hist['next'] = df_hist['slope1'].shift(-1)
    
    # Mark apex bins
    df_hist['apex'] = df_hist.apply(lambda x: 1 if (x['slope1']<0 and x['previous']>0)
                                                else 0, axis = 1)
    
    ### TEST ###
    df_hist['peak_begin'] = df_hist.apply(lambda x: 1 if (abs(x['slope1'])>slope and x['slope1']>0 and abs(x['previous'])<slope) #beginning
                                                else 0, axis = 1)
    df_hist['peak_end'] = df_hist.apply(lambda x: 1 if (abs(x['slope1'])>slope and x['slope1']<0 and abs(x['next'])<slope) #end
                                        else 0, axis = 1)
    df_hist['slope_threshold'] = df_hist.apply(lambda x: 1 if (x['peak_begin'] == 1 and abs(x['slope1'])>slope)
                                                           or (x['peak_end'] == 1 and abs(x['slope1'])>slope)
                                                           else 0, axis = 1)
            
    df_hist['peak_group'] = 0
    begin_list = df_hist.index[(df_hist['peak_begin'] == 1) & (df_hist['slope_threshold'] == 1)].tolist()
    for i in begin_list:
        around_apex = [i]
        new_index = i
        in_peak = True
        while in_peak == True:
            new_index += 1
            if new_index >= len(df_hist)-1:
                in_peak = False
                break
            if df_hist.loc[new_index]['peak_begin'] != 0 and df_hist.loc[new_index]['slope_threshold'] != 0:
                in_peak == False
                break
            if df_hist.loc[new_index]['peak_end'] != 0 and df_hist.loc[new_index]['slope_threshold'] != 0: #TODO: we never reach here?
                in_peak == False
                around_apex.extend(range(i, new_index+1))
                break
        for k in around_apex:
            df_hist.at[k, 'peak_group'] = 1
            
    df_hist = df_hist.drop(columns='previous')
    df_hist = df_hist.drop(columns='next')
    
    ### FILTER PEAKS ###
    
    grouped_hist = df_hist.groupby((df_hist['peak_group'].shift() != df_hist['peak_group']).cumsum())
    
    apex_bin_list = []
    for position, peak in grouped_hist:
        peak_df = peak
        if all(peak_df['peak_group'] != 0): #groups marked as peaks
            if any(peak_df['count'] >= frequency): #TODO fix for several apexes ## HERE I STOPPED
                if 1 in peak_df['apex'].value_counts().index and peak_df['apex'].value_counts()[1] == 1: #one apex
                    for i in peak_df['midpoint'].loc[peak_df['apex'] == 1]:
                        apex_bin_list.append(i)
                if 1 in peak_df['apex'].value_counts().index and peak_df['apex'].value_counts()[1] > 1: #more than one potential apex
                    #apex_bin_list.extend(multipleApex(list(peak_df['midpoint'].loc[peak_df['apex'] == 1]), apex_massdiff))
                    apex_bin_list.extend(firstAndLastApex(list(peak_df['midpoint'].loc[peak_df['apex'] == 1])))

    ### CALCULATE APEX ###
    apex_list = []
    before = apex_points//2
    after = (apex_points//2) - 1
    for apex_bin in apex_bin_list:
        bin_subset = df_hist.loc[df_hist['midpoint'] == apex_bin]
        try:
            for i in range(1, before + 1):
                bin_subset = pd.concat([bin_subset, pd.DataFrame(df_hist.loc[df_hist['midpoint'].shift(-i) == apex_bin])], ignore_index=True)
            for i in range(1, after + 1):
                bin_subset = pd.concat([bin_subset, pd.DataFrame(df_hist.loc[df_hist['midpoint'].shift(i) == apex_bin])], ignore_index=True)
            bin_subset.sort_values(by=['midpoint'], inplace=True)
            bin_subset.reset_index(drop=True, inplace=True)
            apex = interpolateApex(bin_subset)
            apex_list.append(round(apex, decimal_places))
        except:
            logging.info("Not enough bins to interpolate apex at" + str(apex_bin))
    
    return apex_list
    
def filterPeaks(df_hist, slope, frequency):
    '''
    Find peaks that are above the thresholds for slope and PSMs.
    '''
    # TODO: allow specify slope and count columns in INI?
    df_hist['apex'] = 0
    df_hist['previous'] = df_hist['slope1'].shift()
    df_hist['next'] = df_hist['slope1'].shift(-1)
    df_hist['apex'] = df_hist.apply(lambda x: 1 if (x['slope1']<0 and x['previous']>0) or (x['slope1']>0 and x['next']<0) else 0, axis = 1)
    df_hist = df_hist.drop(columns='previous')
    df_hist = df_hist.drop(columns='next')
    
    df_hist1 = df_hist[abs(df_hist['slope1']) >= slope] # keep those whose slope1 is over the threshold
    df_hist2 = df_hist[df_hist['apex'] == 1] # keep those where there is a sign change
    
    df_hist = pd.concat([df_hist1, df_hist2])
    df_hist.drop_duplicates(subset ="midpoint", keep = "first", inplace = True) 
    df_hist.sort_values(by=['midpoint'], inplace=True)
    df_hist.reset_index(drop=True, inplace=True)
    
    df_hist = df_hist[df_hist['count'] >= frequency]
    
    # outfile = r'D:\CNIC\SHIFTS-4\testCris\Cris\df_hist.txt'
    # df_hist.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    # print("Done filtering")
    
    return df_hist

def parseInterval(bins_df):
    '''
    Read 'bin' column as an interval.
    '''
    for i in range(0, len(bins_df)):
        to_interval = bins_df.loc[i, 'bin']
        left = float(re.findall(r'-?\d+\.\d+', to_interval)[0]) 
        right = float(re.findall(r'-?\d+\.\d+', to_interval)[1])
        bins_df.loc[i, 'bin'] = pd.Interval(left, right, closed='right')
    return bins_df

def areValid(intervals):
    '''
    Check whether intervals in a list are contiguous, have a change in
    sign of the slope, and the central point is the closest to 0.
    '''
    cont = 0
    slope1_list = intervals['slope1'].tolist()
    zero_bin = min(slope1_list, key=abs)
    zero_index = slope1_list.index(zero_bin)
    
    first_half = slope1_list[:len(slope1_list)//2]
    second_half = slope1_list[(len(slope1_list)//2)+1:]
    if all([x > 0 for x in first_half]) and all([x < 0 for x in second_half]):
        if zero_index == len(intervals)//2: # Central point is closest to 0
           if (intervals.loc[zero_index, 'slope1'] >= 0 and intervals.loc[zero_index+1, 'slope1'] < 0) or (intervals.loc[zero_index, 'slope1'] <= 0 and intervals.loc[zero_index-1, 'slope1'] > 0): # Change in sign  
               bin_list = intervals['bin'].tolist()
               cont = 1
               for i in range(1, len(bin_list)):
                   if bin_list[i-1].right != bin_list[i].left: # Not contiguous
                       cont = 0
    if cont == 0:
        return False
    else:
        return True
    
def interpolateApex(bin_subset):
    x_list = bin_subset['midpoint'].tolist()
    y_list = bin_subset['slope1'].tolist()
    sum1, sum2 = 0, 0
    for i in range(len(x_list)):
        sum1 += (x_list[i] - np.mean(x_list)) * (y_list[i] - np.mean(y_list))
        sum2 += (x_list[i] - np.mean(x_list)) ** 2
    working_slope = sum1 / sum2
    intercept = np.mean(y_list) - working_slope*np.mean(x_list)
    apex = -intercept / working_slope # x where y=0
    return apex

def peakApex(bins_df, apex_points):
    '''
    Calculate apex for each peak.
    '''
    apex_list = []
    # for i in range(1, len(bins_df)):
    #     if bins_df.loc[i, 'slope2'] is not None:
    #         i1 = bins_df.loc[i-1, 'bin'] # TODO parse interval
    #         i2 = bins_df.loc[i, 'bin']
    #         # Check intervals are consecutive, and there is a change in sign of slope2
    #         if i1.right == i2.left and bins_df.loc[i, 'slope2'] < 0 and bins_df.loc[i-1, 'slope2'] >= 0:
    #             peak = pd.Series([bins_df.loc[i-1, 'midpoint'], np.nan, bins_df.loc[i, 'midpoint']],
    #                               index=[bins_df.loc[i-1, 'slope2'], 0, bins_df.loc[i, 'slope2']])
    #             peak = peak.interpolate(method='index')
    #             apex_list.append(peak[0])
    for i in range(apex_points//2, len(bins_df)-apex_points//2):
        # Check there is a change of sign
        intervals = []
        for j in range(i-apex_points//2, i+apex_points//2+1):
            intervals.append(bins_df.loc[j])
        intervals = pd.DataFrame(intervals)
        intervals.reset_index(drop=True, inplace=True)
        if areValid(intervals):
            peak = interpolateApex(intervals)
            apex_list.append(peak)
    return apex_list

def main(args):
    '''
    Main function
    '''
    
    # Main variables
    decimal_places = int(config._sections['General']['decimal_places'])
    slope = float(config._sections['PeakSelector']['slope'])
    frequency = int(config._sections['PeakSelector']['frequency'])
    apex_points = int(config._sections['PeakSelector']['apex_points'])
    #apex_massdiff = float(config._sections['PeakSelector']['apex_massdiff'])
    
    # Read DM Histogram
    logging.info("Reading input file...")
    df_hist = readHistogram(args.infile)
    df_hist.reset_index(drop=True, inplace=True)
    
    # Filter by slope and frequency, calculate apexes
    logging.info("Filtering...")
    logging.info("Slope threshold = " + str(slope))
    logging.info("Frequency threshold = " + str(frequency))
    logging.info("Number of points to use for apex calculation = " + str(apex_points))
    # df_hist = filterPeaks(df_hist, slope, frequency)
    # df_hist.reset_index(drop=True, inplace=True)
    # df_hist = parseInterval(df_hist)
    # apex_list = peakApex(df_hist, apex_points)
    apex_list = peakSelector(df_hist, slope, frequency, apex_points, decimal_places)
    apex_info = str(len(apex_list)) + " peaks"
    logging.info(apex_info)
    
    # Write apex list
    logging.info("Writing apex list...")
    outfile = args.infile[:-15] + 'ApexList.txt'
    with open(outfile, 'w') as f:
        for apex in apex_list:
            f.write("%s\n" % apex)
    logging.info("Peak Selection finished")

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Selector',
        epilog='''
        Example:
            python PeakSelector.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='DMHistogram to be filtered')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    parser.add_argument('-s', '--slope', help='Threshold for slope of DM peak')
    parser.add_argument('-f', '--frequency', help='Threshold for number of PSMs')
    parser.add_argument('-p', '--apex_points', help='Number of points (bins) to use for apex calculation')
    #parser.add_argument('-a', '--apex_massdiff', help='Threshold for distance between apexes')

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.slope is not None:
        config.set('PeakSelector', 'slope', str(args.slope))
        config.set('Logging', 'create_ini', '1')
    if args.frequency is not None:
        config.set('PeakSelector', 'frequency', str(args.frequency))
        config.set('Logging', 'create_ini', '1')
    if args.apex_points is not None:
        config.set('PeakSelector', 'apex_points', str(args.apex_points))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-15] + 'ApexList_log.txt'
    log_file_debug = outfile = args.infile[:-15] + 'ApexList_log_debug.txt'
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