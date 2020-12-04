#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Andrea Laguillo Gómez"
__email__ = "jmrodriguezc@cnic.es;andrea.laguillo@cnic.es"
__status__ = "Development"

# import modules
import argparse
import configparser
import concurrent.futures
from itertools import repeat
import logging
import numpy as np
import os
import pandas as pd
import sys
import math
pd.options.mode.chained_assignment = None  # default='warn'

# TODO recom parameter: if no recom column specified, proceed without spire FDR

###################
# Local functions #
###################

def read_experiments(experiments_table):
    '''
    Read input file containing groups and filenames in tab-separated format.
    '''
    df = pd.read_csv(experiments_table, sep="\t", names=['Experiment', 'Filename'])
    df['Experiment'] = df['Experiment'].astype('string')
    df['Filename'] = df['Filename'].astype('string')
    if df['Filename'].duplicated().any(): # Check no repeats
        sys.exit("Experiments table contains repeat values in the filename column")
    #exp_groups = exp_df.groupby(by = exp_df.columns[0], axis = 0)
    #for position, exp in exp_groups:
        #TODO: read filepath or everything in folder
    return df

def make_groups(df, groups):
    '''
    Add group column to input file with the peak assignation.
    '''
    def _match_file(groups, filename):
        # if filename in groups['Filename'].unique():
            # group = df.loc[df['Filename'] == filename]['Experiment']
            # group.reset_index(drop=True, inplace=True)
            # group = group[0]
        # if filename in [x for v in group_dict.values() for x in v]:
        if filename in group_dict:
            group = group_dict.get(filename)[0]
        else:
            group = 'N/A'
        return group
    df['Experiment'] = 'N/A'
    #df['Experiment'] = df.apply(lambda x: _match_file(groups, x['Filename']), axis = 1)
    group_dict = {}
    for x in range(len(groups)):
        currentid = groups.iloc[x,1]
        currentvalue = groups.iloc[x,0]
        group_dict.setdefault(currentid, [])
        group_dict[currentid].append(currentvalue)
    df['Experiment'] = np.vectorize(_match_file)(group_dict, df['Filename'])
    if 'N/A' in df['Experiment'].unique():
        logging.info('Warning: ' + str(df['Experiment'].value_counts()['N/A']) + ' rows could not be assigned to an experiment!') # They will all be grouped together for FDR calculations
    return df

def get_spire_FDR(df, score_column, col_Peak, xcorr_type): #TODO: we don't have xcorr_type, we have recom_data, take out column names
    #This will be for the group of scans in a peak that are contained within 
    #one recom-assignated theoretical deltamass. Then, when we do peak_FDR, we
    #include these as well as the rest of the values in the peak.
    #Cuando el peak y el spire se solapen, esos escanes lo sometería a ambas FDR
    # How we will handle filtering is still to be determined.
    '''
    Calculate spire FDR for each spire in one bin (1 Da)
    '''
    df['SpireFDR'] = 100
    df['Rank'] = -1
    df['Spire_Rank_T'] = -1
    df['Spire_Rank_D'] = -1
    # TODO: Operations
    # identify spires (filter by recom HERE OR IN RECOM?)
    peaks = df[df['Peak'] == 'PEAK'] # filter by Peak
    recom_peaks = peaks[peaks['XcorType'] == 'RECOM'] # TODO: XcorType needs to be created
    #recom-identified scans are not necessarily peaks, what to do?
    grouped_recom_peaks = recom_peaks.groupby(['ClosestPeak']) # group by ClosestPeak
    for group in grouped_recom_peaks:
        group_index = group[1].index.values
        df.loc[group_index] # repeat steps of local_FDR
        # sort bin
        if xcorr_type == 0: # by Comet Xcorr
            df.loc[group_index].sort_values(by=['Xcor', 'Label'], inplace=True, ascending=False)
        else: # by Comet cXcorr
            df.loc[group_index].sort_values(by=['CorXcor', 'Label'], inplace=True, ascending=False) # TODO: Fix SHIFTS cXcorr
        # count targets and decoys
        df.loc[group_index]['Rank'] = df.loc[group_index].groupby('Label').cumcount()+1 # This column can be deleted later
        df.loc[group_index]['Spire_Rank_T'] = np.where(df.loc[group_index]['Label']=='Target', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Spire_Rank_T'] = df.loc[group_index]['Spire_Rank_T'].replace(to_replace=0, method='ffill')
        df.loc[group_index]['Spire_Rank_D'] = np.where(df.loc[group_index]['Label'] == 'Decoy', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Spire_Rank_D'] =  df.loc[group_index]['Spire_Rank_D'].replace(to_replace=0, method='ffill')
        # calculate local FDR
        df.loc[group_index]['SpireFDR'] = df.loc[group_index]['Spire_Rank_D']/df.loc[group_index]['Spire_Rank_T']
    # TODO: End Operations
    df.drop(['Rank'], axis = 1, inplace = True)
    return df

def get_peak_FDR(df, score_column, col_Peak, closestpeak_column, recom_data):
    '''
    Calculate peak FDR for each peak in one bin (1 Da)
    '''
    df['PeakFDR'] = 100
    df['Rank'] = -1
    df['Peak_Rank_T'] = -1
    df['Peak_Rank_D'] = -1
    # identify peaks
    peaks = df[df[col_Peak] == 'PEAK'] # filter by Peak
    grouped_peaks = peaks.groupby([closestpeak_column]) # group by ClosestPeak
    # df.get_group("group")
    #grouped_peaks.groups # group info
    # for group in grouped_peaks:
    #     group_index = group[1].index.values
    #     df.loc[group_index] # repeat steps of local_FDR
    #     # sort bin
    #     # if recom_data == 0: # by Comet Xcorr
    #     #     df.loc[group_index].sort_values(by=['Xcor', 'Label'], inplace=True)
    #     # else: # by Comet cXcorr
    #     #     df.loc[group_index].sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
    #     df.loc[group_index].sort_values(by=[score_column, 'Label'], inplace=True)
    #     # count targets and decoys
    #     df.loc[group_index]['Rank'] = df.loc[group_index].groupby('Label').cumcount()+1 # This column can be deleted later
    #     df.loc[group_index]['Peak_Rank_T'] = np.where(df.loc[group_index]['Label']=='Target', df.loc[group_index]['Rank'], 0)
    #     df.loc[group_index]['Peak_Rank_T'] = df.loc[group_index]['Peak_Rank_T'].replace(to_replace=0, method='ffill')
    #     df.loc[group_index]['Peak_Rank_D'] = np.where(df.loc[group_index]['Label'] == 'Decoy', df.loc[group_index]['Rank'], 0)
    #     df.loc[group_index]['Peak_Rank_D'] =  df.loc[group_index]['Peak_Rank_D'].replace(to_replace=0, method='ffill')
    #     # calculate peak FDR
    #     df.loc[group_index]['PeakFDR'] = df.loc[group_index]['Peak_Rank_D']/df.loc[group_index]['Peak_Rank_T']
        
    ###
    def _peak_FDR(group, score_column):
        group.sort_values(by=[score_column, 'Label'], inplace=True, ascending=False)
        group['Rank'] = group.groupby('Label').cumcount()+1 # This column can be deleted later
        group['Peak_Rank_T'] = np.where(group['Label']=='Target', group['Rank'], 0)
        group['Peak_Rank_T'] = group['Peak_Rank_T'].replace(to_replace=0, method='ffill')
        group['Peak_Rank_D'] = np.where(group['Label'] == 'Decoy', group['Rank'], 0)
        group['Peak_Rank_D'] =  group['Peak_Rank_D'].replace(to_replace=0, method='ffill')
        # calculate peak FDR
        group['PeakFDR'] = group['Peak_Rank_D']/group['Peak_Rank_T']
        return group
    peaks_df = []
    for group in grouped_peaks:
        peak_df = _peak_FDR(group[1], score_column)
        peaks_df.append(peak_df)
    if len(peaks_df) > 0:
        final_peaks_df = pd.concat(peaks_df)
        # join with df
        common_index = df.index.intersection(final_peaks_df.index)
        common_columns = df.columns.intersection(final_peaks_df.columns)
        df.loc[common_index, common_columns] = final_peaks_df.loc[common_index, common_columns]
    ###
    
    df.drop(['Rank'], axis = 1, inplace = True)
    return df

def get_local_FDR(df, score_column, recom_data):
    '''
    Calculate local FDR for one bin (1 Da)
    '''
    # sort bin
    #if recom_data == 0: # by Comet Xcorr
        #df.sort_values(by=['Xcor', 'Label'], inplace=True, ascending=False)
    #else: # by Comet cXcorr
        #df.sort_values(by=['CorXcor', 'Label'], inplace=Tru, ascending=False) # TODO: Fix SHIFTS cXcorr
    df.sort_values(by=[score_column, 'Label'], inplace=True, ascending=False)
        
    # count targets and decoys
    df['Rank'] = df.groupby('Label').cumcount()+1 # This column can be deleted later
    df['Local_Rank_T'] = np.where(df['Label']=='Target', df['Rank'], 0)
    df['Local_Rank_T'] = df['Local_Rank_T'].replace(to_replace=0, method='ffill')
    df['Local_Rank_D'] = np.where(df['Label'] == 'Decoy', df['Rank'], 0)
    df['Local_Rank_D'] =  df['Local_Rank_D'].replace(to_replace=0, method='ffill')
    df.drop(['Rank'], axis = 1, inplace = True)
    
    # calculate local FDR
    df['LocalFDR'] = df['Local_Rank_D']/df['Local_Rank_T']
    return df

def get_global_FDR(df, score_column, recom_data, peak_label, col_Peak, closestpeak_column, dm_column, dm_region_limit, n_workers):
    '''
    Calculate global FDR
    '''
    # get the EXPERIMENT value from the input tuple df=(experiment,df)
    (experiment_value, df) = df[0], df[1]
    print("\t\t\t\t\tCalculating Global/Local/Peak FDR for: " + experiment_value)
    # sort by score
    # if recom_data == 0: # by Comet Xcorr
    #     df.sort_values(by=['Xcor', 'Label'], inplace=True, ascending=False)
    # else: # by Comet cXcorr
    #     df.sort_values(by=['CorXcor', 'Label'], inplace=True, ascending=False) # TODO: Fix SHIFTS cXcorr
    
    #### TODO: make two regions separated by dm_region_limit ####
    #df.sort_values(by=[dm_column], inplace=True)
    #df.reset_index(drop=True, inplace=True)
    df_below =  df.loc[df[dm_column] < dm_region_limit]
    df_above = df.loc[df[dm_column] >= dm_region_limit]
    df_list = [df_below, df_above]
    #############################################################
    
    for each_df in df_list:
        each_df.sort_values(by=[score_column, 'Label'], inplace=True, ascending=False)
            
        # count targets and decoys
        each_df['Rank'] = each_df.groupby('Label').cumcount()+1 # This column can be deleted later
        each_df['Global_Rank_T'] = np.where(each_df['Label']=='Target', each_df['Rank'], 0)
        each_df['Global_Rank_T'] = each_df['Global_Rank_T'].replace(to_replace=0, method='ffill')
        each_df['Global_Rank_D'] = np.where(each_df['Label'] == 'Decoy', each_df['Rank'], 0)
        each_df['Global_Rank_D'] =  each_df['Global_Rank_D'].replace(to_replace=0, method='ffill')
        each_df.drop(['Rank'], axis = 1, inplace = True)
        
        # calculate global FDR
        each_df['GlobalFDR'] = each_df['Global_Rank_D']/each_df['Global_Rank_T']
    
    df = pd.concat([df_below, df_above])
    
    # calculate local and peak FDR
    #logging.info("Calculating Local and Peak FDR for: " + experiment_value)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:        
        df = executor.map(bin_operations, list(df.groupby('LocalBin')), repeat(score_column),
                                                                   repeat(recom_data), 
                                                                   repeat(peak_label),
                                                                   repeat(col_Peak),
                                                                   repeat(closestpeak_column)) 
    df = pd.concat(df)
    
    return df

def filtering(df, fdr_filter, target_filter): # This goes on a separate module now
    if target_filter: # =! 0
        df[df['Label'] == 'Target']
    if fdr_filter: # =! 0
        df[df['GlobalFDR'] >= fdr_filter]
    return df

def bin_operations(df, score_column, recom_data, peak_label, col_Peak, closestpeak_column):
    '''
    Main function that handles the operations by BIN
    '''
    
    # get the BIN value from the input tuple df=(bin,df)
    (bin_value, df) = df[0], df[1]
    
    # calculate local FDR
    df = get_local_FDR(df, score_column, recom_data)
    
    # calculate peak FDR
    df = get_peak_FDR(df, score_column, col_Peak, closestpeak_column, recom_data)
    
    # calculate spire FDR
    #if recom_data: #recom_data =! 0
    #df = get_spire_FDR(df, score_column, col_Peak, recom_data)
    
    return df

def make_bins(col_CalDeltaMH):
    '''
    Make bins for local FDR, centered at .5 Da
    '''
    bin_width = 1 #Da
    decimal, deltamass = math.modf(float(col_CalDeltaMH))
    if deltamass >= 0:
        if abs(decimal) >= 0.5:
            local_bin = deltamass + 0.5
        else:
            local_bin = deltamass -0.5
    else:
        if abs(decimal) >= 0.5:
            local_bin = deltamass - 1.5
        else:
            local_bin = deltamass - 0.5
    local_bin_str = str(local_bin) + " to " + str(local_bin + bin_width)
    return local_bin_str
    

#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    # Main variables
    n_workers = args.n_workers
    score_column = config._sections['PeakFDRer']['score_column']
    dm_column = config._sections['PeakFDRer']['dm_column']
    dm_region_limit = float(config._sections['PeakFDRer']['dm_region_limit'])
    recom_data = config._sections['PeakFDRer']['recom_data']
    peak_label = config._sections['PeakAssignator']['peak_label']
    col_Peak = config._sections['PeakAssignator']['peak_column']
    col_CalDeltaMH = config._sections['PeakAssignator']['caldeltamh_column']
    closestpeak_column = config._sections['PeakAssignator']['closestpeak_column']
    # fdr_filter = config._sections['PeakFDRer']['fdr_filter']
    # target_filter = config._sections['PeakFDRer']['target_filter']
    
    # try:
    #     if not os.path.exists(args.output):
    #         os.makedirs(args.output)
    #         logging.info("Create output directory at %s " % args.output)
    # except OSError:
    #     sys.exit("Could not create output directory at %s" % args.output)
    
    # Read input file
    logging.info('Read input file')
    #df = pd.read_feather(args.infile)
    df = pd.read_csv(args.infile, sep="\t", float_precision='high')
    
    # Add groups
    logging.info('Read experiments table')
    groups = read_experiments(args.experiment_table)
    df = make_groups(df, groups)
    # Return info
    group_dict = {a: b['Filename'].tolist() for a,b in groups.groupby('Experiment')}
    for key in group_dict:
        logging.info('\t' + key + ': ' + str(len(group_dict[key])) + ' files')
    
    logging.info("Binning")
    df['LocalBin'] = np.vectorize(make_bins)(df[col_CalDeltaMH])
    # logging.info("Calculate Peak and Local FDR")
    # with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:        
    #     df = executor.map(bin_operations, list(df.groupby('LocalBin')), repeat(score_column),
    #                                                                repeat(recom_data), 
    #                                                                repeat(peak_label),
    #                                                                repeat(col_Peak),
    #                                                                repeat(closestpeak_column)) 
    # df = pd.concat(df)
    #df.drop(['LocalBin'], axis = 1, inplace = True)
    
    logging.info("Calculate FDR")
    # df = get_global_FDR(df, score_column, recom_data)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
        df = executor.map(get_global_FDR, list(df.groupby('Experiment')), repeat(score_column),
                                                                          repeat(recom_data),
                                                                          repeat(peak_label),
                                                                          repeat(col_Peak),
                                                                          repeat(closestpeak_column),
                                                                          repeat(dm_column),
                                                                          repeat(dm_region_limit),
                                                                          repeat(n_workers))
    df = pd.concat(df)
    
    logging.info("Sort by calibrated DeltaMass")
    df.sort_values(by=[col_CalDeltaMH], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    # TODO: groups?????
    
    # Filtering # This goes on a separate module now
    # df = filtering(df, fdr_filter, target_filter)
    # df.reset_index(drop=True, inplace=True)

    # d_h = df.head()
    # d_t = df.tail()
    # d_h.to_csv("kk_head.tsv", sep="\t")
    # d_t.to_csv("kk_tail.tsv", sep="\t")
    
    # Clean filename
    df['Filename'] = df.apply(lambda x: x['Filename'].split('\\')[-1].split('/')[-1][:-4], axis=1)
    # Split in folders by Experiment
        
    logging.info("Write output files")
    dfs = df.groupby('Experiment')
    for group in list(dfs.groups.keys()):
        group_path = args.output + group
        if not os.path.exists(group_path):
            os.mkdir(group_path)
        outfile = os.path.join(group_path, args.infile.split('\\')[-1].split('/')[-1][:-4] + '_' + group + '_FDR.txt')
        group_df = dfs.get_group(group)
        group_df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    
    #outfile = args.infile[:-4] + '_FDR.txt'
    #df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    

    

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak FDRer',
        epilog='''
        Example:
            python PeakFDRer.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")
    
    parser.add_argument('-i',  '--infile', required=True, help='Input file with the peak assignation')
    parser.add_argument('-e',  '--experiment_table', required=True, help='Tab-separated file containing experiment names and file paths')
    parser.add_argument('-c',  '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-o',  '--output', required=True, help='Output directory. Will be created if it does not exist')
    
    parser.add_argument('-s',  '--score_column', help='Name of column with score for FDR calculation')
    #parser.add_argument('-f',  '--fdr_filter', help='FDR value to filter by')
    #parser.add_argument('-t',  '--target_filter', help='Filter targets, 0=no 1=yes')
    parser.add_argument('-r',  '--recom_data', help='Score for FDR calculation: 0=Xcorr, 1=cXcorr (default: %(default)s)')

    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.score_column is not None:
        config.set('PeakFDRer', 'score_column', str(args.score_column))
        config.set('Logging', 'create_ini', '1')
    # if args.fdr_filter is not None:
    #     config.set('PeakFDRer', 'fdr_filter', str(args.fdr_filter))
    #     config.set('Logging', 'create_ini', '1')
    # if args.target_filter is not None:
    #     config.set('PeakFDRer', 'target_filter', str(args.target_filter))
    #     config.set('Logging', 'create_ini', '1')
    if args.recom_data is not None:
        config.set('PeakFDRer', 'recom_data', str(args.recom_data))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)
    
    created = 0
    try:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
            created = 1
    except OSError:
        sys.exit("Could not create output directory at %s" % args.output)

    # logging debug level. By default, info level
    #log_file = args.infile[:-4] + '_FDR_log.txt'
    log_file = os.path.join(args.output, args.infile.split('\\')[-1].split('/')[-1][:-4] + '_FDR_log.txt')
    log_file_debug = os.path.join(args.output, args.infile.split('\\')[-1].split('/')[-1][:-4] + '_FDR_log_debug.txt')
    #log_file_debug = args.infile[:-4] + '_FDR_log_debug.txt'
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
    if created == 1:
        logging.info("Created output directory at %s " % args.output)
    main(args)
    logging.info('end script')