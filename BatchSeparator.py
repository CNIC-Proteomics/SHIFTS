# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 10:32:39 2021

@author: Andrea
"""

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.3.0"
__maintainer__ = "Andrea Laguillo Gómez"
__email__ = "jmrodriguezc@cnic.es;andrea.laguillo@cnic.es"
__status__ = "Development"

import argparse
import logging
import numpy as np
import os
import pandas as pd
import sys

def read_experiments(experiments_table):
    '''
    Read input file containing groups and filenames in tab-separated format.
    '''
    #df = pd.read_csv(experiments_table, sep="\t", names=['Batch', 'Experiment', 'Filename'])
    df = pd.read_csv(experiments_table, sep="\t", names=['Batch', 'Experiment'])
    df['Batch'] = df['Batch'].astype('string')
    df['Batch'] = df['Batch'].str.strip()
    df['Experiment'] = df['Experiment'].astype('string')
    df['Experiment'] = df['Experiment'].str.strip()
    if df['Experiment'].duplicated().any(): # Check no repeats
        sys.exit('ERROR: Batch table contains repeat values in the experiments column')
    #exp_groups = exp_df.groupby(by = exp_df.columns[0], axis = 0)
    #for position, exp in exp_groups:
        #TODO: read filepath or everything in folder
    return df

def make_groups(df, groups):
    '''
    Add Batch and Experiment columns to input file with the peak assignation.
    '''
    def _match_file(groups, experiment):
        if experiment in group_dict:
            group = group_dict.get(experiment)[0]
        else:
            group = 'N/A'
        return group
    df['Batch'] = 'N/A'
    #df['Experiment'] = df.apply(lambda x: _match_file(groups, x['Filename']), axis = 1)
    group_dict = {}
    for x in range(len(groups)):
        currentid = groups.iloc[x,1]
        currentvalue = groups.iloc[x,0]
        group_dict.setdefault(currentid, [])
        group_dict[currentid].append(currentvalue)
    df['Batch'] = np.vectorize(_match_file)(group_dict, df['Experiment'])
    if 'N/A' in df['Batch'].unique():
        unassigned = str(df['Batch']=='N/A'['Experiment'].unique().str.join(" ,"))
        logging.info('Warning: Experiment(s) ' +unassigned + ' could not be assigned to a batch!')
    return df

def main(args):
    '''
    Main function
    '''
    # Main variables
    logging.info('Read input file')
    df = pd.read_csv(args.infile, sep="\t", float_precision='high')
    
    logging.info('Read experiments table')
    groups = read_experiments(args.batch)
    df = make_groups(df, groups)
    
    logging.info("Write output files by batch:")
    dfs = df.groupby('Batch')
    for group in list(dfs.groups.keys()):
        group_path = os.path.join(args.output, group)
        if group == 'N/A':
            group_path = os.path.join(args.output, 'Unassigned')
        if not os.path.exists(group_path):
            os.mkdir(group_path)
        if group == 'N/A':
            outfile = os.path.join(group_path, args.infile.split('\\')[-1].split('/')[-1][:-4] + '_Unassigned_FDR.txt')
        else:
            outfile = os.path.join(group_path, args.infile.split('\\')[-1].split('/')[-1][:-4] + '_' + group + '_FDR.txt')
        group_df = dfs.get_group(group)
        logging.info('\t' + group + ': ' + str(outfile))
        group_df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')


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
    parser.add_argument('-b',  '--batch', required=True, help='Tab-separated file containing batch and experiment names')
    parser.add_argument('-o',  '--output', required=True, help='Output directory. Will be created if it does not exist')
    
    #parser.add_argument('-f',  '--fdr_filter', help='FDR value to filter by')
    #parser.add_argument('-t',  '--target_filter', help='Filter targets, 0=no 1=yes')
    #parser.add_argument('-r',  '--recom_data', help='Score for FDR calculation: 0=Xcorr, 1=cXcorr (default: %(default)s)')
   
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
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