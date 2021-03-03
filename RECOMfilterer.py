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
import logging
import numpy as np
import os
import pandas as pd
from pathlib import Path
import sys

def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high')
    return df

def labelTargetDecoy(df, proteincolumn, decoyprefix):
    '''
    Label targets and decoys according to protein ID column.
    '''
    if 'Label' not in df:
        df.insert(df.columns.get_loc(proteincolumn)+1, 'Label', np.nan)
    df['Label'] = df.apply(lambda x: 'Decoy' if (x[proteincolumn][0:5]==decoyprefix) else 'Target', axis = 1)
    return df

def labelAD(df):
    '''
    Label increases and decreases according to DiffScore column.
    '''
    df['DiffType'] = df.apply(lambda x: 'D' if x['DiffScore']<1E-5 else ('A' if x['DiffScore']>1E-5 else 'N'), axis = 1)
    return df

def modelDecoys(df):
    '''
    Model increases as a function of decreases in an "FDR-like" way, using rank vs. DiffScore.
    '''
    # Count As and Ds
    df['Rank'] = df.groupby('DiffType').cumcount()+1 # This column can be deleted later
    df['Rank_A'] = np.where(df['DiffType']=='A', df['Rank'], 0)
    df['Rank_A'] = df['Rank_A'].replace(to_replace=0, method='ffill')
    df['Rank_D'] = np.where(df['DiffType'] == 'D', df['Rank'], 0)
    df['Rank_D'] =  df['Rank_D'].replace(to_replace=0, method='ffill')
    df.drop(['Rank'], axis = 1, inplace = True)
    
    # Calculate A/D ratio ("FDR")
    df['A/D_Ratio'] = df['Rank_A']/df['Rank_D']
    
    
    return df

def main(args):
    '''
    Main function
    '''
    # Variables
    t_decoy = float(config._sections['RECOMfilterer']['decoy_threshold'])
    t_target = float(config._sections['RECOMfilterer']['target_threshold'])
    t_increase = float(config._sections['RECOMfilterer']['target_threshold'])
    proteincolumn = config._sections['DMcalibrator']['proteincolumn']
    decoyprefix = config._sections['DMcalibrator']['decoyprefix']
    recom_score = config._sections['RECOMfilterer']['recom_score']
    comet_score = config._sections['RECOMfilterer']['comet_score']
    
    # Read infile
    logging.info("Reading input file...")
    df = readInfile(Path(args.infile))
    
    # Separate targets and decoys
    df = labelTargetDecoy(df, proteincolumn, decoyprefix)
    df['DiffScore'] = df[recom_score] - df[comet_score]
    df['DiffScoreAbs'] = abs(df['DiffScore'])
    df = labelAD(df)
    targets = df[df['Label']=="Target"]
    decoys = df[df['Label']=="Decoy"]
    logging.info("Targets: " + str(targets.shape[0]) + " | Decoys: " + str(decoys.shape[0]))
    
    # True decoys
    true_decoys = decoys[decoys[recom_score]<=t_decoy]
    logging.info("Decoys <= " + str(t_decoy) + ": " + str(true_decoys.shape[0]))
    true_decoys = true_decoys[true_decoys['DiffType']!='N'] # Don't use values close to 0
    true_decoys.sort_values(by=['DiffScoreAbs'], ascending=False, inplace=True) # Sort by descending abs. DiffScore
    true_decoys.reset_index(drop=True, inplace=True)
    # Make model
    
    
    # Fake targets
    fake_targets = targets[targets[recom_score]<=t_decoy]
    # Check model (graph?)
    
    
    # True targets
    true_targets = targets[targets[recom_score]>=t_target]
    logging.info("Targets >= " + str(t_target) + ": " + str(true_targets.shape[0]))
    # Apply model
    
    # Apply DiffScoreCutOff (0.05 or t_increase) to whole df
    
    
    # Choose which Recom improvements to keep: only those that pass DiffScoreCutOff. Otherwise we keep Comet
    
    
    # Write to file    
    
    
if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='RECOM Filterer',
        epilog='''
        Example:
            python RECOMfilterer.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")
    
    parser.add_argument('-i',  '--infile', required=True, help='Input file')  
    
    parser.add_argument('-d', '--decoy', help='Decoy threshold')
    parser.add_argument('-t', '--target', help='Target threshold')
    parser.add_argument('-s', '--increase', help='Significant increase threshold')

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.percentage is not None:
        config.set('RECOMfilterer', 'decoy_threshold', str(args.percentage))
        config.set('Logging', 'create_ini', '1')
    if args.cometcolumn is not None:
        config.set('RECOMfilterer', 'target_threshold', str(args.cometcolumn))
        config.set('Logging', 'create_ini', '1')
    if args.recomcolumn is not None:
        config.set('RECOMfilterer', 'increase_threshold', str(args.recomcolumn))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1: #TODO: check that other modules will not overwrite
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_RECOMfilterer_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_RECOMfilterer_log_debug.txt'
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