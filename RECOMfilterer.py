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
import matplotlib.pyplot as plt
import scipy.optimize
import sys

#test file: r'D:\CNIC\asdf\RECOMFilterer\TMT2_HDL_ALL.txt'

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
    df['Label'] = df.apply(lambda x: 'Decoy' if (str(x[proteincolumn])[0:5]==decoyprefix) else 'Target', axis = 1)
    return df

def labelAD(df):
    '''
    Label increases and decreases according to DiffScore column. A = Increase, D = Decrease, N = close to 0 or not rescored.
    '''
    df['DiffType'] = df.apply(lambda x: 'D' if x['DiffScore']<1E-5 else ('A' if x['DiffScore']>1E-5 else 'N'), axis = 1)
    df['DiffType'] = df.apply(lambda x: 'N' if x['Closest_Xcorr']==0 else x['DiffType'], axis = 1)
    return df

def expFunction(x, b, c, a):
    '''
    Exponential decay function, Y=a+b*exp(-c*X).
    '''
    return b * np.exp(-c * x) + a

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
    
    # Exponential curve fitting
    xs = df['DiffScoreAbs'].to_numpy()
    ys = df['A/D_Ratio'].to_numpy()
    # perform the fit #TODO: let user choose starting parameters, using those as default
    popt, pcov = scipy.optimize.curve_fit(expFunction, xs, ys, p0=[0.5,5,0], maxfev=5000) # 5000 iterations should be enough but this limit might need to be adjusted
    
    # Estimate As
    df['Fitted_Curve'] = expFunction(df['DiffScoreAbs'], *popt)
    df['A_Est'] = df['Rank_D'] * df['Fitted_Curve']
    
    # Plot results
    fig1, plt1 = plt.subplots()
    plt1.plot(xs, ys, '.', label="A/D Ratio")
    plt1.plot(xs, expFunction(xs, *popt), 'r-', label="Fitted Curve") # The * in front of popt when you plot will expand out the terms into the a, b, and c that expFunction is expecting
    plt1.set_xlabel('DiffScoreAbs')
    plt1.set_ylabel('A/D Ratio')
    plt1.set_title("A/D Ratio vs DiffScore and theoretical curve")
    plt1.legend()
    #fig1.savefig(os.path.join(Path(args.output), Path(args.infile).stem + '_Ratio_vs_DiffScore.png'))
    fig1.savefig(args.infile[:-4] + '_Ratio_vs_DiffScore.png')
    
    fig2, plt2 = plt.subplots()
    plt2.plot(xs, df['Rank_D'].to_numpy(), '-', label="Rank D")
    plt2.plot(xs, df['Rank_A'].to_numpy(), '-', label="Rank A")
    plt2.plot(xs, df['A_Est'].to_numpy(), '-', label="Est. A")
    plt2.set_xlabel('DiffScoreAbs')
    plt2.set_ylabel('Rank')
    plt2.set_title("Rank D and A vs DiffScore (Decoys)")
    plt2.legend()
    #fig2.savefig(os.path.join(Path(args.output), Path(args.infile).stem + '_Rank_vs_DiffScore.png'))
    fig2.savefig(args.infile[:-4] + '_Rank_vs_DiffScore_Decoys.png')
    
    return popt

def checkTargets(df, popt):
    '''
    Check model using targets below score threshold.
    '''
    # Count As and Ds
    df['Rank'] = df.groupby('DiffType').cumcount()+1 # This column can be deleted later
    df['Rank_A'] = np.where(df['DiffType']=='A', df['Rank'], 0)
    df['Rank_A'] = df['Rank_A'].replace(to_replace=0, method='ffill')
    df['Rank_D'] = np.where(df['DiffType'] == 'D', df['Rank'], 0)
    df['Rank_D'] =  df['Rank_D'].replace(to_replace=0, method='ffill')
    df.drop(['Rank'], axis = 1, inplace = True)
    
    # Estimate As
    df['Fitted_Curve'] = expFunction(df['DiffScoreAbs'], *popt)
    df['A_Est'] = df['Rank_D'] * df['Fitted_Curve']
    
    # Plot results
    fig1, plt1 = plt.subplots()
    plt1.plot(df['DiffScoreAbs'].to_numpy(), df['Rank_D'].to_numpy(), '-', label="Rank D")
    plt1.plot(df['DiffScoreAbs'].to_numpy(), df['Rank_A'].to_numpy(), '-', label="Rank A")
    plt1.plot(df['DiffScoreAbs'].to_numpy(), df['A_Est'].to_numpy(), '-', label="Est. A")
    plt1.set_xlabel('DiffScoreAbs')
    plt1.set_ylabel('Rank')
    plt1.set_title("Rank D and A vs DiffScore (Low-scoring targets)")
    plt1.legend()
    #fig1.savefig(os.path.join(Path(args.output), Path(args.infile).stem + '_Rank_vs_DiffScore_Targets_Below_Threshold.png'))
    fig1.savefig(args.infile[:-4] + '_Rank_vs_DiffScore_Targets_Below_Threshold.png')
    
    return

def getDiffScoreCutOff(df, popt, t_increase):
    '''
    Use targets above score threshold to calculate DiffScoreCutOff that achieves "A est. FDR" < increase_threshold
    '''
    # Count As and Ds
    df['Rank'] = df.groupby('DiffType').cumcount()+1 # This column can be deleted later
    df['Rank_A'] = np.where(df['DiffType']=='A', df['Rank'], 0)
    df['Rank_A'] = df['Rank_A'].replace(to_replace=0, method='ffill')
    df['Rank_D'] = np.where(df['DiffType'] == 'D', df['Rank'], 0)
    df['Rank_D'] =  df['Rank_D'].replace(to_replace=0, method='ffill')
    df.drop(['Rank'], axis = 1, inplace = True)
    
    # Estimate As
    df['Fitted_Curve'] = expFunction(df['DiffScoreAbs'], *popt)
    df['A_Est'] = df['Rank_D'] * df['Fitted_Curve']
    
    # Plot results
    fig1, plt1 = plt.subplots()
    plt1.plot(df['DiffScoreAbs'].to_numpy(), df['Rank_D'].to_numpy(), '-', label="Rank D")
    plt1.plot(df['DiffScoreAbs'].to_numpy(), df['Rank_A'].to_numpy(), '-', label="Rank A")
    plt1.plot(df['DiffScoreAbs'].to_numpy(), df['A_Est'].to_numpy(), '-', label="Est. A")
    plt1.set_xlabel('DiffScoreAbs')
    plt1.set_ylabel('Rank')
    plt1.set_title("Rank D and A vs DiffScore (High-scoring targets)")
    plt1.legend()
    #fig1.savefig(os.path.join(Path(args.output), Path(args.infile).stem + '_Rank_vs_DiffScore_Targets_Above_Threshold.png'))
    fig1.savefig(args.infile[:-4] + '_Rank_vs_DiffScore_Targets_Above_Threshold.png')
    
    # Calculate "FDR"
    df['A_Est/A_Obs'] = df['A_Est'] / df['Rank_A']
    # TODO: calculate min DiffScore that still meets the threshold for FDR
    cutoff = df[df['A_Est/A_Obs']<=t_increase].tail(1)
    DiffScoreCutOff = cutoff.iloc[0]['DiffScoreAbs']
    return DiffScoreCutOff

def filterRECOM(df, dsco, a_dm, r_dm):
    '''
    Keep only RECOM IDs that pass the DiffScore threshold.
    '''
    df['RECOMfiltered_DM'] = df.apply(lambda x: x[r_dm] if x['DiffScore']<=dsco else x[a_dm], axis = 1)
    df['RECOMfiltered_type'] = df.apply(lambda x: 'RECOM' if x['DiffScore']<=dsco else 'COMET', axis = 1)
    df = df.drop(['DiffScore', 'DiffScoreAbs'], 1)
    recomized = df['RECOMfiltered_type'].value_counts()['RECOM']
    recomized_t = df[df['Label']=="Target"]['RECOMfiltered_type'].value_counts()['RECOM']
    recomized_d = df[df['Label']=="Decoy"]['RECOMfiltered_type'].value_counts()['RECOM']
    recomized = [recomized, recomized_t, recomized_d]
    return df, recomized

def main(args):
    '''
    Main function
    '''
    # Variables
    t_decoy = float(config._sections['RECOMfilterer']['decoy_threshold'])
    t_target = float(config._sections['RECOMfilterer']['target_threshold'])
    t_increase = float(config._sections['RECOMfilterer']['increase_threshold'])
    proteincolumn = config._sections['RECOMfilterer']['protein_column']
    assigneddm = config._sections['RECOMfilterer']['assigned_deltamass']
    recomdm = config._sections['RECOMfilterer']['recom_deltamass']
    decoyprefix = config._sections['RECOMfilterer']['decoyprefix']
    recom_score = config._sections['RECOMfilterer']['recom_score']
    comet_score = config._sections['RECOMfilterer']['comet_score']
    decimal_places = int(config._sections['General']['decimal_places'])
    
    # Read infile
    logging.info("Reading input file...")
    df = readInfile(Path(args.infile)) # df = readInfile(Path(infile))
    
    # Separate targets and decoys
    df = labelTargetDecoy(df, proteincolumn, decoyprefix) # df = labelTargetDecoy(df, 'protein', 'DECOY')
    df['DiffScore'] = df[recom_score] - df[comet_score] # df['DiffScore'] = df['Closest_Xcorr'] - df['xcorr']
    df['DiffScoreAbs'] = abs(df['DiffScore'])
    df = labelAD(df)
    targets = df[df['Label']=="Target"]
    decoys = df[df['Label']=="Decoy"]
    logging.info("Targets: " + str(targets.shape[0]) + " | Decoys: " + str(decoys.shape[0]))
    
    # True decoys
    true_decoys = decoys[decoys[recom_score]<=t_decoy] # TODO: are we filtering by Recom or Comet score?
    true_decoys = true_decoys[true_decoys['DiffType']!='N'] # Don't use values close to 0 or not rescored by Recom
    logging.info("Decoys <= " + str(t_decoy) + ": " + str(true_decoys.shape[0])) # true_decoys = decoys[decoys['Closest_Xcorr']<=1]
    true_decoys.sort_values(by=['DiffScoreAbs'], ascending=False, inplace=True) # Sort by descending abs. DiffScore
    true_decoys.reset_index(drop=True, inplace=True)
    # Make model
    popt = modelDecoys(true_decoys) # popt is an array
    logging.info("Model parameters for Y = a+b*exp(-c*X): a = " + str(round(popt[0], decimal_places)) + " b = " + str(round(popt[1], decimal_places)) + " c = " + str(round(popt[2], decimal_places)))
    
    # Fake targets
    fake_targets = targets[targets[recom_score]<=t_decoy]
    fake_targets = fake_targets[fake_targets['DiffType']!='N'] # Don't use values close to 0 or not rescored by Recom
    logging.info("Targets <= " + str(t_decoy) + ": " + str(fake_targets.shape[0])) # true_decoys = decoys[decoys['Closest_Xcorr']<=1]
    fake_targets.sort_values(by=['DiffScoreAbs'], ascending=False, inplace=True) # Sort by descending abs. DiffScore
    fake_targets.reset_index(drop=True, inplace=True)
    # Check model
    checkTargets(fake_targets, popt)
    
    # True targets
    true_targets = targets[targets[recom_score]>=t_target]
    true_targets = true_targets[true_targets['DiffType']!='N'] # Don't use values close to 0 or not rescored by Recom
    logging.info("Targets >= " + str(t_target) + ": " + str(true_targets.shape[0]))
    true_targets.sort_values(by=['DiffScoreAbs'], ascending=False, inplace=True) # Sort by descending abs. DiffScore
    true_targets.reset_index(drop=True, inplace=True)
    # Apply model
    dsco = getDiffScoreCutOff(true_targets, popt, t_increase)
    logging.info("Increase threshold: " + str(t_increase) + " | DiffScoreCutOff: " + str(round(dsco, decimal_places)))
    
    # Apply DiffScoreCutOff (0.05 or t_increase) to whole df
    # Choose which Recom improvements to keep: only those that pass DiffScoreCutOff. Otherwise we keep Comet
    logging.info("Filtering by DiffScoreCutOff...")
    df, recomized = filterRECOM(df, dsco, assigneddm, recomdm)
    logging.info("RECOMized " + str(recomized[0]) + " PSMs (" + str(round((recomized[0]/df.shape[0])*100, 2)) + " % of total PSMs)")
    logging.info("\t\t" + str(recomized[1]) + " Targets")
    logging.info("\t\t" + str(recomized[2]) + " Decoys")

    # Write to file
    logging.info("Writing output file...")
    outfile = args.infile[:-4] + '_RECOMfiltered.txt'
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')   
    logging.info("Recom filtering finished.")
    
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
    #parser.add_argument('-o',  '--output', required=True, help='Output directory')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    parser.add_argument('-d', '--decoy', help='Decoy threshold')
    parser.add_argument('-t', '--target', help='Target threshold')
    parser.add_argument('-s', '--increase', help='Significant increase threshold')

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.decoy is not None:
        config.set('RECOMfilterer', 'decoy_threshold', str(args.percentage))
        config.set('Logging', 'create_ini', '1')
    if args.target is not None:
        config.set('RECOMfilterer', 'target_threshold', str(args.cometcolumn))
        config.set('Logging', 'create_ini', '1')
    if args.increase is not None:
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