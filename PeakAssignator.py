#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Jose Rodriguez"
__credits__ = ["Jose Rodriguez", "Andrea Laguillo Gómez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.3.0"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es;andrea.laguillo@cnic.es"
__status__ = "Development"

# import modules
import os
import sys
import argparse
import configparser
import logging
import re
import pandas as pd
import concurrent.futures
from itertools import repeat
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

# TODO: allow user to set output column names in the INI

#feather_data = pd.read_feather(r"C:\Users\Andrea\Desktop\data.ftr")
#apex_list = _extract_ApexList(r"C:\Users\Andrea\Desktop\SHIFTS-4\testing\bin002_w7_SL4500Target_ApexList_mass.txt")

###################
# Local functions #
###################
def concatInfiles(peakpickingfile):
    '''    
    Concat input files...
    adding Experiment column (dirname of input file), and adding a FWHM columns by Experiment
    '''
    
    def _extract_FWHM(file):
        with open(file) as f:
            data = f.read()
            m = re.findall(r'FWHM:\s*([^\n]*)', data)
            if m and len(m)>0:
                return m[0]
            else:
                sys.exit("ERROR! FWHM is not defined for {}".format(file))
        
    # read input file
    # use high precision with the floats
    df = pd.read_csv(peakpickingfile, sep="\t", float_precision='high')
    # add folder name into column
    foldername = os.path.dirname(peakpickingfile)
    df['Experiment'] = foldername
    # add filename column
    df['Filename'] = os.path.basename(peakpickingfile)
    # add fwhm column
    #fwhm_file = "{}/{}".format(foldername, fwhm_fname)
    #fwhm = _extract_FWHM(fwhm_file)
    #fwhm = _extract_FWHM(fwhm_fname)
    #df['FWHM'] = float(fwhm)
    # assign type to categorical columns
    df['Experiment'] = df['Experiment'].astype('category')
    df['Filename'] = df['Filename'].astype('category')
    #df['Label'] = df['Label'].astype('category')
    #df['IsotpicJump'] = df['IsotpicJump'].astype('category')
    return df

def closest_peak(apex_list, delta_MH):
    '''
    Assign a delta_MH value to the closest apex in a list
    '''
    peak = min(apex_list, key = lambda x : abs(x - delta_MH))
    return peak

    
def find_orphans(ppm_max, theo_mass, peak, delta_MH, peak_label, orphan_label):
    '''
    Identify orphans and peaks
    '''
    # window = float(nsigma) * fwhm / 2
    distance = abs(peak - delta_MH)
    distance_ppm = (distance / (theo_mass + peak)) * 1e6
    #max_distance = abs(float(nsigma) * fwhm / 2)
    if distance_ppm <= ppm_max:
        ID = peak_label 
    else:
        ID = orphan_label
    # TODO: return ppm
    return ID, distance_ppm

# def get_deltamod(col_CalDeltaMH, col_Peak, orphan_label, col_ClosestPeak):
#     '''
#     Identify orphans and peaks
#     '''
#     if col_Peak == orphan_label:
#         deltamod = col_CalDeltaMH
#     else:
#         deltamod = col_ClosestPeak
#     return deltamod

def bin_operations(df, apex_list, ppm_max, peak_label, orphan_label,
                   col_ClosestPeak, col_CalDeltaMH, col_Peak, col_DM, col_TheoMass, col_ppm):
    '''
    Main function that handles the operations by BIN
    '''
    # get the BIN value from the input tuple df=(bin,df)
    (bin_value, df) = float(df[0]), df[1]
    
    # assign to peaks
    df[col_ClosestPeak] = df.apply(lambda x: closest_peak(apex_list, float(x[col_CalDeltaMH])), axis = 1)

    # identify orphans
    df[col_Peak] = df.apply(lambda x: find_orphans(ppm_max, float(x[col_TheoMass]), float(x[col_ClosestPeak]), float(x[col_CalDeltaMH]), peak_label, orphan_label)[0], axis = 1)
    df[col_ppm] = df.apply(lambda x: find_orphans(ppm_max, float(x[col_TheoMass]), float(x[col_ClosestPeak]), float(x[col_CalDeltaMH]), peak_label, orphan_label)[1], axis = 1)
    df[col_Peak] = df[col_Peak].astype('category')
    
    # calculate FDR
    
    # create deltamass column # TODO: Recom
    df[col_DM] = df.apply(lambda x: x[col_CalDeltaMH] if (x[col_Peak]==orphan_label) else x[col_ClosestPeak], axis = 1)
    #df[col_DM] = df.apply(lambda x: get_deltamod(x[col_CalDeltaMH], x[col_Peak], orphan_label, x[col_ClosestPeak]), axis = 1)
    
    # def peak_FDR():
      # for each peak sort by xcorr (comet) # should we separate recom peaks?
      # for each peak separate targets and decoys
      # make new fdr column rank_D/rank/T
    # def local_FDR():
      # for each bin sort by xcorr (comet)
      # for each bin separate targets and decoys
      # make new fdr_column rank_D/rank_T
    #def recom_FDR():
      # for each recom peak sort by xcorr (recom)
    
    # TO CHECK, print by BIN
    # outfile = os.path.join("D:/tmp/kk/", bin+"_kk.tsv")
    # df.to_csv(outfile, sep="\t", index=False)
    return df  

def format_seq(seqdm, dm, decimal_places):
    '''
    Make column with sequence and deltamass.    
    '''
    #df.apply(lambda x: x[seqdmcolumn].split('[')[0] + '[' + str(round(x[col_DM], decimal_places)) + ']' + x[seqdmcolumn].split(']')[1], axis = 1)
        
    if '[' in str(seqdm):
        formatseq = str(seqdm).split('[')[0] + '[' + str(round(float(dm), decimal_places)) + ']' + str(seqdm).split(']')[1]
    elif '_' in str(seqdm):
        formatseq = str(seqdm).split('_')[0] + '_' + str(round(float(dm), decimal_places))
    else:
        sys.exit("Unrecognized sequence format in '" + str(config._sections['PeakAssignator'][seqdm]) + "' column!")
    return formatseq
    
#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    # Variables
    ppm_max = abs(float(config._sections['PeakAssignator']['ppm_max']))
    col_TheoMass = config._sections['PeakAssignator']['theomh_column']
    col_CalDeltaMH = config._sections['PeakAssignator']['caldeltamh_column']
    col_ClosestPeak = config._sections['PeakAssignator']['closestpeak_column']
    col_Peak = config._sections['PeakAssignator']['peak_column']
    col_DM = config._sections['PeakAssignator']['deltamass_column']
    col_ppm = config._sections['PeakAssignator']['ppm_column']
    peak_label = config._sections['PeakAssignator']['peak_label']
    orphan_label = config._sections['PeakAssignator']['orphan_label']
    seqdmcolumn = config._sections['PeakAssignator']['mod_peptide_column']
    assignseqcolumn = config._sections['PeakAssignator']['assign_seq_column']
    decimal_places = int(config._sections['General']['decimal_places'])
    
    # logging.info("get the list of files with the inputs")
    # with open(args.infile) as f:
    #     infiles = f.readlines()
    # # you may also want to remove whitespace characters like `\n` at the end of each line
    # infiles = [x.strip() for x in infiles] 
    # logging.debug(infiles)
    
     # read apex list
    def _extract_ApexList(file):
        with open(file) as f:
            data = f.read().split('\n')
            data = [x for x in data if x.strip()]
            data = np.array(data, dtype=np.float64)
            return data
    
    #foldername = os.path.dirname(args.appfile)
    #apex_file = "{}/{}".format(foldername, args.appfile)
    #apex_list = _extract_ApexList(apex_file)
    apex_list = _extract_ApexList(args.appfile)

    #logging.info("concat input files...")
    #logging.info("adding Experiment column (dirname of input file),")
    #logging.info("and adding a FWHM column by Experiment")
    # with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
    #     df = executor.map(concatInfiles, infiles)
    # df = pd.concat(df)
    # df.reset_index(drop=True, inplace=True)
    if args.infile.strip().lower().endswith('.feather'):
        logging.info("Reading input feather file...")
        try:
            df = pd.read_feather(args.infile)
        except Exception as e:
            logging.error("ERROR: Reading the input feather file")
            sys.exit(1)
    else:
        logging.info("Reading input feather file...")
        try:
            df = pd.read_csv(args.infile, sep="\t", float_precision='high', low_memory=False)
        except Exception as e:
            logging.error("ERROR: Reading the input tabular file")
            sys.exit(1)

    logging.info("Create a column with the bin")
    df['bin'] = df[col_CalDeltaMH].astype(str).str.extract(r'^([^\.]*)')
    df['theo_mh'] = df['theo_mh'].astype(float)
    df[col_CalDeltaMH] = df[col_CalDeltaMH].astype(float)

    logging.info("Parallel the operations by bin")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:        
        df = executor.map(bin_operations, list(df.groupby("bin")), repeat(apex_list),
                                                                   repeat(ppm_max),
                                                                   repeat(peak_label),
                                                                   repeat( orphan_label),
                                                                   repeat(col_ClosestPeak),
                                                                   repeat(col_CalDeltaMH),
                                                                   repeat(col_Peak),
                                                                   repeat(col_DM),
                                                                   repeat(col_TheoMass),
                                                                   repeat(col_ppm)) 
    df = pd.concat(df)
    #logging.info("calculate gobal FDR")
    #df = get_global_FDR(df, args.xcorr)
    #logging.info("sort by DeltaMax cal")
    #df.sort_values(by=[col_CalDeltaMH], inplace=True)
    df.reset_index(drop=True, inplace=True)

    # d_h = df.head()
    # d_t = df.tail()
    # d_h.to_csv("kk_head.tsv", sep="\t")
    # d_t.to_csv("kk_tail.tsv", sep="\t")
    
    # Make assignseqcolumn # TODO: make new seqcolumnXXXXX_xxxx
    # flag = 0
    # if recomseqdmcolumn == "" or recomseqdmcolumn not in df:
    #     df["RECOMfiltered_type"] = "COMET"
    #     flag = 1
    #     recomseqdmcolumn = "placeholder"
    #     df[recomseqdmcolumn] = "placeholder"
    if assignseqcolumn not in df:
        df.insert(df.columns.get_loc(seqdmcolumn)+2, assignseqcolumn, np.nan) 
    df[assignseqcolumn] = df.apply(lambda x: format_seq(x[seqdmcolumn], x[col_DM], decimal_places), axis = 1)
    #df[assignseqcolumn] = df.apply(lambda x: x[seqdmcolumn].split('[')[0] + '[' + str(round(x[col_DM], decimal_places)) + ']' + x[seqdmcolumn].split(']')[1], axis = 1)
    # if flag == 1:
    #     df = df.drop(columns=["RECOMfiltered_type",recomseqdmcolumn])
    logging.info("Writing output files...")
    # https://towardsdatascience.com/the-best-format-to-save-pandas-data-414dca023e0d
    # begin:printHDF5
    # Note: Explote the Memory!!!
    # assign NumExpr for the tables module
    # tables.parameters.MAX_NUMEXPR_THREADS = args.n_workers
    # df.to_hdf('data.h5', key='df', mode='w')
    # end:printHDF5
    if args.infile.strip().lower().endswith('.feather'):
        outfile = args.infile[:-8] + '_PeakAssignation.feather'
        df.to_feather(outfile)
    else:
        outfile = args.infile[:-8] + '_PeakAssignation.tsv'
        df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    logging.info("Peak assignation finished.")
    

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Assignator',
        epilog='''
        Example:
            python PeakAssignator.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")
    
    parser.add_argument('-i',  '--infile', required=True, help='Input file with the peak picking')
    parser.add_argument('-a',  '--appfile', required=True, help='File with the apex list of Mass')
    #parser.add_argument('-f',  '--fwhm_filename', default='MAD_and_FWHM_calculations.txt', help='File name with the FWHM value (default: %(default)s)')    
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    #parser.add_argument('-mn', '--mindelta', help='Minimum Delta Mass (default: %(default)s)')
    #parser.add_argument('-mx', '--maxdelta', help='Maximum Delta Mass (default: %(default)s)')
    parser.add_argument('-p',  '--ppm', help='Maximum ppm difference for peak assignation')

    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    #if args.mindelta is not None:
        #config.set('PeakAssignator', 'mindelta', str(args.mindelta))
        #config.set('Logging', 'create_ini', '1')
    #if args.maxdelta is not None:
        #config.set('PeakAssignator', 'maxdelta', str(args.maxdelta))
        #config.set('Logging', 'create_ini', '1')
    if args.ppm is not None:
        config.set('PeakAssignator', 'ppm_max', str(args.ppm))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-8] + 'PeakAssignation_log.txt'
    log_file_debug = outfile = args.infile[:-8] + 'PeakAssignation_log_debug.txt'
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