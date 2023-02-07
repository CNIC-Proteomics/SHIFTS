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
import argparse
import glob
import logging
import pandas as pd
import sys

def main(args):
    '''
    Main function
    '''
    # Main variables
    scan = args.scan
    num = args.num
    score = args.score
    spscore = args.spscore
    
    logging.info('Reading input file ' + str(args.infile))
    df = pd.read_feather(args.infile)
    
    logging.info('Removing duplicates')
    # grouped_df = df.groupby([scan])
    # for group in grouped_df:
    #     scan = group[1]
    #     scan = scan[scan.num != 1]
    #     scan.sort_values(by=[num], inplace=True, ascending=True)
    df.sort_values([scan, num, score, spscore], ascending=[True, True, False, False], inplace=True)
    df.drop_duplicates(subset=[scan], keep='first', inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    outfile = args.infile[:-8] + '_Unique.txt'
    logging.info('Writing output file ' + str(outfile))
    df.to_feather(outfile)
    
    logging.info('Done')
    

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='DuplicateRemover',
        epilog='''
        Example:
            python DuplicateRemover.py

        ''')
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-s', '--scan', required=True, help='ScanID column')
    parser.add_argument('-n', '--num', required=True, help='Rank column')
    parser.add_argument('-x', '--score', required=True, help='Score column')
    parser.add_argument('-p', '--spscore', required=True, help='SpScore column')

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args() 
    
    if '*' in args.infile: # wildcard
        flist = glob.glob(args.infile)
        for f in flist:
            args.infile = f
            # logging debug level. By default, info level
            log_file = args.infile[:-8] + '_log.txt'
            log_file_debug = args.infile[:-8] + '_log_debug.txt'
            # Logging debug level. By default, info level
            log_file = args.infile[:-8] + '_log.txt'
            log_file_debug = args.infile[:-8] + '_log_debug.txt'
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
            logging.info('start script: '+'{0}'.format(' '.join([x for x in sys.argv])))
            main(args)
        logging.info('end script')
    else:
        # logging debug level. By default, info level
        log_file = args.infile[:-8] + '_log.txt'
        log_file_debug = args.infile[:-8] + '_log_debug.txt'
        # Logging debug level. By default, info level
        log_file = args.infile[:-8] + '_log.txt'
        log_file_debug = args.infile[:-8] + '_log_debug.txt'
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
        logging.info('start script: '+'{0}'.format(' '.join([x for x in sys.argv])))
        main(args)
        logging.info('end script')