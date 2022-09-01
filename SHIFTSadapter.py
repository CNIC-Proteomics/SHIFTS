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
import os
import pandas as pd
from pathlib import Path
import sys

def main(args):
    '''
    Main function
    '''
    # Main variables
    logging.info('Reading input file ' + str(args.infile))
    with open(args.infile) as f:
        first_line = f.readline().strip().split('\t')
    df = pd.read_csv(args.infile, sep='\t', skiprows=1, float_precision='high', low_memory=False, index_col=False)
    
    logging.info('Search Engine: ' + first_line[0])
    logging.info('Raw: ' + first_line[1])
    logging.info('Date: ' + first_line[2])
    logging.info('Database: ' + first_line[3])
    
    df["Raw"] = os.path.basename(Path(args.infile))
    
    outfile = args.infile[:-4] + '_SHIFTS.txt'
    logging.info('Writing output file ' + str(outfile))
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    
    logging.info('Done')
    

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='SHIFTSadapter',
        epilog='''
        Example:
            python SHIFTSadapter.py

        ''')
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()
    
    if '*' in args.infile: # wildcard
        flist = glob.glob(args.infile)
        for f in flist:
            args.infile = f
            # logging debug level. By default, info level
            log_file = args.infile[:-4] + '_log.txt'
            log_file_debug = args.infile[:-4] + '_log_debug.txt'
            # Logging debug level. By default, info level
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
            logging.info('start script: '+'{0}'.format(' '.join([x for x in sys.argv])))
            main(args)
        logging.info('end script')
    else:
        # logging debug level. By default, info level
        log_file = args.infile[:-4] + '_log.txt'
        log_file_debug = args.infile[:-4] + '_log_debug.txt'
        # Logging debug level. By default, info level
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
        logging.info('start script: '+'{0}'.format(' '.join([x for x in sys.argv])))
        main(args)
        logging.info('end script')
        