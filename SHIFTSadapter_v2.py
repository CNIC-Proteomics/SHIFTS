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
import pyarrow
import sys


def main(ifile, ofile):
    '''
    Main function
    '''
    # obtain the first line
    logging.info('Giving the input file ' + str(ifile))
    with open(ifile) as f:
        first_line = f.readline().strip().split('\t')
    
    # read the data depending on the type of search engine
    if 'CometVersion' in first_line[0]:
        logging.info('Reading the "comet" data file')
        df = pd.read_csv(ifile, sep='\t', skiprows=1, float_precision='high', low_memory=False, index_col=False)
    else:
        logging.info('Reading the "msfragger" data file')
        df = pd.read_csv(ifile, sep='\t', float_precision='high', low_memory=False, index_col=False)
    
    # add the file name without extension into 'Raw' column
    df["Raw"] = '.'.join(os.path.basename(Path(ifile)).split(".")[:-1])
    
    # write file
    outfile = ofile[:-4] + '_SHIFTS.feather'
    logging.info('Writing output file ' + str(outfile))
    # df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    df = df.reset_index(drop=True)
    df.to_feather(outfile)
        

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='SHIFTSadapter v2',
        epilog='''
        Example:
            python SHIFTSadapter_v2.py

        ''')
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-o', '--outdir', help='Output dir')

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    # getting input parameters
    ifile = args.infile

    # get the output file
    # if output directory is not defined, get the folder from given file
    outdir = args.outdir if args.outdir else os.path.dirname( ifile )
    ofile = os.path.join( outdir, os.path.basename(ifile) )

    if '*' in ifile: # wildcard
        flist = glob.glob(ifile)
        for f in flist:
            # create ofile
            of = os.path.join( outdir, os.path.basename(f) )
            # logging debug level. By default, info level
            log_file = ofile[:-4] + '_log.txt'
            log_file_debug = ofile[:-4] + '_log_debug.txt'
            # Logging debug level. By default, info level
            log_file = ofile[:-4] + '_log.txt'
            log_file_debug = ofile[:-4] + '_log_debug.txt'
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
            main(f,of)
        logging.info('end script')
    else:
        # logging debug level. By default, info level
        log_file = ofile[:-4] + '_log.txt'
        log_file_debug = ofile[:-4] + '_log_debug.txt'
        # Logging debug level. By default, info level
        log_file = ofile[:-4] + '_log.txt'
        log_file_debug = ofile[:-4] + '_log_debug.txt'
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
        main(ifile,ofile)
        logging.info('end script')
        