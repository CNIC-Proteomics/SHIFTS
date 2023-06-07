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
import logging
import pandas as pd
import sys

def main(args):
    '''
    Main function
    '''
    # Main variables
    logging.info('Reading input file')
    with open(args.infile) as f:
        first_line = f.readline().strip().split('\t')
    df = pd.read_csv(args.infile, sep='\t', skiprows=0, float_precision='high', low_memory=False, index_col=False)
    
    logging.info('Cleaning up filename')
    df['Filename'] = df.apply(lambda x: str(x['Filename']).replace(str(args.remove), ''), axis = 1)
    
    logging.info('Generating ScanID')
    df['ScanID'] = str(df['Filename']) + '-' + str(df['scan']) + '-' + str(df['charge'])
    df['ScanID'] = df.apply(lambda x: str(x['Filename']) + '-' +
                                      str(x['scan']) + '-' +
                                      str(x['charge']), axis = 1)
    
    logging.info('Writing output file')
    outfile = args.infile[:-4] + '_ScanID.txt'
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    
    logging.info('Done')
    

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='ScanIDgenerator',
        epilog='''
        Example:
            python ScanIDgenerator.py

        ''')
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-r', '--remove', required=False, default='_SHIFTS_Unique_calibrated', help='String to remove from Filename')

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args() 

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.StreamHandler()])

    # start main function
    logging.info('start script: '+'{0}'.format(' '.join([x for x in sys.argv])))
    main(args)
    logging.info('end script')