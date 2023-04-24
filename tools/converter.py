#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 11:39:39 2023

@author: alaguillog
"""

# import modules
import argparse
import glob
import logging
import pandas as pd
from pathlib import Path

#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    infile = Path(args.infile)
    ext = infile.suffix.lower()
    # Convert
    if ext == '.feather':
        logging.info('Reading feather file...')
        df = pd.read_feather(args.infile)
        logging.info('Writing tab-separated text file...')
        # infile.rename(infile.with_suffix('.tsv'))
        outfile = args.infile.split('.')
        outfile[-1] = 'tsv'
        outfile = '.'.join(outfile)
        df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    elif ext in ['.txt', '.tsv']:
        logging.info('Reading tab-separated text file...')
        df = pd.read_csv(args.infile, sep="\t")
        logging.info('Writing feather file...')
        # infile.rename(infile.with_suffix('.feather'))
        outfile = args.infile.split('.')
        outfile[-1] = 'feather'
        outfile = '.'.join(outfile)
        df.to_feather(outfile)
    elif ext in ['.csv']:
        logging.info('Reading comma-separated text file...')
        df = pd.read_csv(args.infile, sep=",")
        logging.info('Writing feather file...')
        outfile = args.infile.split('.')
        outfile[-1] = 'feather'
        outfile = '.'.join(outfile)
        # infile.rename(infile.with_suffix('.feather'))
        df.to_feather(outfile)
    else:
        logging.error('Input file format "' + ext +'" not recognized. Skipping...')

if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser(
        description='DMcalibrator',
        epilog='''
        Example:
            python DMcalibrator.py

        ''')
    parser.add_argument('-i', '--infile', required=True, help='Path to input file(s)')
    args = parser.parse_args()

    if '*' in args.infile: # wildcard
        flist = glob.glob(args.infile)
        for f in flist:
            args.infile = f
            logging.basicConfig(level=logging.INFO,
                                format='%(asctime)s - %(levelname)s - %(message)s',
                                datefmt='%m/%d/%Y %I:%M:%S %p',
                                handlers=[logging.StreamHandler()])
            # start main function
            logging.info('start script')
            main(args)
        logging.info('end script')
    else:
        logging.info('start script')
        main(args)
        logging.info('end script')
