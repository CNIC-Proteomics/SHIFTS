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
import configparser
import glob
import logging
import pandas as pd
import sys
import os


def main(args):
    '''
    Main function
    '''
    # Main variables
    scan = config._sections['DuplicateRemover']['scancolumn']
    num = config._sections['DuplicateRemover']['rankcolumn']
    score = config._sections['DuplicateRemover']['scorecolumn']
    spscore = config._sections['DuplicateRemover']['spscorecolumn']
    
    logging.info('Reading input file2 ' + str(args.infile))
    df = pd.read_feather(args.infile)
    
    logging.info('Removing duplicates')
    df.sort_values([scan, num, score, spscore], ascending=[True, True, False, False], inplace=True)
    df.drop_duplicates(subset=[scan], keep='first', inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    outfile = args.infile[:-8] + '_Unique.feather'
    logging.info('Writing output file ' + str(outfile))
    df.to_feather(outfile)
    
    logging.info('Done')
    

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='DuplicateRemover v2',
        epilog='''
        Example:
            python DuplicateRemover.py

        ''')

    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")

    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')

    # these will overwrite the config if specified
    parser.add_argument('-s', '--scan', default=None, help='ScanID column')
    parser.add_argument('-n', '--num', default=None, help='Rank column')
    parser.add_argument('-x', '--score', default=None, help='Score column')
    parser.add_argument('-p', '--spscore', default=None, help='SpScore column')

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    # parse config
    mass_config = configparser.ConfigParser(inline_comment_prefixes='#')
    mass_config.read(os.path.join(os.path.dirname(__file__), "config/MassMod.ini"))
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.scan is not None:
        config.set('DuplicateRemover', 'scancolumn', str(args.scan))
        config.set('Logging', 'create_ini', '1')
    if args.num is not None:
        config.set('DuplicateRemover', 'scancolumn', str(args.num))
        config.set('Logging', 'create_ini', '1')
    if args.score is not None:
        config.set('DuplicateRemover', 'scorecolumn', str(args.score))
        config.set('Logging', 'create_ini', '1')
    if args.spscore is not None:
        config.set('DuplicateRemover', 'spscorecolumn', str(args.spscore))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)


    
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