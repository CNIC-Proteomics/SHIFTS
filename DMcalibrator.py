#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "andrea.laguillo@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# import modules
import os
import sys
from pathlib import Path
import configparser
import argparse
import logging
import math
import pandas as pd
from scipy.special import erfinv
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

# os.chdir(r"C:\Users\Andrea\Desktop\SHIFTS-4")

###################
# Local functions #
###################

# Calibrate mass separately for each raw file.

def readInfile(infile, scorecolumn, mzcolumn, zcolumn, seqcolumn, proteincolumn):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high', low_memory=False) # TODO: option for header/no header
    #df = pd.read_csv(infile, sep="\t", float_precision='high')
    # Cleanup
    df = df[df[scorecolumn].notna()]
    df[scorecolumn] = pd.to_numeric(df[scorecolumn])
    df = df[df[mzcolumn].notna()]
    df[mzcolumn] = pd.to_numeric(df[mzcolumn])
    df = df[df[zcolumn].notna()]
    df[zcolumn] = pd.to_numeric(df[zcolumn])
    df = df[df[seqcolumn].notna()]
    #df = df[df[config._sections['Input']['dmcolumn']].notna()]
    #df[config._sections['Input']['dmcolumn']] = pd.to_numeric(df[config._sections['Input']['dmcolumn']])
    df = df[df[proteincolumn].notna()]
    return df

def getTheoMZ(df, mzcolumn, zcolumn, seqcolumn):
    '''    
    Calculate theoretical MZ using the PSM sequence.
    '''
    AAs = dict(mass_config._sections['Aminoacids'])
    MODs = dict(mass_config._sections['Fixed Modifications'])
    m_proton = mass_config.getfloat('Masses', 'm_proton')
    m_hydrogen = mass_config.getfloat('Masses', 'm_hydrogen')
    m_oxygen = mass_config.getfloat('Masses', 'm_oxygen')
    if 'theo_mz' not in df:
        df.insert(df.columns.get_loc(mzcolumn)+1, 'theo_mz', np.nan)
    if 'theo_mh' not in df:
        df.insert(df.columns.get_loc('theo_mz'), 'theo_mh', np.nan)
    
    def _PSMtoMZ(sequence, charge):
        total_aas = 2*m_hydrogen + m_oxygen
        total_aas += charge*m_proton
        total_aas += float(MODs['nt']) + float(MODs['ct'])
        for aa in sequence:
            if aa.lower() in AAs:
                total_aas += float(AAs[aa.lower()])
            #else: # aminoacid not in list (ask for user input?)
                # TODO
            if aa.lower() in MODs:
                total_aas += float(MODs[aa.lower()])
        MH = total_aas - (charge-1)*m_proton
        #MZ = (total_aas + int(charge)*m_proton) / int(charge)
        MZ = total_aas / int(charge)
        return MZ, MH
    
    df['theo_mz'] = df.apply(lambda x: _PSMtoMZ(x[seqcolumn], x[zcolumn])[0], axis = 1)
    df['theo_mh'] = df.apply(lambda x: _PSMtoMZ(x[seqcolumn], x[zcolumn])[1], axis = 1)
    #df['theo_mh'] = df.apply(lambda x: (x['theo_mz'] * x[zcolumn]) - (m_proton * (x[zcolumn]-1)), axis = 1)
    return df

def getErrors(df, mzcolumn, calibrated):
    '''    
    Calculate absolute (in m/z) and relative (in ppm) errors.
    '''
    if calibrated:
        abs_error = 'cal_dm_mz' #cal_abs_error
        rel_error = 'cal_ppm'
        i = 2
    else:
        abs_error = 'abs_error'
        rel_error = 'ppm'
        i = 1
        
    if abs_error not in df:
        df.insert(df.columns.get_loc('theo_mz')+i, abs_error, np.nan)
        
    if calibrated:
        if rel_error not in df:
            df.insert(df.columns.get_loc(abs_error)+1, rel_error, np.nan)
        df[abs_error] = df['cal_exp_mz'] - df['theo_mz']
        #df[rel_error] = (df[abs_error] / df['theo_mz']) * 1e6
        df[rel_error] = df[abs_error] / df[mzcolumn] * 1e6
    else:
       # if 'exp_mh' not in df:
            #df.insert(df.columns.get_loc(mzcolumn)+1, rel_error, np.nan)
        #df['exp_mh'] = df[mzcolumn] * df[config._sections['Input']['zcolumn']] + df[config._sections['Input']['zcolumn']] * mass_config.getfloat('Masses', 'm_proton')
        df[abs_error] = df[mzcolumn] - df['theo_mz']
        #df[rel_error] = df[abs_error] / df['theo_mz'] * 1e6
    return df

def filterPeptides(df, scoremin, ppmmax, scorecolumn, chargecolumn, mzcolumn,
                   seqcolumn, proteincolumn, abscolumn, decoyprefix):
    '''    
    Filter and keep target peptides that match Xcorrmin and PPMmax conditions.
    This high-quality subpopulation will be used for calibration.
    '''
    
    # def _correctXcorr(charge, xcorr, length):
    #     if charge < 3:
    #         cxcorr = math.log10(xcorr) / math.log10(2*length)
    #     else:
    #         cxcorr = math.log10(xcorr/1.22) / math.log10(2*length)
    #     return cxcorr
    
    #keep targets
    df_filtered = df[~df[proteincolumn]
                     .str.startswith(decoyprefix)]
    #keep score > scoremin
    df_filtered = df_filtered[df_filtered[scorecolumn]
                              >=scoremin]
    #keep abs_error <= ppmmax
    df_filtered['abs_error_ppm'] = df_filtered[abscolumn]/df_filtered[mzcolumn] * 1e6
    df_filtered = df_filtered[df_filtered['abs_error_ppm']
                              <=ppmmax]
    df_filtered = df_filtered[df_filtered['abs_error_ppm']
                              >=-ppmmax]
    df_filtered = df_filtered.drop('abs_error_ppm', axis = 1)
    logging.info("Number of PSMs before filtering: " + str(df.shape[0]))
    logging.info("Number of PSMs after filtering: " + str(df_filtered.shape[0]))
    return df_filtered

def getSysError(df_filtered, mzcolumn, calibrated):
    '''
    Calculate systematic error and average PPM error.
    '''
    if calibrated:
        abs_error = 'cal_dm_mz' #cal_abs_error
    else:
        abs_error = 'abs_error'
        
    sys_error = df_filtered[abs_error].median()
    alpha = (df_filtered[abs_error]/df_filtered[mzcolumn]).median()
    
    if calibrated:
        phi = math.sqrt(2) * erfinv(0.5)
        mad = df_filtered['cal_ppm'].mad()
        avg_ppm_error = (mad / phi) 
        logging.info("Systematic error after calibration: " + "{:.4e}".format(sys_error))
        logging.info("Alpha after calibration: " + "{:.4e}".format(alpha))
        logging.info("StdDevMAD_ppm: " + "{:.4e}".format(avg_ppm_error))
        return sys_error, alpha, avg_ppm_error
    else:
        logging.info("Systematic error: " + "{:.4e}".format(sys_error))
        logging.info("Alpha: " + "{:.4e}".format(alpha))
        return sys_error, alpha

def rawCorrection(df, mzcolumn, alpha):
    '''
    Correct exp_mz values from infile using the systematic error.
    '''
    if 'cal_exp_mz' not in df:
        df.insert(df.columns.get_loc(mzcolumn)+1, 'cal_exp_mz', np.nan)
        df.insert(df.columns.get_loc('cal_exp_mz')+1, 'cal_exp_mh', np.nan)
    #if 'exp_mh_cal' not in df:
        #df.insert(df.columns.get_loc('cal_exp_mz')+1, 'exp_mh_cal', np.nan)
    
    def _correct(exp_mz, abs_error, alpha):
        cal_exp_mz = exp_mz * (1  - alpha)
        return cal_exp_mz
    
    #df['cal_exp_mz'] = df[config._sections['Input']['mzcolumn']] - sys_error
    df['cal_exp_mz'] = df.apply(lambda x: _correct(x[mzcolumn], x['abs_error'], alpha), axis = 1)
    df['cal_exp_mh'] = df.apply(lambda x: (x['cal_exp_mz'] * x[config._sections['DMcalibrator']['zcolumn']]) - ((x[config._sections['DMcalibrator']['zcolumn']]-1) * mass_config.getfloat('Masses', 'm_proton')), axis = 1)
    return df

def getDMcal(df, mzcolumn, calmzcolumn, zcolumn):
    '''
    Calculate calibrated DM values.
    '''
    # Before calibration
    if 'dm_mz' not in df:
        df.insert(df.columns.get_loc(mzcolumn)+1,
                  'dm_mz',
                  np.nan)
    df['dm_mz'] = df[mzcolumn] - df['theo_mz']
    if 'dm_mh' not in df:
        df.insert(df.columns.get_loc(mzcolumn)+1,
                  'dm_mh',
                  np.nan)
    #df['dm_mh'] = df['dm_mz'] * df[zcolumn] - (df[zcolumn] * mass_config.getfloat('Masses', 'm_proton'))
    #df['dm_mh'] = df.apply(lambda x: (x['dm_mz'] * x[zcolumn] - (x[zcolumn] * m_proton)) if (x['dm_mz'] >= 0) else (x['dm_mz'] * x[zcolumn] + (x[zcolumn] * m_proton)), axis = 1)
    df['dm_mh'] = df['dm_mz'] * df[zcolumn]
    # After calibration
    if 'cal_dm_mz' not in df:
        df.insert(df.columns.get_loc(calmzcolumn)+1,
                  'cal_dm_mz',
                  np.nan)
    df['cal_dm_mz'] = (df[calmzcolumn] - df['theo_mz'])
    if 'cal_dm_mh' not in df:
        df.insert(df.columns.get_loc(calmzcolumn)+1,
                  'cal_dm_mh',
                  np.nan)
    #df['cal_dm_mh'] = (df['cal_dm_mz'] * df[zcolumn]) - (df[zcolumn] * mass_config.getfloat('Masses', 'm_proton'))
    #df['cal_dm_mh'] = df.apply(lambda x: (x['cal_dm_mz'] * x[zcolumn] - (x[zcolumn] * m_proton)) if (x['cal_dm_mz'] >= 0) else (x['cal_dm_mz'] * x[zcolumn] + (x[zcolumn] * m_proton)), axis = 1)
    df['cal_dm_mh'] = (df['cal_dm_mz'] * df[zcolumn])
    return df

def labelTargetDecoy(df, proteincolumn, decoyprefix):
    '''
    Label targets and decoys according to protein ID column.
    '''
    if 'Label' not in df:
        df.insert(df.columns.get_loc(proteincolumn)+1, 'Label', np.nan)
    df['Label'] = df.apply(lambda x: 'Decoy' if (x[proteincolumn][0:len(decoyprefix)]==decoyprefix) else 'Target', axis = 1)
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
        sys.exit("Unrecognized sequence format in '" + str(config._sections['DMcalibrator']['seqdmcolumn']) + "' column!")
    return formatseq

#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    # Main variables 
    score_min = float(config._sections['DMcalibrator']['score_min'])
    ppm_max = float(config._sections['DMcalibrator']['ppm_max'])
    scorecolumn = config._sections['DMcalibrator']['scorecolumn']
    zcolumn = config._sections['DMcalibrator']['zcolumn']
    mzcolumn = config._sections['DMcalibrator']['mzcolumn']
    seqcolumn = config._sections['DMcalibrator']['seqcolumn']
    seqdmcolumn = config._sections['DMcalibrator']['seqdmcolumn']
    #dmcolumn = config._sections['DMcalibrator']['dmcolumn']
    proteincolumn = config._sections['DMcalibrator']['proteincolumn']
    decoyprefix = config._sections['DMcalibrator']['decoyprefix']
    abscolumn = 'abs_error'
    calabscolumn = 'cal_dm_mh'
    calmzcolumn = 'cal_exp_mz'
    calseqcolumn = config._sections['DMcalibrator']['calseqcolumn']
    decimal_places = int(config._sections['General']['decimal_places'])
    
    log_str = "Calibrating file: " + str(Path(args.infile))
    logging.info(log_str)
    # Read infile
    df = readInfile(Path(args.infile),
                    scorecolumn,
                    mzcolumn,
                    zcolumn,
                    seqcolumn,
                    proteincolumn)
    # Label targets and decoys
    df = labelTargetDecoy(df, proteincolumn, decoyprefix)
    # Calculate theoretical MZ
    df = getTheoMZ(df, mzcolumn, zcolumn, seqcolumn)
    # Calculate errors
    df = getErrors(df, mzcolumn, 0)
    # Filter identifications
    logging.info("Filtering by score_min = " + str(score_min))
    logging.info("Filtering by ppm_max = " + str(ppm_max))
    df_filtered = filterPeptides(df,
                                 score_min,
                                 ppm_max,
                                 scorecolumn,
                                 zcolumn,
                                 mzcolumn,
                                 seqcolumn,
                                 proteincolumn,
                                 abscolumn,
                                 decoyprefix)
    if len(df_filtered) == 0:
        logging.info("ERROR: Cannot calibrate without a filtered dataset. Please lower the score or ppm thresholds and try again.")
        sys.exit()

    # Use filtered set to calculate systematic error
    sys_error, alpha = getSysError(df_filtered, mzcolumn, 0)
    # Use systematic error to correct infile
    df = rawCorrection(df, mzcolumn, alpha)
    # Recalculate systematic error using calibrated masses
    df = getErrors(df, calmzcolumn, 1)
    df_filtered = filterPeptides(df,
                                 score_min,
                                 ppm_max,
                                 scorecolumn,
                                 zcolumn,
                                 mzcolumn,
                                 seqcolumn,
                                 proteincolumn,
                                 abscolumn,
                                 decoyprefix)
    cal_sys_error, cal_alpha, avg_ppm_error = getSysError(df_filtered, mzcolumn, 1)
    # Calculate DMCal 
    df = getDMcal(df, mzcolumn, calmzcolumn, zcolumn)
    # Make calseqcolumn
    df.insert(df.columns.get_loc(seqdmcolumn)+1, calseqcolumn, np.nan)
    #df[calseqcolumn] = df.apply(lambda x: x[seqdmcolumn].split('[')[0] + '[' + str(round(x['cal_dm_mh'], decimal_places)) + ']' + x[seqdmcolumn].split(']')[1], axis = 1)
    df[calseqcolumn] = df.apply(lambda x: format_seq(x[seqdmcolumn], x['cal_dm_mh'], decimal_places), axis = 1)
    #Write to txt file
    logging.info("Writing output file...")
    outfile = args.infile[:-4] + '_calibrated.txt'
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    logging.info("Calibration finished")

    
if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='DMcalibrator',
        epilog='''
        Example:
            python DMcalibrator.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    # these will overwrite the config if specified
    parser.add_argument('-s', '--scoremin', default=None, help='Minimum score')
    parser.add_argument('-p', '--ppmmax', default=None, help='Maximum PPM error')
    parser.add_argument('-sc', '--scorecolumn', default=None, help='Name of the column containing the score')
    parser.add_argument('-zc', '--chargecolumn', default=None, help='Name of the column containing the charge')
    parser.add_argument('-mc', '--mzcolumn', default=None, help='Name of the column containing the experimental m/z')
    parser.add_argument('-se', '--seqcolumn', default=None, help='Name of the column containing the sequence')
    #parser.add_argument('-dm', '--dmcolumn', default=None, help='Name of the column containing the deltamass')

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    mass_config = configparser.ConfigParser(inline_comment_prefixes='#')
    mass_config.read(os.path.join(os.path.dirname(__file__), "config/MassMod.ini"))
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.scoremin is not None:
        config.set('DMcalibrator', 'score_min', str(args.scoremin))
        config.set('Logging', 'create_ini', '1')
    if args.ppmmax is not None:
        config.set('DMcalibrator', 'ppm_max', str(args.ppmmax))
        config.set('Logging', 'create_ini', '1')
    if args.scorecolumn is not None:
        config.set('DMcalibrator', 'scorecolumn', str(args.scorecolumn))
        config.set('Logging', 'create_ini', '1')
    if args.mzcolumn is not None:
        config.set('DMcalibrator', 'mzcolumn', str(args.mzcolumn))
        config.set('Logging', 'create_ini', '1')
    if args.chargecolumn is not None:
        config.set('DMcalibrator', 'zcolumn', str(args.zcolumn))
        config.set('Logging', 'create_ini', '1')
    if args.seqcolumn is not None:
        config.set('DMcalibrator', 'seqcolumn', str(args.seqcolumn))
        config.set('Logging', 'create_ini', '1')
    #if args.dmcolumn is not None:
        #config.set('Input', 'dmcolumn', str(args.dmcolumn))
        #config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)
        

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_log_debug.txt'
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