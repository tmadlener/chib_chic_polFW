"""
Module for data handling and related things.
"""
import sys
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import pandas as pd

def store_dataframe(dfr, outfile, tname='chi2_values'):
    """
    Store the dataframe either into a pkl file or into a root file via
    root_pandas.

    Args:
        dfr (pandas.DataFrame): The dataframe that should be stored
        outfile (str): The filename to which the DataFrame should be stored.
            If this ends with .pkl, a pkl file will be created, if it ends on
            .root a root file will be created (if root_pandas is available),
            Otherwise a .pkl file will be created with .root replaced with .pkl
        tname (str, optional): Name of the TTree to be used for storing the
            DataFrame is stored to a root file
    """
    logging.debug('Storing DataFrame to {}'.format(outfile))
    if not outfile.endswith('.pkl') and not outfile.endswith('.root'):
        logging.warning('Output file doesnot have .root or .pkl format. '
                        'Creating a .pkl file instead')
        logging.debug('Output filename before substitution: {}'.format(outfile))
        import re
        outfile = re.sub(r'(.*\.)(\w*)$', r'\1pkl', outfile)
        logging.debug('Output filename after substitution: {}'.format(outfile))

    logging.info('Writing resulting DataFramet to: {}'.format(outfile))
    # if .root is requested check if root_pandas is here, otherwise go to .pkl
    if outfile.endswith('.root'):
        try:
            from root_pandas import to_root
            # current version of to_root doesn't support the store_index argument
            to_root(dfr, outfile, tname, mode='w'# , store_index=False
            )
        except ImportError:
            logging.warning('Output to .root file was requested, but root_pandas'
                            ' was not found. Creating a .pkl file instead')
            outfile = outfile.replace('.pkl', '.root')

    if outfile.endswith('.pkl'):
        dfr.to_pickle(outfile)


def get_dataframe(infile, treename=None):
    """
    Get the dataframe from the input file.

    Args:
        infile (str): Name of the inputfile from which the dataframe should be
            read. Must either be a pkl or a root file. Which format will be read
            depends entirely on the ending of the filename.
        treename (str, optional): The TTree in the TFile that should be read.
            Since it is possible to store multiple trees in one file it can be
            necessary to specify which on to read. Option is only used for reads
            from .root files.

    Returns:
        pandas.DataFrame: The dataframe read from the file.
    """
    logging.debug('Getting DataFrame from {}'.format(infile))
    if not infile.endswith('.pkl') and not infile.endswith('.root'):
        logging.error('Infile does not have a parseable format: {}'
                      ' Valid formats are .root and .pkl'.format(infile))

    if infile.endswith('.pkl'):
        return pd.read_pickle(infile)
    if infile.endswith('.root'):
        try:
            from root_pandas import read_root
            return read_root(infile, key=treename)
        except ImportError:
            # log and bail out
            logging.error('Requested to read DataFrame from {}, but could not '
                          'import root_pandas'.format(infile))
    sys.exit(1)
