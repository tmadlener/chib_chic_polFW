"""
Module for data handling and related things.
"""
import sys
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import pandas as pd
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from root_numpy import array2root

from utils.misc_helpers import make_iterable

def check_branch_available(tree, branch):
    """
    Check if a branch with the passed name is available in the TTree

    Args:
        tree (ROOT.TTree): tree for which the check should be performed
        branch (str): branch name

    Returns:
        bool: True if branch is in tree, else False
    """
    logging.debug('Checking if {} is available in {}'
                  .format(branch, tree.GetName()))

    all_branches = [b.GetName() for b in tree.GetListOfBranches()]
    if branch in all_branches:
        return True
    logging.warning('Could not find branch \'{}\' in TTree \'{}\''
                    .format(branch, tree.GetName()))
    return False


def store_dataframe(dfr, outfile, tname='chi2_values', **kwargs):
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
            DataFrame if stored to a root file

    Keyword Args:
        Forwarded to root_pandas.to_root

    See Also: root_pandas.to_root
    """
    logging.debug('Storing DataFrame to {}'.format(outfile))
    if not outfile.endswith('.pkl') and not outfile.endswith('.root'):
        logging.warning('Output file doesnot have .root or .pkl format. '
                        'Creating a .pkl file instead')
        logging.debug('Output filename before substitution: {}'.format(outfile))
        import re
        outfile = re.sub(r'(.*\.)(\w*)$', r'\1pkl', outfile)
        logging.debug('Output filename after substitution: {}'.format(outfile))

    logging.info('Writing resulting DataFrame to: {}'.format(outfile))
    # if .root is requested check if root_pandas is here, otherwise go to .pkl
    if outfile.endswith('.root'):
        try:
            from root_pandas import to_root
            # current version of to_root doesn't support the store_index argument
            to_root(dfr, outfile, tname, **kwargs# , store_index=False
            )
        except ImportError:
            logging.warning('Output to .root file was requested, but root_pandas'
                            ' was not found. Creating a .pkl file instead')
            outfile = outfile.replace('.pkl', '.root')

    if outfile.endswith('.pkl'):
        dfr.to_pickle(outfile)


def add_branch(arr, bname, rfile, tname):
    """
    Add the passed array to an existing TTree in an existing TFile

    Args:
        arr (numpy.array): 1D numpy array that will be stored under a new branch
        bname (str): Branch name for the values
        rfile (str): Filename to which the new branch should be added
        tname (str): Name of the TTree to which the values should be added
    """
    arr = np.array(arr, dtype=[(bname, np.find_common_type(arr, []))])
    array2root(arr, rfile, treename=tname, mode='update')


def get_dataframe(infile, treename=None, **kwargs):
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
    Keyword Args:
         Forwarded to root_pandas.read_root

    See also: root_pandas.read_root

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
            return read_root(infile, key=treename, **kwargs)
        except ImportError:
            # log and bail out
            logging.error('Requested to read DataFrame from {}, but could not '
                          'import root_pandas'.format(infile))
    sys.exit(1)


def apply_selections(dataframe, selections, negate=False):
    """
    Apply all selections and return the reduced dataframe.

    Args:
        dataframe (pandas.DataFrame): The data to which the selections should be
            applied
        selections (list of functions, function or numpy.ndarray): List of
            functions taking the DataFrame as single argument and returning a
            list of booleans (with the same number) of rows as the DataFrame,
            where the elements with True will be selected or a selection array
            that can be used to index into a DataFrame and select certain events
        negate (Boolean): Instead of returning all events fulfilling the
            selection return all events not fulfilling the selection

    Returns:
         pandas.DataFrame: New DataFrame with only the elements of the passed
             DataFrame that pass the selection
    """
    if selections is None:
        return dataframe

    if isinstance(selections, np.ndarray) or isinstance(selections, pd.Series):
        sum_selection = selections
    else:
        selections = make_iterable(selections)
        # Check if all selections are actually functions. If not sipmly log as
        # this will fail in the next few lines anyway
        if not all(callable(f) for f in selections):
            logging.error('Passed selections are not all functions and also not'
                          ' an array of boolean indices')
        sum_selection = np.ones(dataframe.shape[0], dtype=bool)
        for sel in selections:
            sum_selection &= sel(dataframe)

    if negate:
        sum_selection = np.invert(sum_selection)

    # NOTE: since this indexing uses an array of bools this will always return a
    # copy
    return dataframe[sum_selection]
