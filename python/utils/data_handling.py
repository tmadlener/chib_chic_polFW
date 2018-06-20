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

from root_numpy import fill_hist

from utils.hist_utils import set_hist_opts
from utils.misc_helpers import create_random_str, make_iterable

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

    logging.info('Writing resulting DataFrame to: {}'.format(outfile))
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


def get_dataframe(infile, treename=None, columns=None):
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
        columns (str or sequence of str, optional): Only read the specified
            branches into the DataFrame (only used when reading from a root file)

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
            return read_root(infile, key=treename, columns=columns)
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
        selections (list or function): List of functions taking the DataFrame as
            single argument and returning a list of booleans (with the same
            number) of rows as the DataFrame, where the elements with True will
            be selected
        negate (Boolean): Instead of returning all events fulfilling the
            selection return all events not fulfilling the selection

    Returns:
         pandas.DataFrame: New DataFrame with only the elements of the passed
             DataFrame that pass the selection
    """
    if selections is None:
        return dataframe

    sum_selection = np.ones(dataframe.shape[0], dtype=bool)
    for sel in make_iterable(selections):
        sum_selection &= sel(dataframe)

    if negate:
        sum_selection = np.invert(sum_selection)

    # NOTE: since this indexing uses an array of bools this will always return a
    # copy
    return dataframe[sum_selection]


def create_histogram(var, hist_sett, **kwargs):
    """
    Create a ROOT histogram from the passed variable(s)

    Args:
        var (np.array): Array with maximum of 3 columns containing the variables
            to plot.
        hist_set (tuple): Histogram settings, that are directly unpacked into
            the constructor of the ROOT histogram

    Keyword Args:
        name (str, optional): Name to be used for the histogram
        weights (np.array, optional): weight array with the same number of
             events as the var array. Each entry corresponds to the weight of
             the event
        {x,y,z}_axis (str): axis labels to be set for the histogram

    Returns:
         ROOT.TH{1,2,3}D: The histogram with the dimension corresponding to the
             number of columns of var
    """
    name = kwargs.pop('name', '')
    if not name:
        name = create_random_str()
    # use the number of dimensions from the var to determine which sort of
    # histogram to use
    ndim = var.shape
    if len(ndim) == 1:
        ndim = 1
    else:
        ndim = ndim[1]

    if ndim > 3 or ndim < 0:
        logging.error('Dimension of histogram is {}. Cannot create histogram'
                      .format(ndim))
        raise TypeError('Invalid number of dimensions in create_histograms')

    hist_type = 'TH{}D'.format(ndim)
    try:
        hist = getattr(r, hist_type)(name, '', *hist_sett)
    except TypeError as exc:
        logging.error('Could not construct TH{}D with passed hist_sett: {}'
                      .format(ndim, hist_sett))
        raise exc

    set_hist_opts(hist)

    # set axis labels
    xax, yax, zax = (kwargs.pop(a, '') for a in ['x_axis', 'y_axis', 'z_axis'])
    if xax:
        hist.SetXTitle(xax)
    if yax:
        hist.SetYTitle(yax)
    if zax:
        hist.SetZTitle(zax)

    fill_hist(hist, var, weights=kwargs.pop('weights', None))

    return hist
