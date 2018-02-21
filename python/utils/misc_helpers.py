"""
Module for small miscellaneous helper functions that are intended to handle
things concerning os or python built-ins.
"""

import os
import re
import numpy as np
from random import choice
from string import ascii_letters, digits

import logging
logging.basicConfig(level=logging.WARNING,
                    format='%(levelname)s - %(funcName)s: %(message)s')

def cond_mkdir(path):
    """
    Conditionally make the directory with the passed path if it doesn't already
    exist. Raises an exception if something goes wrong (which is not an already
    existing folder)

    Implementation following: http://stackoverflow.com/a/14364249/3604607
    """
    try:
        os.makedirs(path)
    except OSError as err:
        if not os.path.isdir(path):
            print('Caught error \'{}\' while trying to create directory \'{}\''
                  .format(err.strerror, path))
            raise


def cond_mkdir_file(filename):
    """
    Conditionally make the folder that should hold the passed filename, so that
    after a call to this the filename can be used to store something into it.

    The filename is removed (i.e. everything after the last "/") and the
    directory is then conditionally made.
    """
    path, _ = os.path.split(filename)
    cond_mkdir(path)


def replace_all(string, repl_pairs, reverse=False):
    """
    Replace all occurrences of different sub strings with other specified
    substrings.

    Args:
        string (str): input string
        repl_pairs (list of tuples): list containing tuples where the first
            element is the string to be replaced and the second element is the
            replacement
        reverse (bool): Reverse the ordering of the repl_pair list entries (i.e.
            when applying the function a second time, with reverse set to True
            the original string will be returned)

    Returns:
        str: string with all sub strings in the first element of the repl_pairs
            replaced with the accompanying second element
    """
    ret_str = string
    if reverse:
        for rep, sym in repl_pairs:
            ret_str = ret_str.replace(sym, rep)
    else:
        for sym, rep in repl_pairs:
            ret_str = ret_str.replace(sym, rep)

    return ret_str



def stringify(selection, reverse=False):
    """Get a string representation of a selection string

    Get a string representation without any special characters from a selection
    as it is used e.g. in TTree::Draw or in RooDataSet::reduce.

    Args:
        selection (str): Selection string
        reverse (bool, optional): If true, the original representation can be
            obtained (with stripped whitespace)

    Returns:
        str
    """
    repl_pairs = (
        ('>=', '_ge_'),
        ('<=', '_le_'),
        ('<', '_lt_'),
        ('>', '_gt_'),
        ('==', '_eq_'),
        ('!=', '_ne_'),
        ('&&', '_AND_'),
        ('||', '_OR_'),
        ('.', 'p'),
        ('&', '_BITAND_'),
        ('(', '_PAROP_'), # Rather ugly and hacky solution to handle parenthesis
        (')', '_PARCL_')  # without usage of regular expressions
    )

    ret_str = selection.replace(' ', '')
    return replace_all(ret_str, repl_pairs, reverse)


def create_random_str(length=32):
    """
    Create a random string containing upper and lower case ASCII characters and
    numbers

    Implementation following: https://stackoverflow.com/questions/2257441/
    Args:
        length (int, optional): length of the output string

    Returns:
        str: Non-cryptographically safe random string
    """
    return ''.join(choice(ascii_letters + digits) for _ in xrange(length))


def get_full_trigger(subpath):
    """
    Get the full name of the trigger path from a subexpression

    Args:
        subpath (str): A uniquely identifying subpath of a a trigger

    Returns:
        str: The full name of the trigger identified by subpath or empty string
            if a trigger cannot by uniquely identified
    """
    # define a map from paths to the identifiers that are at least necessary
    # to identify a given path. The indices are the indices of the group in the
    # trg_rgx
    unique_path_ids = (
        ('HLT_Dimuon20_Jpsi_Barrel_Seagulls', ((2, '20'), (4, 'Barrel'))),
        ('HLT_Dimuon10_Jpsi_Barrel', ((2, '10'),)),
        ('HLT_Dimuon8_Jpsi', ((2, '8'),)),
        ('HLT_Dimuon16_Jpsi', ((2, '16'),)),
        ('HLT_Dimuon20_Jpsi', ((2, '20'),)),
        ('HLT_Dimuon25_Jpsi', ((2, '25'),)),
    )

    logging.debug('Trying to find full path for subexpression: {}'
                  .format(subpath))
    trg_rgx = r'(HLT_)?(Dimuon)?(\d{1,2})(_Jpsi)?(_Barrel)?(_Seagulls)?'
    trg_m = re.search(trg_rgx, subpath)
    if trg_m:
        logging.debug('Got a regex match, Now checking if we can uniquely '
                      'identify the path')

        match_ids = [(i, m) for i, m in enumerate(trg_m.groups())
                     if m is not None]
        for path, path_ids in unique_path_ids:
            if all([ids in match_ids for ids in path_ids]):
                return path

    logging.warning('Could not uniquely identify a trigger path from the sub '
                    'expression: {}'.format(subpath))
    return ''


def get_storable_name(name, reverse=False):
    """
    Get a name that can be stored as a TBranch such that it can also be read
    afterwards

    E.g. '-' and '.' are characters that cannot appear in the branch name if
    the TTree::Draw() command should be usable with the branch

    Args:
        name (str): Original name of the branch
        reverse (bool): If passed a storable name, return the original one,
           making this the inverse of the function

    Returns:
        str: storeable name or original name, depending on argument of reverse
    """
    repl_pairs = (
        ('-', '_m_'),
        ('.', '_p_')
    )
    return replace_all(name, repl_pairs, reverse)


def get_equi_pop_bins(dfr, get_var, n_bins):
    """
    Get equi-populated bins for the data

    Args:
        dfr (pandas.DataFrame): DataFrame containing all the data that should
            be considered for the binning
        get_var (function): Function taking a DataFrame as only argument and
            returning one value for every row in the DataFrame. The return value
            of this function will be the variable in which the binning is done.
        n_bins (int): The number of bins

    Returns:
        list: List of tuples with two elements, where the [0] element is the
            lower bound of the bin and [1] the upper bound of the bin
    """
    bin_var = np.sort(get_var(dfr))
    dbin = bin_var.shape[0] / n_bins # number of events per bin
    # bin borders are the elements at the positions i * dbin - 1
    binb = [bin_var[0]] + [bin_var[i * dbin - 1] for i in xrange(1, n_bins)]
    binb += [bin_var[-1]] # the last bin has to include everything up to the max

    return zip(binb[:-1], binb[1:])


def get_costh_binning(dfr, n_bins, selection=None):
    """
    Get an equi-populated binning in abs(costh_HX)

    Args:
        dfr (pandas.DataFrame): DataFrame containing all the data that should be
            considered for the binning
        n_bins (int): Number of bins
        selection (numpy.array, optional): Selection array that can be used in a
            DataFrame indexing to select certain events, defaults to None, when
            all entries are used

    Returns:
        list: list of tuples with the bin borders for all the bins, where the
            first bin starts at 0 and the last one ends at 1, regardless of the
            exact values for the bin borders
    """
    if selection is not None:
        sel_dfr = dfr[selection]
    else:
        sel_dfr = dfr

    abs_costh = lambda d: np.abs(d.costh_HX)
    binning = get_equi_pop_bins(sel_dfr, abs_costh, n_bins)
    # replace the lowest and highest bin border
    binning[0] = (0, binning[0][1])
    binning[-1] = (binning[-1][0], 1)
    return binning
