"""
Module for data handling and related things.
"""
import sys
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')
import re

import pandas as pd
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from root_numpy import array2root

from utils.misc_helpers import make_iterable

def check_branch_available(tree, branch, nowarn=False):
    """
    Check if a branch with the passed name is available in the TTree

    Args:
        tree (ROOT.TTree): tree for which the check should be performed
        branch (str): branch name
        nowarn (Boolean, optional): Suppress the warning log output

    Returns:
        bool: True if branch is in tree, else False
    """
    logging.debug('Checking if {} is available in {}'
                  .format(branch, tree.GetName()))

    all_branches = [b.GetName() for b in tree.GetListOfBranches()]
    if branch in all_branches:
        return True
    log_func = logging.info if nowarn else logging.warning
    log_func('Could not find branch \'{}\' in TTree \'{}\''
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
             DataFrame that pass the selection, unless all elements pass the
             selection, than the original dataframe will be returned
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

    if np.sum(sum_selection) == dataframe.shape[0]:
        logging.debug('Sum of selections (after possible negation) selects all '
                      'elements from passed DataFrame.')
        return dataframe

    # NOTE: since this indexing uses an array of bools this will always return a
    # copy
    return dataframe[sum_selection]


def get_treename(filename):
    """
    Try to get the ONLY tree in the passed root file

    Args:
        filename (str): Filename of the root file

    Returns:
        treename (str) or None: If there is only one tree in the root file, than
            the name of that tree will be returned, if there is more than one
            or no TTree at all, then None will be returned
    """
    rfl = r.TFile.Open(filename)
    contents = {k.GetName(): rfl.Get(k.GetName()) for k in rfl.GetListOfKeys()}
    trees = [n for n, o in contents.iteritems() if o.InheritsFrom('TTree')]

    if len(trees) == 1:
        return trees[0]

    if len(trees) > 1:
        logging.warning('Found more than one TTrees in {}: {}'
                        .format(filename, trees))
    if len(trees) == 0:
        logging.warning('Found no TTrees in {}'.format(filename))

    return None


def list_obj(rfile, obj_type='', filter_str=None):
    """
    Get a list (i.e. generator) of the names of all the objects in a root files

    Args:
        rfile (ROOT.TFile): root file
        obj_type (str, optional): Name of a root class from which the objects
            have to inherit. If no string is passed than all classes are valid.
        filter_str (str, regex, or None, optional): if not None, the names have
            to match this to be considered
    Returns:
        generator: A generator yielding all the names of the objects in the file
            that satisfy the criteria stated by the arguments
    """
    if filter_str is None:
        match_f = lambda x: True
    else:
        filter_rgx = re.compile(filter_str)
        match_f = lambda x: re.search(filter_rgx, x)

    inherits = lambda k: rfile.Get(k.GetName()).InheritsFrom(obj_type)
    if not obj_type:
        inherits = lambda k: True

    return (k.GetName() for k in rfile.GetListOfKeys()
            if match_f(k.GetName()) and inherits(k))


def common_obj(files, obj_type='', filter_str=None):
    """
    Get a list of the names of all objects that are common to all root files.

    Args:
        files (list of ROOT.TFile): root files
        obj_type (str, optional): Name of a root class from which the objects
            have to inherit. If no string is passed than all classes are valid.
        filter_str (str, regex or None, optional): If not None, the names have
            to match this to be considered

    Returns:
        list: A list with all names that are common to all files after
            considering the additional constraints that are posed by the input
            arguments.
    """
    file_obj = [list_obj(f, obj_type, filter_str) for f in files]
    common_obj = set(file_obj[0])
    for obj in file_obj[1:]:
        common_obj.intersection_update(obj)
    return list(common_obj)
