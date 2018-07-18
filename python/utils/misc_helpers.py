"""
Module for small miscellaneous helper functions that are intended to handle
things concerning os or python built-ins.
"""

import os
import re
import numpy as np
from random import choice
from string import ascii_letters, digits
from collections import Iterable
from decorator import decorator

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


def make_iterable(possibly_iterable):
    """
    Make the passed element iterable if it is a single element or return it
    unchanged if it already is.

    Args:
        possibly_iterable: Any object or built-in or anything that is
            already iterable

    Returns:
        list or possibly_iterable: If input is not iterable this returns the
            passed object wrapped inside a list, otherwise it simply returns
            the passed in object unchanged.
    """
    if not isinstance(possibly_iterable, Iterable):
        return [possibly_iterable]
    return possibly_iterable


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
        ('HLT_Dimuon10_Upsilon_Barrel_Seagulls', ((2, '10'), (3, 'Upsilon'))),
        ('HLT_Dimuon12_Upsilon_eta1p5', ((2, '12'),)),
        ('HLT_Dimuon10_Jpsi_Barrel', ((2, '10'), (3, 'Jpsi'))),
        ('HLT_Dimuon8_Upsilon_Barrel', ((2, '8'), (3, 'Upsilon'))),
        ('HLT_Dimuon8_Jpsi', ((2, '8'), (3, 'Jpsi'))),
        ('HLT_Dimuon16_Jpsi', ((2, '16'),)),
        ('HLT_Dimuon20_Jpsi', ((2, '20'),)),
        ('HLT_Dimuon25_Jpsi', ((2, '25'),)),
    )

    logging.debug('Trying to find full path for subexpression: {}'
                  .format(subpath))
    trg_rgx = r'(HLT_)?(Dimuon)?(\d{1,2})_?(Jpsi|Upsilon)?_?(Barrel)?_?(Seagulls)?'
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
    n_overflow = bin_var.shape[0] % n_bins

    if n_overflow == 0:
        # bin borders are the elements at the positions i * dbin - 1
        edge_idcs = [0] + [i * dbin - 1 for i in xrange(1, n_bins)] + [-1]
    else: # distribute the overflowing events evenly over the whole range
        edge_idcs = [0]
        for i in xrange(1, n_bins):
            if i < n_overflow:
                edge_idcs.append(i * dbin + i)
            else:
                edge_idcs.append(i * dbin + n_overflow)

        edge_idcs.append(-1)

    binb = [bin_var[i] for i in edge_idcs]

    return zip(binb[:-1], binb[1:])


def get_costh_binning(dfr, n_bins, full_range=False, selection=None):
    """
    Get an equi-populated binning in abs(costh_HX)

    Args:
        dfr (pandas.DataFrame): DataFrame containing all the data that should be
            considered for the binning
        n_bins (int): Number of bins
        full_range (bool): Make the last bin span the full range up to 1 (True)
            or stop at the maximum observed value (False)
        selection (numpy.array, optional): Selection array that can be used in a
            DataFrame indexing to select certain events, defaults to None, where
            all entries are used

    Returns:
        list: list of tuples with the bin borders for all the bins, where the
            first bin starts at 0 and the last one ends at 1, regardless of the
            exact values for the bin borders
    """
    from utils.data_handling import apply_selections
    sel_dfr = apply_selections(dfr, selection)

    abs_costh = lambda d: np.abs(d.costh_HX)
    binning = get_equi_pop_bins(sel_dfr, abs_costh, n_bins)
    # replace the lowest and highest bin border
    binning[0] = (0, binning[0][1])
    if full_range:
        binning[-1] = (binning[-1][0], 1)
    return binning


def get_bin_means(dfr, get_var, bins, selection=None):
    """
    Get the the mean value in all bins from the data

    Args:
        dfr (pandas.DataFrame): DataFrame containing the data
        get_var (function): Function taking a DataFrame as only argument and
            returning one value for every row in the DataFrame. The return
            value of this will be used to calculate the mean in each bin
        bins (list of tuples): list of tuples containing the bin borders for
            each bin
        selection (numpy.array, optional): Selection array that can be used in a
            DataFrame indexing to select certain events, defaults to None, where
            all entries are used

    Returns:
        list: list of mean values for each bin
    """
    from utils.data_handling import apply_selections
    sel_dfr = apply_selections(dfr, selection)

    var = get_var(sel_dfr)
    means = []
    for low, high in bins:
        var_bin = var[(var > low) & (var < high)]
        means.append(np.mean(var_bin))

    return means


def get_bin_cut_df(dfr, bin_var, bin_low, bin_up):
    """
    Get the binning selection for a variable in a format suitable for dataframes

    Args:
        dfr (pandas.DataFrame): The dataframe from which the selection is done
        bin_var (str or function): The variable for which the binning should be
            done. If it is a function it has to take the dfr as single argument
            and return a value for each row in the dataframe, such that it can
            be used to index into the dataframe
        bin_low (float): lower edge of bin
        bin_up (float): upper edge of bin

    Returns:
        pandas.Series: A series that can be used to index into the originally
            passed DataFrame to select only events fullfilling the binn criteria
    """
    if hasattr(bin_var, '__call__'):
        # if bin_Var is a function call it on the dataframe to decide the cut
        cut = (bin_var(dfr) > bin_low) & (bin_var(dfr) < bin_up)
    else:
        cut = (dfr[bin_var] > bin_low) & (dfr[bin_var] < bin_up)
    return cut


def get_bin_cut_root(bin_var, bin_low, bin_up):
    """
    Get the binning cut for a variable in a format root understands

    Args:
        bin_var (str): Name of the branch that is used for binning
        bin_low (float): lower edge of bin
        bin_up (float): upper edge of bin

    Returns:
        str: The selection string that specifies the cut for the bin
    """
    return '{0} > {1} && {0} < {2}'.format(bin_var, bin_low, bin_up)


def get_vals_from_rwbuffer(rw_buffer, n_points):
    """
    Get the first n_points from a read-write buffer

    TGraphAsymmErrors returns read-write buffers without an appropriate size.
    Thus it is necessary to read the points using the knowledge of how many
    values to read.

    Args:
        rw_buffer (read-write buffer): read-write buffer returned by ROOT
            functions returning a Double_t* in C++
        n_points (int): Number of values to read from buffer

    Returns:
        numpy.array: numpy array containing the values
    """
    vals = [0] * n_points
    for i in xrange(n_points):
        vals[i] = rw_buffer[i]

    return np.array(vals)


@decorator
def log_key_error(func, *args):
    """
    Decorator for logging a KeyError, where such an exception is not critical to
    the operation of the program
    """
    def try_catch_access(*args):
        """The closure holding the try-catch block"""
        try:
            return func(*args)
        except KeyError:
            logging.warning('Caught KeyError')

    return try_catch_access(*args)


def flatten(iterable):
    """
    Flatten any list or iterable into a 1D-list.
    Taken from here: http://stackoverflow.com/a/2158532/3604607

    Args:
        iterable (iterable): Possibly nested list or iterable that should be
            converted into a flat (1D) list

    Returns:
        generator: The generator that yields all the elements of the passed in
            list as a flat iterable
    """
    for elem in iterable:
        # only checking for the __iter__ attribute here should allow to enter
        # sub lists but not tear apart strings since they don't have it
        if hasattr(elem, "__iter__"):
            for sub in flatten(elem):
                yield sub
        else:
            yield elem


def chunks(iterable, chunk_size):
    """
    Split the iterable into chunks and return the list of chunks

    Args:
        iterable (iterable): Flat list of n elements that will be split into
            n / chunk_size chunks
        chunk_size (int): Number of elements per chunk

    Returns:
        generator: The generator that yields the list of chunks, where each
            chunk contains chunk_size elements of the input and the last chunk
            contains the remainder of elements
    """
    for ichunk in xrange(0, len(iterable), chunk_size):
        yield iterable[ichunk:ichunk + chunk_size]


def find_common_binning(bin_borders):
    """
    Try to find a common binning in all the passed bin borders.

    Will first find all possible bins and then check if there is an overlap in
    them by comparing the lower boundaries with the upper boundaries of the
    preceding bin. If no overlap is found the possible bins will be returned.

    Args:
        bin_borders (list of np.arrays): list of 2D np.arrays where each array
            holds the bins (low, high) as columns for the different bins (rows).

    Return:
        np.array or None: 2D array with the bins (low, high) as columns for the
            different bins (rows). None if no common binning can be found.
    """
    # First get all possible bins from all the passed bin_borders
    # Then see if they overlap
    get_bins = lambda bb: [(bb[i, 0], bb[i, 1]) for i in xrange(len(bb))]
    poss_bins = [
        b for b in set(get_bins(bin_borders[0])).union(
            *(get_bins(bb) for bb in bin_borders))
    ]
    poss_bins.sort(key=lambda b: b[0]) # sort bins according to low bin boundary
    poss_bins = np.array(poss_bins) # for easier slicing

    # If any of the upper boundaries of the preceding bin are above the lower
    # boundary of the current bin then they overlap
    low_bounds = poss_bins[1:, 0]
    high_bounds = poss_bins[:-1, 1]
    if np.all(low_bounds <= high_bounds):
        return poss_bins
    else:
        logging.warning('Cannot find a common binning without overlaps. All '
                        'possible bins are: {}'.format(poss_bins))
        return None


def get_bin(binning, value):
    """
    Get the bin index of the given value in the binning
    """
    for ibin in xrange(len(binning) - 1):
        if value > binning[ibin] and value <= binning[ibin + 1]:
            return ibin

    return -1


def get_bin_edges(bins):
    """
    Get the bin edges from a list of bins (low, high)
    """
    uniq_bin_bord = sorted({b for b in flatten(bins)})
    return np.array(uniq_bin_bord)


def longest_match(string, poss_strings):
    """
    Find the longest match of a string in a set of possibly matching strings

    Args:
        string (str): String that should be matched
        poss_strings (list): List of strings that should be matched against
            string

    Returns:
        str: The element of poss_strings that has the longest match in string
    """
    match_strs = [ps for ps in poss_strings if ps in string]
    match_strs.sort(key=lambda x: len(x))

    return match_strs[-1]
