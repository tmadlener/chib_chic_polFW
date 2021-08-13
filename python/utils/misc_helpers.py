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
from decorator import decorator, decorate

import logging
logger = logging.getLogger()

try:
    basestring
except NameError:
    basestring = str

try:
    xrange
except:
    xrange = range

def cond_mkdir(path):
    """
    Conditionally make the directory with the passed path if it doesn't already
    exist. Raises an exception if something goes wrong (which is not an already
    existing folder)

    Implementation following: http://stackoverflow.com/a/14364249/3604607
    """
    try:
        logger.debug('trying to make directory: \'{}\''.format(path))
        if path:
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
            replacement. Longer matches get precedent, for patterns with equal
            length they are used as they appear in the list.
        reverse (bool): Reverse the ordering of the repl_pair list entries (i.e.
            when applying the function a second time, with reverse set to True
            the original string will be returned)

    Returns:
        str: string with all sub strings in the first element of the repl_pairs
            replaced with the accompanying second element
    """
    # swap the patterns and replacement here if we want to reverse the list
    if reverse:
        repl_pairs = [(rep, sym) for sym, rep in repl_pairs]

    # sort according to length of replacement pattern
    repl_pairs = sorted(repl_pairs, key=lambda p: len(p[0]), reverse=True)

    # fill a dict for easier replacing of the patterns using a regex
    subst = {p[0]: p[1] for p in repl_pairs}
    rgx = re.compile('|'.join(map(lambda p: re.escape(p[0]), repl_pairs)))
    return rgx.sub(lambda match: subst[match.group(0)], string)


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
    return ''.join(choice(ascii_letters + digits) for _ in range(length))


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

    logger.debug('Trying to find full path for subexpression: {}'
                  .format(subpath))
    trg_rgx = r'(HLT_)?(Dimuon)?(\d{1,2})_?(Jpsi|Upsilon)?_?(Barrel)?_?(Seagulls)?'
    trg_m = re.search(trg_rgx, subpath)
    if trg_m:
        logger.debug('Got a regex match, Now checking if we can uniquely '
                      'identify the path')

        match_ids = [(i, m) for i, m in enumerate(trg_m.groups())
                     if m is not None]
        for path, path_ids in unique_path_ids:
            if all([ids in match_ids for ids in path_ids]):
                return path

    logger.warning('Could not uniquely identify a trigger path from the sub '
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


def _get_var(dfr, get_var, np_func=None):
    """
    Get the variable from the dataframe depending on if is a one argument
    function, a string or already a numpy array
    """
    if np_func is None and isinstance(get_var, tuple):
        if len(get_var) == 2:
            get_var, np_func = get_var # unpack the tuple

    if np_func is not None:
        if isinstance(np_func, basestring):
            np_func = getattr(np, np_func)

    if hasattr(get_var, '__call__'):
        var = get_var(dfr)
    elif isinstance(get_var, basestring) or \
       (isinstance(get_var, list) and
        all(isinstance(v, basestring) for v in get_var)):
        var = dfr.loc[:, get_var]
    else:
        var = get_var

    return var if np_func is None else np_func(var)


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
    bin_var = np.sort(_get_var(dfr, get_var))
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


def get_bin_means(dfr, get_var, bins, selection=None, weights=None):
    """
    Get the the mean value in all bins from the data

    Args:
        dfr (pandas.DataFrame): DataFrame containing the data
        get_var (function): Function taking a DataFrame as only argument and
            returning one value for every row in the DataFrame. The return
            value of this will be used to calculate the mean in each bin
        bins (list of tuples): list of tuples containing the bin borders for
            each bin
        selection (numpy.array or list of functions, optional): Selection that
            is valid in a call to apply_selections
        weights(numpy.array, optional): If not none, the weighted average will
            be returned

    Returns:
        list: list of mean values for each bin

    See also:
        apply_selections, numpy.average
    """
    from utils.data_handling import apply_selections
    sel_dfr = apply_selections(dfr, selection)

    var = _get_var(sel_dfr, get_var)
    means = []
    for low, high in bins:
        bin_sel = (var > low) & (var < high)
        var_bin = var[bin_sel]
        w_bin = weights[bin_sel] if weights is not None else None
        # NOTE: here we use the fact that np.sum(None) does't evaluate to 0
        if np.sum(w_bin) != 0 and np.sum(bin_sel) != 0:
            # Have to do this to avoid np.average raising a ZeroDivisionError
            means.append(np.average(var_bin, weights=w_bin))
        else:
            means.append(np.nan)

    return means


def deprecated_soon(replacement=None):
    """
    Decorator to inform the user that this function will soon be deprecated
    offering the possibility to point out a replacement.
    """
    @decorator
    def _deprecated(func, *args, **kwargs):
        """
        Inform the user that this function will soon be deprecated.
        """
        message = '\'{}\' will soon be deprecated.'.format(func.__name__)
        if replacement is not None:
            message += ' You should start to use the replacement \'{}\''.format(replacement)
        logger.warn(message)

        return func(*args, **kwargs)

    return _deprecated


@deprecated_soon('select_bin')
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
    var = _get_var(dfr, bin_var)
    return (var > bin_low) & (var < bin_up)


def select_bin(var, low, high):
    """
    Create a function that selects events within a given bin for a variable

    Args:
        var (str or function): The variable for which the bin should be defined.
            This can either be a function that returns the variable from the
            DataFrame, or a string that can be parsed by parse_func_val
        low (float): lower edge of bin
        high (float): upper edge of bin

    Returns:
        function: A function that takes a DataFrame as only argument and returns
            an array of boolean values that can then be used to retrieve only
            the selected elements from the DataFrame, e.g. in apply_selections

    See also:
        parse_func_val, apply_selections
    """
    class BinSelection(object):
        """
        Internal helper class that captures everything that is necessary to get
        the selection function
        """
        def __init__(self, var_exp, b_lo, b_hi):
            # have to check here whether we have a function or a string to parse
            if hasattr(var_exp, '__call__'):
                self.var = (var_exp, None)
            else:
                self.var = parse_func_var(var_exp)

            self.low = b_lo
            self.high = b_hi

            # allow it to be used in collect_requirements environment
            self.requires = [self.var[0]]

        def __call__(self, dfr):
            var = _get_var(dfr, *self.var)
            return (var > self.low) & (var < self.high)


    return BinSelection(var, low, high)


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


def combine_cuts(cuts, comb=' && '):
    """
    Combine cut strings to one so that they can be safely used in the ROOT cut
    string interface.

    Args:
        cuts (list of str): Cuts to be combined
        comb (str, defaults): How the cuts should be combined. Defaults to
            returning the AND of all cuts
    """
    safe_cuts = ["".join(["(", c, ")"]) for c in cuts]
    return comb.join(safe_cuts)


def get_vals_from_rwbuffer(rw_buffer, n_points):
    """
    Get the first n_points from a read-write buffer

    TGraphAsymmErrors returns read-write buffers without an appropriate size.
    Thus it is necessary to read the points using the knowledge of how many
    values to read.

    Throws and IndexError exception if the rw_buffer is empty (i.e. null-buffer)

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


def get_np_from_tmatrix(tmatrix):
    """
    Convert a TMatrixT<T> to a numpy array

    Args:
        tmatrix (TMatrixT<T>): The ROOT matrix to convert into a numpy array

    Returns:
        numpy.array: a 2D numpy array with the values of the tmatrix
    """
    n_cols = tmatrix.GetNcols()
    n_rows = tmatrix.GetNrows()

    raw_vals = get_vals_from_rwbuffer(tmatrix.GetMatrixArray(), n_cols * n_rows)
    return np.reshape(raw_vals, (n_rows, n_cols))


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
            logger.warning('Caught KeyError')

    return try_catch_access(*args)


def flatten(iterable):
    """
    Flatten any list or iterable into a 1D-list. Keep dictionaries and strings
    intact.
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
        if hasattr(elem, "__iter__") and not isinstance(elem, dict):
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
        logger.warning('Cannot find a common binning without overlaps. All '
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


def get_bin_centers(binning):
    """
    Get the centers of the passed binning
    """
    n_bins = len(binning) - 1
    centers = [0.5 * (binning[i] + binning[i + 1]) for i in xrange(n_bins)]
    return np.array(centers)


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


def unique_w_key(elements, key):
    """
    Get unique elements of an iterable where uniqueness is defined by the key

    Args:
        elements (iterable): Elements that are possibly non-unique
        key (callable): Function that takes one element of the iterable as
            argument and returns a hashable object which will be used to compare
            different elements

    Returns:
        list of unique elements. If there are multiple elements that evaluate to
            the same value according to the key function the first element of
            the elements will be used
    """
    seen = set()
    # set.add returns None so that the first expression in the comprehension has
    # two effects: adding it to the seen dict, and also using it in the list
    # comprehension, but only if it is not already in seen
    # see: https://stackoverflow.com/a/10024750
    return [seen.add(key(e)) or e for e in elements if key(e) not in seen]


def parse_binning(binning_str):
    """
    Parse the variable binning and return the desired binning.

    Args:
        binning_str (str): The string that should be parsed

    Returns:
         numpy.array: The bin edge values that were obtained from parsing the
             string

    Valid strings are:
    * comma separated list of floating point numbers
    * 'BEGIN:END:DELTA', where the interval from BEGIN to END will be
      divided into intervals of length DELTA. END is not necessarily included
      as it uses numpy.arange
    * 'BEGIN:END,NSTEPS', where the interval from BEGIN to END will be
      divided into NSTEPS intervals of equal length. This will include END as a
      bin edge as it uses numpy.linspace

    See also:
        numpy.linspace, numpy.arange
    """
    # check if one of the two special strings is present and handle them
    # appropriately
    flt_rgx = r'(-?\d+(?:\.\d*)?)' # match a floating point number

    # Match the linspace case
    rmatch =  re.match(flt_rgx + ':' + flt_rgx + r',(\d+)$', binning_str)
    if rmatch:
        return np.linspace(float(rmatch.group(1)), float(rmatch.group(2)),
                           int(rmatch.group(3)))

    # match the arange case
    rmatch = re.match(r':'.join([flt_rgx]*3), binning_str)
    if rmatch:
        return np.arange(float(rmatch.group(1)), float(rmatch.group(2)),
                         float(rmatch.group(3)))

    # check if we have a comma separated list of values
    if ',' in binning_str and not ':' in binning_str:
        try:
            return np.array([float(v) for v in binning_str.split(',') if v])
        except ValueError as parse_err:
            logger.error('Could not parse binning from \'{}\' because \'{}\''
                          .format(binning_str, parse_err.message))
            return np.array([])

    # if we fall through here then we got a case that we don't handle
    logger.error('Cannot handle \'{}\' in trying to parse a binning'
                  .format(binning_str))

    return np.array([])


def float_rgx(char_separated=False):
    rgx = r'([-+]?\d*\.?\d+)'
    if char_separated:
        return rgx.replace(r'\.', 'p')
    return rgx


def parse_func_var(expr):
    """
    Parse the functional form of an expression

    Args:
        expr (string): The expression that holds the variable and the function
            on that variable in the form func(var). The function has to be
            callable on numpy.arrays (NOTE: This is not checked here, but might
            fail later)

    Returns:
        tuple: The name of the variable as string and the numpy function. If no
            no function is called, then the function return will be None. If the
            parsing fails None will be returned
    """
    func_rgx = r'^(\w+)\((\w+)\)$'
    match = re.match(func_rgx, expr)
    if match:
        var = match.group(2)
        func = match.group(1)
        if hasattr(np, func):
            return var, getattr(np, func)
        else:
            logger.error('Could not find a numpy function named \'{}\' that '
                          'was parsed from \'{}\''.format(func, expr))

    var_rgx = r'^\w+$'
    match = re.match(var_rgx, expr)
    if match:
        return expr, None


def parse_sel_expr(expr):
    """
    Parse a selection expression of the form XY < funcExpr < YZ (where one of
    the two values XY or YZ can also be omitted and funcExpr is an expression
    that can be parsed by parse_func_var):

    Args:
        expr (string): Expression that defines a selection and should be parsed

    Returns:
        list or False: list of expression functions that can be passed to
            apply_selections and select what is specified in the passed expr.
            Returning False in case of failure should ensure that at least the
            apply_selections call fails when this is used inline
    """
    parts = [s.strip() for s in expr.split('<')]

    flt_rgx = re.compile(float_rgx())
    if len(parts) == 3:
        vmatch1, vmatch2 = flt_rgx.match(parts[0]), flt_rgx.match(parts[2])
        if vmatch1 and vmatch2:
            return select_bin(parts[1], float(parts[0]), float(parts[2]))

    if len(parts) == 2:
        # Need to find out whether the value or the expression comes first
        vmatch1, vmatch2 = flt_rgx.match(parts[0]), flt_rgx.match(parts[1])
        if vmatch1 and not vmatch2:
            func = parse_func_var(parts[1])
            if func is not None:
                return lambda d: _get_var(d, *func) > float(parts[0])
        if not vmatch1 and vmatch2:
            func = parse_func_var(parts[0])
            if func is not None:
                return lambda d: _get_var(d, *func) < float(parts[1])

    # If we are still here then we could not parse the expression
    logger.error('Could not parse expression \'{}\' to extract a selection '
                  'function'.format(expr))
    return False


def is_divisable(dividend, divisor):
    """
    Check if the dividend is (evenly) divisable by the passed divisor and return
    the result if it is

    Args:
        dividend (int)
        divisor int)

    Returns:
        ratio (int) or None: If number is evenly divisable by factor then the
            ratio is returned otherwise None
    """
    ratio = dividend / divisor
    if ratio * divisor == dividend:
        return ratio
    logger.warning('{} is not a divisor of {}'.format(divisor, dividend))


def memoize(func):
    """
    Rather simple memoize implementation supporting arbitrary signature function
    calls. Since it simply appends a caching dict to the function, it is neither
    optimal in performance, nor does it cover all edge cases.

    Implementation is shamelessly stolen from:
    https://decorator.readthedocs.io/en/latest/tests.documentation.html#the-solution
    """
    def _memoize(func, *args, **kwargs):
        if kwargs:
            key = args, frozenset(kwargs.items())
        else:
            key = args
        cache = func.cache
        if key not in cache:
            cache[key] = func(*args, **kwargs)
        return cache[key]

    func.cache = {}
    return decorate(func, _memoize)


# for formatting, the absolute value of the decimal exponent above which the
# decimal formatting will be used
MAX_ABS_EXP = 2

# The maximum number of the digits after the decimal point
MAX_N_DIGITS = 4

def fmt_float(number, use_exp=None):
    """
    Format a floating point number such that it is (more or less) nicely written
    using decimal exponents if necessary

    Args:
        number (float):
        use_exp (int, optional): Use this exponent to format the number

    Returns:
        str: The formatted number
    """
    if number == 0:
        return '0'

    if number < 0:
        exp = np.floor(np.log10(-number))
    else:
        exp = np.floor(np.log10(number))

    if use_exp is not None:
        if not isinstance(use_exp, int):
            logger.warning('Rounding use_exp={} to nearest int'.format(use_exp))
        # float to avoid "Integers to negative integer powers are not allowed"
        exp = np.round(use_exp).astype(float)

    if np.abs(exp) > MAX_ABS_EXP:
        fmt_str = r'{{:.{}f}} \cdot 10^{{{{{{:.0f}}}}}}'.format(MAX_N_DIGITS - 1)
        return fmt_str.format(number * 10**(-exp), exp)

    fmt_str = '{{:.{}f}}'.format(MAX_N_DIGITS)
    return fmt_str.format(number)


def quantile(vals, quant, weights=None):
    """
    Get the quantile from the values using weights

    Calculates the quantiles of the possibly weighted distribution of vals, by
    building the empirical CDF and then doing an interpolation on the inverse to
    get the values corresponding to a given quantile.

    Args:
        vals (np.array): Values for which the quantiles should be calculated
        quant (array like): The quantiles (between 0 and 1) which should be
            determined
        weights (np.array, optional): Array of weights of the same length as
            vals. Each entry will be used to weight the corresponding value in
            vals.

    Returns:
         quantile(s): float or nd.array. The quantiles of the possibly weighted
             distribution of vals.
    """
    sort_idx = np.argsort(vals)
    sort_vals = vals[sort_idx]
    if weights is not None:
        sort_weights = weights[sort_idx]
    else:
        sort_weights = np.ones_like(sort_vals)

    part_sum = np.cumsum(sort_weights)

    # Calcualte the empirical distribution function by dividing the partial sums
    # By the total sum
    ecdf = part_sum / part_sum[-1]

    # In principle it would be enough now to find the element for which the ecdf
    # exceeds the desired quantile(s), but it is easier to do this for more
    # than one quantile by simply interpolating on the "inverse" ecdf to find
    # the values at which the quantiles are.
    return np.interp(quant, ecdf, sort_vals)
