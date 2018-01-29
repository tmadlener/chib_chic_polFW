"""
Module for small miscellaneous helper functions that are intended to handle
things concerning os or python built-ins.
"""

import os
import re
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
    if reverse:
        for rep, sym in repl_pairs:
            ret_str = ret_str.replace(sym, rep)
    else:
        for sym, rep in repl_pairs:
            ret_str = ret_str.replace(sym, rep)

    return ret_str


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
