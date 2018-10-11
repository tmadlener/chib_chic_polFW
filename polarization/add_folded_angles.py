#!/usr/bin/env python
"""
Add the folded angular variables to the root file
"""

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import numpy as np

from itertools import product

from utils.data_handling import get_dataframe, add_branch, get_treename
from utils.misc_helpers import flatten

def calc_folded_vars(costh, phi):
    """
    Calculate the folded variales from the unfolded ones.

    The folded values are calculated as follows:
    if phi >= -90 and phi < 0:
        phi = -1 * phi
    elif phi >= 90 and phi < 180::
        phi = 180 - phi
        costh = -1 * costh
    elif phi >= -180 and phi < -90:
        phi = 180 + phi
        costh = -1 * costh

    Args:
        costh (numpy.array): unfolded costh values
        phi (numpy.array): unfolded phi values

    Returns:
        (costh_f, phi_f): tuple of numpy.arrays containing the folded
        values
    """
    phi_fold = np.copy(phi)
    costh_fold = np.copy(costh)

    quad_4 = (phi >= -90) & (phi < 0)
    phi_fold[quad_4] *= -1

    quad_2 = (phi >= 90) & (phi < 180)
    phi_fold[quad_2] = 180 - phi_fold[quad_2]
    costh_fold[quad_2] *= -1

    quad_3 = (phi >= -180) & (phi < -90)
    phi_fold[quad_3] += 180
    costh_fold[quad_3] *= -1

    return costh_fold, phi_fold


def get_folded_frame(df, frame, costhn='costh', phin='phi'):
    """
    Add the folded variables to the dataframe for the passed frame

    Args:
        df (pandas.DataFrame): The DataFrame holding all the data for which the
            folded angles should be added
        frame (str): The reference frame for which the folded angles should be
            added.
        costhn, phin (str, optional): The base names of the costh and phi
            variables in the DataFrame (Default to 'costh' and 'phi' so that
            e.g. the complete column names for the HX frame are generated to be
            'costh_HX' and 'phi_HX')

    Returns:
        (costh_f, phi_f): tuple of numpy.arrays containing the folded
        values
    """
    logging.debug('Adding folded variables for frame {}'.format(frame))
    # get the values from the dataframe as np.arrays
    phi = df['_'.join([phin, frame])].values
    costh = df['_'.join([costhn, frame])].values
    return calc_folded_vars(costh, phi)


def main(args):
    """Main"""
    # In order to not have to load the whole file first determine the names of
    # the branches that are necessary
    var_names = [('costh', 'phi')]
    if args.genlevel:
        logging.info('Also adding generator level folding')
        var_names.append(('gen_costh', 'gen_phi'))
    frames = args.frames.split(',')

    load_variables = ['_'.join(p) for p in product(flatten(var_names), frames)]

    for infile in args.inputfiles:
        logging.info('Processing file {}'.format(infile))
        if not args.treename:
            treename = get_treename(infile)
        else:
            treename = args.treename

        df = get_dataframe(infile, treename, columns=load_variables)

        for var_pair in var_names:
            for frame in frames:
                costh_f, phi_f = get_folded_frame(df, frame, *var_pair)
                add_branch(costh_f, '_'.join([var_pair[0], frame, 'fold']),
                           infile, treename)
                add_branch(phi_f, '_'.join([var_pair[1], frame, 'fold']),
                           infile, treename)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script that adds the folded '
                                     'angular distribuion variables to the root'
                                     ' file in a separate branch')
    parser.add_argument('inputfiles', nargs='+', help='Files to process')
    parser.add_argument('-t', '--treename', default='', type=str,
                        help='name of the tree in which the original variables are '
                        'stored. (Only necessary if more than one TTree is in '
                        ' the input file')
    parser.add_argument('-f', '--frames', type=str, default='HX,PX,CS',
                        help='reference frames for which to add the folded variables'
                        ' (comma separated list of two char abbreviations)')
    parser.add_argument('-g', '--genlevel', default=False, action='store_true',
                        help='add the folding also for the generator level angles')

    clargs = parser.parse_args()
    main(clargs)
