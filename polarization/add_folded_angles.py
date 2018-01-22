#!/usr/bin/env python
"""
Add the folded angular variables to the root file
"""

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import numpy as np

from utils.data_handling import store_dataframe, get_dataframe


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
    phi_fold[quad_4] = -phi_fold[quad_4]

    quad_2 = (phi >= 90) & (phi < 180)
    phi_fold[quad_2] = 180 - phi_fold[quad_2]
    costh_fold[quad_2] = -costh_fold[quad_2]

    quad_3 = (phi >= -180) & (phi < -90)
    phi_fold[quad_3] = 180 + phi_fold[quad_3]
    costh_fold[quad_3] = -costh_fold[quad_3]

    return costh_fold, phi_fold


def add_folded_frame(df, frame, costhn='costh', phin='phi'):
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
    """
    logging.debug('Adding folded variables for frame {}'.format(frame))
    try:
        # get the values from the dataframe as np.arrays
        phi = df['_'.join([phin, frame])].values
        costh = df['_'.join([costhn, frame])].values
        costh_f, phi_f = calc_folded_vars(costh, phi)

        df['_'.join([phin, frame, 'fold'])] = phi_f
        df['_'.join([costhn, frame, 'fold'])] = costh_f
    except KeyError:
        logging.warning('Could not add folded variables for frame {}, because '
                        'unfolded angular variables were not available'
                        .format(frame))


def main(args):
    """Main"""
    for infile in args.inputfiles:
        logging.info('Processing file {}'.format(infile))
        var_names = [('costh', 'phi')]
        if args.genlevel:
            logging.info('Also adding generator level folding')
            var_names.append(('gen_costh', 'gen_phi'))
        df = get_dataframe(infile, args.treename)
        for var_pair in var_names:
            for frame in args.frames.split(','):
                add_folded_frame(df, frame, *var_pair)

        store_dataframe(df, infile, args.treename)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script that adds the folded '
                                     'angular distribuion variables to the root'
                                     ' file in a separate branch')
    parser.add_argument('inputfiles', nargs='+', help='Files to process')
    parser.add_argument('-t', '--treename', default='chic_tuple', type=str,
                        help='name of the tree in which the original variables are '
                        '(used for storing output file).')
    parser.add_argument('-f', '--frames', type=str, default='HX,PX,CS',
                        help='reference frames for which to add the folded variables'
                        ' (comma separated list of two char abbreviations)')
    parser.add_argument('-g', '--genlevel', default=False, action='store_true',
                        help='add the folding also for the generator level angles')

    clargs = parser.parse_args()
    main(clargs)
