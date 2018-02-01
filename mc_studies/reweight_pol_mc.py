#!/usr/bin/env python
"""
Script for reweighting MC to any desired polarization
"""

import sys

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.data_handling import get_dataframe, store_dataframe
from utils.pol_utils import ang_dist_lth
from utils.misc_helpers import get_storable_name

def calc_weights(costh, lth):
    """
    Calculate the weights to reweight the sample to have the desired
    lambda_theta in the frame in which costh is defined.

    Args:
        costh (float, or numpy.array)
        lth (float): lambda_theta of the desired polarization, assuming that it
            is specified in the same frame as costh

    Returns:
        float or numpy.array: The weight that has to be applied to the event to
            get a polarized sample, with lth. Return type is dependent on the
            type of costh
    """
    logging.debug('Adding weights for lambda_theta = {}'.format(lth))
    return ang_dist_lth(costh, lth)


def main(args):
    """Main"""
    weight_name = get_storable_name('wPol_lth_{:.2f}'.format(args.lambda_theta))
    mc_frame = get_dataframe(args.inputfile)
    if weight_name in mc_frame.columns:
        logging.info('weights already present for this polarization scenario')
        sys.exit(1)

    weight = calc_weights(mc_frame['costh_' + args.frame], args.lambda_theta)

    logging.debug('Shape of data frame before adding weights: {}'.
                  format(mc_frame.shape))
    mc_frame[weight_name] = weight
    logging.debug('Shape of data frame after adding weights: {}'
                  .format(mc_frame.shape))

    outfile = args.outfile if args.outfile else args.inputfile

    store_dataframe(mc_frame, outfile, 'chic_mc_tuple')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for reweighting a '
                                     '(unpolarized) MC sample to a desired '
                                     'other polarization')
    parser.add_argument('inputfile', help='Input file containing the '
                        'unpolarized MC sample')
    parser.add_argument('-lth', '--lambda-theta', type=float, required=True,
                        help='lambda_theta in the HX frame of the desired '
                        'polarization scenario')
    parser.add_argument('-f', '--frame', type=str, default='HX',
                        help='reference frame in which the desired pol scenario'
                        'is defined')
    parser.add_argument('-o', '--outfile', default='', type=str,
                        help='Name of output file, if it should not be the '
                        'input file')



    clargs = parser.parse_args()

    main(clargs)
