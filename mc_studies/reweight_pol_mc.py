#!/usr/bin/env python
"""
Script for reweighting MC to any desired polarization
"""

import sys

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.data_handling import get_dataframe, store_dataframe
from utils.pol_utils import ang_dist_lth, ang_dist_2d
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


def calc_weights_2d(costh, phi, lth, lph):
    """
    Calculate the weights to reweight the sample to have the desired
    lambda_theta and lambda_phi in the frame in which costh and phi are defined.

    Args:
        costh (float or numpy.array)
        phi (float or numpy.array)
        lth (float): lambda_theta of the desired polarization, assuming that it
            is specified in the same frame as costh and phi
        lph (float): lambda_phi of the desired polarization, assuming that it is
            specified in the same frame as costh and phi

    Returns:
        float or numpy.array: The weight that has to be applied to the event to
            get a polarized sample, with lth and lph as polarization parameters.
            Return type is dependent on the type of costh and phi
    """
    logging.debug('Adding weights for lambda_theta = {} and lambda_phi = {}'
                  .format(lth, lph))
    return ang_dist_2d(costh, phi, (lth, lph, 0))


def main(args):
    """Main"""
    lth, lph, frame = args.lambda_theta, args.lambda_phi, args.frame
    weight_name = get_storable_name('wPol_{}_lth_{:.2f}'.format(frame, lth))

    if lph is not None:
        weight_name = get_storable_name('wPol_{}_lth_{:.2f}_lph_{:.2f}'
                                        .format(frame, lth, lph))

    mc_frame = get_dataframe(args.inputfile)
    if weight_name in mc_frame.columns:
        logging.info('weights already present for this polarization scenario')
        sys.exit(1)

    if lph is None:
        weight = calc_weights(mc_frame['costh_' + frame], lth)
    else:
        weight = calc_weights_2d(mc_frame['costh_' + frame],
                                 mc_frame['phi_' + frame], lth, lph)

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
    parser.add_argument('-lph', '--lambda-phi', type=float, default=None,
                        help='lambda_phi in the specified frame of the desired '
                        'polarization scenario')
    parser.add_argument('-f', '--frame', type=str, default='HX',
                        help='reference frame in which the desired pol scenario'
                        'is defined')
    parser.add_argument('-o', '--outfile', default='', type=str,
                        help='Name of output file, if it should not be the '
                        'input file')



    clargs = parser.parse_args()

    main(clargs)
