#!/usr/bin/env python
"""
Script to calculate ppd values for a given input ratio
"""

import os

import numpy as np
import pandas as pd

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from scipy.stats import norm

from utils.graph_utils import scale_graph, get_errors
from utils.pol_utils import costh_ratio_1d
from utils.data_handling import store_dataframe


# For now use an "unknown" but constant shift to avoid early unblinding
RAND_DLTH_SHIFT = float(os.environ['RANDOM_DELTA_LAMBDA_SHIFT'])

# The name of the ratio graph in the fit files
RATIO_NAME = 'r_chic2_chic1_v_costh_HX_fold_bin_0'

# The width and center of the Gaussian importance sampling kernel for the norm
NORM_SIGMA = 0.1
NORM_MU = 0.5


def get_ratio_graph(rfile):
    """
    Get the ratio graph from the file and re scale it such that the first point
    lies at 0.5. This allows for an easier importance sampling of the
    normalization, since the absolute value of the normalization doesn't really
    matter.
    """
    graph = rfile.Get(RATIO_NAME)
    return scale_graph(graph, 0.5 / graph.GetY()[0])


def get_scan_points(n_points):
    """
    Get a tuple of np.arrays of scan values that should be used
    """
    norm_vals = np.random.normal(NORM_MU, NORM_SIGMA, n_points)
    lth_vals = np.random.uniform(-0.33334, 1, n_points)
    # Generate the delta lambdas such that they actually are physically allowed
    # delta lambda values, by generating lambda2 and subtracting lambda1 from it
    # The alternative of generating delta lambda flat, leads to coverage where
    # it is impossible physically.
    dlth_vals = np.random.uniform(-0.6, 1, n_points) - lth_vals

    return norm_vals, lth_vals, dlth_vals


def calc_chi2(graph, func):
    """
    Calculate the chi2 between the graph and all passed values of norm, lth and
    dlth
    """
    costh_v = np.array(graph.GetX())
    ratio_v = np.array(graph.GetY())
    _, _, err_lo, err_hi = get_errors(graph)

    func_vals = func(costh_v)

    # Again let numpy broadcasting to the lifting of getting everything into the
    # correct shapes. Since ratio_v and the error arrays are only 1 dimensional
    # this works perfectly
    diff = ratio_v[:, None] - func_vals
    err = (diff < 0) * err_hi[:, None] + (diff > 0) * err_lo[:, None]

    return np.sum( (diff**2) / err**2, axis=0)


def calc_ppd(graph, norm_v, lth_v, dlth_v):
    """
    Calculate the ppd using the given graph and all the passed values
    """
    def fit_func(costh):
        """
        Get the actual fitting function (also correct for the random dlth shift)
        with only costh as free parameter left
        """
        # let numpy broadcasting do its magic to get the costh values in
        return norm_v * costh_ratio_1d(costh[:, None], lth_v + dlth_v, lth_v)

    chi2_vals = calc_chi2(graph, fit_func)

    return np.exp(-0.5 * chi2_vals)


def main(args):
    """Main"""
    ratiofile = r.TFile.Open(args.ratiofile)
    ratio_graph = get_ratio_graph(ratiofile)

    norm_vals, lth_vals, dlth_vals = get_scan_points(args.number_extractions)

    dfr = pd.DataFrame({
        'norm': norm_vals,
        'lth': lth_vals,
        # Before storing shift the delta lambda values by a constant random
        # shift
        'dlth': dlth_vals + RAND_DLTH_SHIFT,
        # Maybe we do not need to store this since we can always calculate it
        # from the norm values if the center and width of the sampling kernel
        # does not change
        'norm_weight': norm.pdf(norm_vals, NORM_MU, NORM_SIGMA),
        'ppd': calc_ppd(ratio_graph, norm_vals, lth_vals, dlth_vals)
    })



    store_dataframe(dfr, args.scanfile, tname='ppd_vals')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to do a scan of the '
                                     'parameter space and calculate the PPD for '
                                     'all the parameter values.')
    parser.add_argument('ratiofile', help='File containing the input ratio')
    parser.add_argument('scanfile', help='Output file to which the scan values '
                        'are written')
    parser.add_argument('-n', '--number-extractions', help='Number of parameters'
                        'values that should be tested', type=int,
                        default=1000000)


    clargs = parser.parse_args()
    main(clargs)
