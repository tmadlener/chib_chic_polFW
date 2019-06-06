#!/usr/bin/env python
"""
Script to calculate ppd values for a given input ratio
"""

import os

from itertools import combinations

import numpy as np
import pandas as pd

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from scipy.stats import norm

from utils.graph_utils import scale_graph, get_errors
from utils.pol_utils import costh_ratio_1d, phi_ratio_1d
from utils.data_handling import store_dataframe
from utils.hist_utils import hist1d, hist2d


# Define a bin width for the 1d histogram that allows to determine the quantiles
# from the histograms deviating from the ones obtained directly from the
# unbinned sample.
# For a bin width of 0.0005 the deviations is in the order of 2.5e-4. This
# scales pretty linearly with the bin-width
BIN_WIDTH_1D = 0.0005

# The number of bins used for the 2d ppd histograms in each direction.
NBINS_2D = 200

# For now use an "unknown" but constant shift to avoid early unblinding
RAND_DLTH_SHIFT = float(os.environ['RANDOM_DELTA_LAMBDA_SHIFT'])

# The name of the ratio graph in the fit files
RATIO_NAME = 'r_chic2_chic1_v_{}_HX_fold_bin_0'

# The width and center of the Gaussian importance sampling kernel for the norm
NORM_SIGMA = 0.025
NORM_MU = 0.5

XRANGES = {
    'dlth': [76.9, 80.1], # Make sure that this covers the whole range # Make sure that this covers the whole range
    'lth': [-0.33333, 1],
    'norm': [NORM_MU - 5 * NORM_SIGMA, NORM_MU + 5 * NORM_SIGMA],
    'kappa': [-1, 1],
    'dkappa': [-2, 2],
}

def get_ratio_graph(rfile, var='costh'):
    """
    Get the ratio graph from the file and re scale it such that the first point
    lies at 0.5. This allows for an easier importance sampling of the
    normalization, since the absolute value of the normalization doesn't really
    matter.
    """
    graph = rfile.Get(RATIO_NAME.format(var))
    return scale_graph(graph, 0.5 / graph.GetY()[0])


def get_scan_points(n_points, var='costh'):
    """
    Get a tuple of np.arrays of scan values that should be used
    """
    norm_vals = np.random.normal(NORM_MU, NORM_SIGMA, n_points)
    if var == 'costh':
        lth_vals = np.random.uniform(-1./3., 1, n_points)
        # Generate the delta lambdas such that they actually are physically
        # allowed delta lambda values, by generating lambda2 and subtracting
        # lambda1 from it The alternative of generating delta lambda flat,
        # leads to coverage where it is impossible physically.
        dlth_vals = np.random.uniform(-0.6, 1, n_points) - lth_vals

        return pd.DataFrame({
            'norm': norm_vals,
            'lth': lth_vals,
            'dlth': dlth_vals + RAND_DLTH_SHIFT
        })
    if var == 'phi':
        kappa = np.random.uniform(-1, 1, n_points)
        dkappa = np.random.uniform(-1, 1, n_points) - kappa
        return pd.DataFrame({
            'norm': norm_vals,
            'kappa': kappa,
            'dkappa': dkappa
        })


def calc_chi2(graph, func):
    """
    Calculate the chi2 between the graph and all passed values of norm, lth and
    dlth (resp. kappa1 and kappa2)
    """
    x_v = np.array(graph.GetX())
    ratio_v = np.array(graph.GetY())
    _, _, err_lo, err_hi = get_errors(graph)

    func_vals = func(x_v)

    # Again let numpy broadcasting to the lifting of getting everything into the
    # correct shapes. Since ratio_v and the error arrays are only 1 dimensional
    # this works perfectly
    diff = ratio_v[:, None] - func_vals
    err = (diff < 0) * err_hi[:, None] + (diff > 0) * err_lo[:, None]

    return np.sum( (diff**2) / err**2, axis=0)


def calc_ppd(graph, dfr, var='costh'):
    """
    Calculate the ppd using the given graph and all the passed values
    """
    def fit_func_costh(costh):
        """
        Get the actual fitting function (also correct for the random dlth shift)
        with only costh as free parameter left
        """
        # let numpy broadcasting do its magic to get the costh values in
        norm_v= dfr.loc[:, 'norm'].values
        lth_v = dfr.loc[:, 'lth'].values
        dlth_v = dfr.loc[:, 'dlth'].values - RAND_DLTH_SHIFT
        return norm_v * costh_ratio_1d(costh[:, None], lth_v + dlth_v, lth_v)

    def fit_func_phi(phi):
        """
        Get the actual fitting function (also correct for the random shift)
        with only phi as parameter left
        """
        norm_v = dfr.loc[:, 'norm'].values,
        kappa_v = dfr.loc[:, 'kappa'].values
        dkappa_v = dfr.loc[:, 'dkappa'].values

        return norm_v * phi_ratio_1d(phi[:, None], kappa_v + dkappa_v, kappa_v)


    if var == 'costh':
        chi2_vals = calc_chi2(graph, fit_func_costh)
    if var == 'phi':
        chi2_vals = calc_chi2(graph, fit_func_phi)


    return np.exp(-0.5 * chi2_vals)


def get_number_bins(vrange, bwidth):
    """
    Get a number of bins that ensures that each bin is at least bwidth wide but
    also makes it possible to easily rebin the histogram in case it is needed
    """
    n_bins = np.diff(vrange)[0] / bwidth
    # Chosing a number of bins that allows to rebin the histogram with a maximum
    # factor of 200. This allows to always get histograms with 200, 100, 80, 50,
    # etc bins
    n_bins /= (25 * 16) # 5**2 * 2**4

    return np.ceil(n_bins).astype(int) * 25 * 16


def ppd_1d(data, var):
    """
    Get the 1d ppd histogram for a given variable
    """
    xran = XRANGES.get(var, None)
    bounds = {'min': xran[0], 'max': xran[1],
              'nbins': get_number_bins(xran, BIN_WIDTH_1D)}

    return hist1d(data.loc[:, var], weights=data.ppd / data.norm_weight,
                  **bounds)


def ppd_2d(data, var1, var2):
    """
    Get the 2d ppd histogram for a given variable
    """
    xran = XRANGES.get(var1)
    yran = XRANGES.get(var2)

    bounds = {'minx': xran[0], 'maxx': xran[1], 'nbinsx': NBINS_2D,
              'miny': yran[0], 'maxy': yran[1], 'nbinsy': NBINS_2D}

    return hist2d(data.loc[:, var1], data.loc[:, var2],
                  weights=data.ppd / data.norm_weight, **bounds)


def produce_ppd_hists(data, direction='costh'):
    """
    Make the PPD histograms from the passed data
    """
    variables = {
        'costh': ['lth', 'dlth', 'norm'],
        'phi': ['kappa', 'dkappa', 'norm']
    }
    hists = []
    for var in variables[direction]:
        vhist = ppd_1d(data, var)
        vhist.SetName('ppd_1d_{}'.format(var))
        hists.append(vhist)

    for var1, var2 in combinations(variables[direction], 2):
        vhist = ppd_2d(data, var1, var2)
        vhist.SetName('ppd_2d_{}_{}'.format(var1, var2))
        hists.append(vhist)

    return hists


def main(args):
    """Main"""
    ratiofile = r.TFile.Open(args.ratiofile)
    ratio_graph = get_ratio_graph(ratiofile, args.direction)

    dfr = get_scan_points(args.number_extractions, args.direction)
    # Maybe we do not need to store this since we can always calculate it from
    # the norm values if the center and width of the sampling kernel does not
    # change. But for now store it for convenience
    dfr['norm_weight'] = norm.pdf(dfr.loc[:, 'norm'], NORM_MU, NORM_SIGMA)
    dfr['ppd'] = calc_ppd(ratio_graph, dfr, args.direction)

    if not args.hists_only:
        store_dataframe(dfr, args.scanfile, tname='ppd_vals')

    outfile = r.TFile.Open(args.scanfile, 'update')
    outfile.cd()
    hists = produce_ppd_hists(dfr, args.direction)
    for hist in hists:
        hist.Write()

    outfile.Close()


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
    parser.add_argument('--hists-only', default=False, action='store_true',
                        help='Do not put the scan values into a TTree, but only '
                        'store the histograms')

    dir_sel = parser.add_mutually_exclusive_group()
    dir_sel.add_argument('--costh', action='store_const', dest='direction',
                         const='costh',
                         help='Assume that the ratio is vs costh')
    dir_sel.add_argument('--phi', action='store_const', dest='direction',
                         const='phi', help='Assume that the ratio is vs phi')
    parser.set_defaults(direction='costh')




    clargs = parser.parse_args()
    main(clargs)
