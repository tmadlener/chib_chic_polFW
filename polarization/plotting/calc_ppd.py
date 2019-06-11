#!/usr/bin/env python
"""
Script to calculate ppd values for a given input ratio
"""

import os

from itertools import combinations

import numpy as np
import pandas as pd

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')


import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from scipy.stats import norm

from utils.graph_utils import scale_graph, get_errors
from utils.pol_utils import costh_ratio_1d, phi_ratio_1d, lambda_tilde
from utils.data_handling import store_dataframe
from utils.hist_utils import hist1d, hist2d
from utils.misc_helpers import _get_var

from common_func import cond_chi1_lth_lph, cond_chi2_lth_lph


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

LTH_1_BOUNDS, LPH_1_BOUNDS = (-1./3., 1), (-1./3., 1./3.)
LTH_2_BOUNDS, LPH_2_BOUNDS = (-0.6, 1), (-np.sqrt(5) / 5, np.sqrt(5) / 5)

XRANGES = {
    'dlth': [76.8, 80.2], # Make sure that this covers the whole range # Make sure that this covers the whole range
    'lth': LTH_1_BOUNDS,
    'norm_costh': [NORM_MU - 5 * NORM_SIGMA, NORM_MU + 5 * NORM_SIGMA],
    'norm_phi': [NORM_MU - 5 * NORM_SIGMA, NORM_MU + 5 * NORM_SIGMA],
    'lph': LPH_1_BOUNDS,
    'dlph': [-1./3 - np.sqrt(5) / 5, 1./3 + np.sqrt(5) / 5],
    'lth2': [78, 79.8],
    'lph2': LPH_2_BOUNDS,
    'ltilde' : [-1, 1],
    # min / max lambda_tilde 2 = -1 / +3
    'ltilde2': [-1.05 + RAND_DLTH_SHIFT, 3.05 + RAND_DLTH_SHIFT],
    'dltilde': [-2.1 + RAND_DLTH_SHIFT, 4.1 + RAND_DLTH_SHIFT]
}


# Variables (resp. how to compute them for which 1d plots will be produced)
PPD_VARIABLES = {
    'lth': 'lth_1', 'lph': 'lph_1', #'ltp': 'ltp_1',
    'lth2': 'lth_2', 'lph2': 'lph_2', #'ltp2': 'ltp_2', # debugging

    'norm_phi': 'norm_phi',
    'norm_costh': 'norm_cth',

    'dlth': lambda d: d.lth_2 - d.lth_1,
    'dlph': lambda d: d.lph_2 - d.lph_1,

    # We have to play some serious things here to have the random shfit in place
    'ltilde': lambda d: lambda_tilde(d.lth_1, d.lph_1),
    'ltilde2': lambda d: lambda_tilde(d.lth_2 - RAND_DLTH_SHIFT, d.lph_2) + \
                         RAND_DLTH_SHIFT,
    'dltilde': lambda d: lambda_tilde(d.lth_2 - RAND_DLTH_SHIFT, d.lph_2) - \
                         lambda_tilde(d.lth_1, d.lph_1) + RAND_DLTH_SHIFT
}


# 2D combinations of available variables from above
PPD_2D_COMBS = (
    ('lth', 'lph'),
    ('lth2', 'lph2'),
    ('lth', 'dlth'),
    ('lph', 'dlph'),
)


def w_norm(dfr):
    """Get the normalization importance sampling correction weight"""
    return 1 / (dfr.norm_w_cth * dfr.norm_w_phi)


def w_ppd(dfr):
    """Get the full ppd weight (2d)"""
    return dfr.w_ppd_costh * dfr.w_ppd_phi * w_norm(dfr)


def sel_const_2d(dfr):
    """
    Get the weight constraining the 2d prior to the physically allowed regions.
    Can also be used as selection in apply_selections
    """
    return (cond_chi1_lth_lph(dfr.lth_1, dfr.lph_1) &
            cond_chi2_lth_lph(dfr.lth_2 - RAND_DLTH_SHIFT, dfr.lph_2))


# Weight (functions) for producing different prior and ppd distributions
# (Used to only have to do the majority of calculations once and then filling
# histograms accordingly to have different definitions of priors)
PPD_WEIGHTS = {
    'prior_{}': None,
    'prior_{}_norm': w_norm,
    'prior_{}_norm_costh': lambda d: 1 / d.norm_w_cth,
    'prior_{}_norm_phi': lambda d: 1 / d.norm_w_phi,
    'prior_{}_2d_const': sel_const_2d,

    'ppd_{}': w_ppd,
    'ppd_{}_2d_const': lambda d: w_ppd(d) * sel_const_2d(d),

    'ppd_{}_costh': lambda d: d.w_ppd_costh / d.norm_w_cth,
    'ppd_{}_costh_2d_const': lambda d: d.w_ppd_costh / d.norm_w_cth * sel_const_2d(d),

    'ppd_{}_phi': lambda d: d.w_ppd_phi / d.norm_w_phi,
    'ppd_{}_phi_2d_const': lambda d: d.w_ppd_phi / d.norm_w_phi * sel_const_2d(d)
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


def generate_lambdas(n_points, lth_bounds, lph_bounds):
    """
    Generate lambdas flat in 2d in box defined by lth_bounds and lph_bounds
    """
    logging.debug('Generating 2d lambda values for bounds: lth = [{:.2f}, '
                  '{:.2f}], lph = [{:.2f}, {:.2f}]'
                  .format(*(lth_bounds + lph_bounds)))

    lth = np.random.uniform(lth_bounds[0], lth_bounds[1], n_points)
    lph = np.random.uniform(lph_bounds[0], lph_bounds[1], n_points)

    return lth, lph


def get_scan_points(n_points):
    """
    Get a DataFrame of scan values that should be used for scanning the PPD
    """
    logging.info('Generating random scan points: n_points = {}'.format(n_points))
    logging.debug('Generating norm values')
    norm_phi_vals = np.random.normal(NORM_MU, NORM_SIGMA, n_points)
    norm_cth_vals = np.random.normal(NORM_MU, NORM_SIGMA, n_points)

    # lth_1, lph_1, ltp_1 = generate_lambdas(n_points, cond_chi1_lambdas, 0.16)
    # lth_2, lph_2, ltp_2 = generate_lambdas(n_points, cond_chi2_lambdas, 0.52)

    lth_1, lph_1 = generate_lambdas(n_points, LTH_1_BOUNDS, LPH_1_BOUNDS)
    lth_2, lph_2 = generate_lambdas(n_points, LTH_2_BOUNDS, LPH_2_BOUNDS)


    return pd.DataFrame({
        'norm_phi': norm_phi_vals,
        'norm_cth': norm_cth_vals,
        'lth_1': lth_1, 'lph_1': lph_1,# 'ltp_1': ltp_1,
        'lth_2': lth_2, 'lph_2': lph_2,# 'ltp_2': ltp_2,
    })


def calc_chi2(graph, func):
    """
    Calculate the chi2 between the graph and all passed values of norm, lth and
    dlth (resp. kappa1 and kappa2)
    """
    logging.debug('Calculating chi2 for graph \'{}\' and function \'{}\''
                  .format(graph.GetName(), func.__name__))
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


def calc_kappa(lth, lph, costh_range=0.625):
    """
    Calculate the kappa values in a given costh range from lth and lph
    """
    return (3 - costh_range**2) / (3 + lth * costh_range**2) * lph


def calc_ppd(cth_graph, phi_graph, dfr):
    """
    Calculate the ppd using the given graphs and all the passed values
    """
    logging.info('Calculating ppd values')
    def fit_func_costh(costh):
        """
        Get the actual fitting function (also correct for the random dlth shift)
        with only costh as free parameter left
        """
        # let numpy broadcasting do its magic to get the costh values in
        norm_v= dfr.loc[:, 'norm_cth'].values
        lth1_v = dfr.loc[:, 'lth_1'].values
        lth2_v = dfr.loc[:, 'lth_2'].values
        return norm_v * costh_ratio_1d(costh[:, None], lth2_v, lth1_v)

    def fit_func_phi(phi):
        """
        Get the actual fitting function (also correct for the random shift)
        with only phi as parameter left
        """
        norm_v = dfr.loc[:, 'norm_phi'].values
        kappa1_v = calc_kappa(dfr.lth_1, dfr.lph_1, 0.625).values
        kappa2_v = calc_kappa(dfr.lth_2, dfr.lph_2, 0.625).values

        return norm_v * phi_ratio_1d(phi[:, None], kappa2_v, kappa1_v)

    chi2_costh = calc_chi2(cth_graph, fit_func_costh)
    chi2_phi = calc_chi2(phi_graph, fit_func_phi)

    dfr['w_ppd_costh'] = np.exp(-0.5 * chi2_costh)
    dfr['w_ppd_phi'] = np.exp(-0.5 * chi2_phi)
    dfr['w_ppd'] = dfr.w_ppd_costh * dfr.w_ppd_phi

    return dfr


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


def ppd_1d(data, var, vfunc, name_fmt, weights=None):
    """
    Get the 1d ppd histogram for a given variable
    """
    xran = XRANGES.get(var, None)
    bounds = {'min': xran[0], 'max': xran[1],
              'nbins': get_number_bins(xran, BIN_WIDTH_1D),
              'name': name_fmt.format('1d_' + var)}

    return hist1d(_get_var(data, vfunc), weights=weights,
                  **bounds)


def ppd_2d(data, var1, vfunc1, var2, vfunc2, name_fmt, weights=None):
    """
    Get the 2d ppd histogram for a given variable
    """
    xran = XRANGES.get(var1)
    yran = XRANGES.get(var2)

    bounds = {'minx': xran[0], 'maxx': xran[1], 'nbinsx': NBINS_2D,
              'miny': yran[0], 'maxy': yran[1], 'nbinsy': NBINS_2D,
              'name': name_fmt.format('_'.join(['2d', var1, var2]))}

    return hist2d(_get_var(data, vfunc1), _get_var(data, vfunc2),
                  weights=weights, **bounds)


def produce_ppd_hists(data):
    """
    Make the PPD histograms from the passed data
    """
    logging.info('Producing prior and posterior probability density histograms')

    # # Calculate the weights only once, since they will not change afterwards
    weights_names = {
        n: f(data) if f is not None else None for n, f in PPD_WEIGHTS.iteritems()
    }

    hists = []


    for var, vfunc in PPD_VARIABLES.iteritems():
        logging.debug('Filling 1d histograms for {}'.format(var))

        for name, weight in weights_names.iteritems():
            hists.append(ppd_1d(data, var, vfunc, name, weight))

    for var1, var2 in PPD_2D_COMBS:
        logging.debug('Filling 2d histograms for {} v {}'.format(var1, var2))
        vfunc1, vfunc2 = PPD_VARIABLES[var1], PPD_VARIABLES[var2]

        for name, weight in weights_names.iteritems():
            hists.append(ppd_2d(data, var1, vfunc1, var2, vfunc2, name, weight))

    return hists


def main(args):
    """Main"""
    cth_file = r.TFile.Open(args.cthratiofile)
    costh_ratio = get_ratio_graph(cth_file, 'costh')
    phi_file = r.TFile.Open(args.phiratiofile)
    phi_ratio = get_ratio_graph(phi_file, 'phi')

    dfr = get_scan_points(args.number_extractions)

    # Maybe we do not need to store this since we can always calculate it from
    # the norm values if the center and width of the sampling kernel does not
    # change. But for now store it for convenience
    dfr['norm_w_cth'] = norm.pdf(dfr.loc[:, 'norm_cth'], NORM_MU, NORM_SIGMA)
    dfr['norm_w_phi'] = norm.pdf(dfr.loc[:, 'norm_phi'], NORM_MU, NORM_SIGMA)
    calc_ppd(costh_ratio, phi_ratio, dfr)

    # For now shift some values by a constant random number
    dfr['lth_2'] += RAND_DLTH_SHIFT

    if not args.hists_only:
        store_dataframe(dfr, args.scanfile, tname='ppd_vals')


    outfile = r.TFile.Open(args.scanfile, 'update')
    outfile.cd()
    hists = produce_ppd_hists(dfr)
    for hist in hists:
        hist.Write()

    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to do a scan of the '
                                     'parameter space and calculate the PPD for '
                                     'all the parameter values.')
    parser.add_argument('cthratiofile', help='File containing the input costh '
                        'ratio')
    parser.add_argument('phiratiofile', help='File containing the input phi '
                        'ratio')
    parser.add_argument('scanfile', help='Output file to which the scan values '
                        'are written')
    parser.add_argument('-n', '--number-extractions', help='Number of parameters'
                        'values that should be tested', type=int,
                        default=1000000)
    parser.add_argument('--hists-only', default=False, action='store_true',
                        help='Do not put the scan values into a TTree, but only '
                        'store the histograms')

    clargs = parser.parse_args()
    main(clargs)
