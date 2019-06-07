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
from utils.pol_utils import costh_ratio_1d, phi_ratio_1d
from utils.data_handling import store_dataframe
from utils.hist_utils import hist1d, hist2d
from utils.misc_helpers import _get_var


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
RAND_DLPH_SHIFT = float(os.environ['RANDOM_DELTA_LAMPHI_SHIFT'])

# The name of the ratio graph in the fit files
RATIO_NAME = 'r_chic2_chic1_v_{}_HX_fold_bin_0'

# The width and center of the Gaussian importance sampling kernel for the norm
NORM_SIGMA = 0.025
NORM_MU = 0.5

XRANGES = {
    'dlth': [76.8, 80.2], # Make sure that this covers the whole range # Make sure that this covers the whole range
    'lth': [-0.33333, 1],
    'norm_costh': [NORM_MU - 5 * NORM_SIGMA, NORM_MU + 5 * NORM_SIGMA],
    'norm_phi': [NORM_MU - 5 * NORM_SIGMA, NORM_MU + 5 * NORM_SIGMA],
    'lph': [-1./3., 1./3.],
    'dlph': [27.5, 29.3],
    'ltp': [-0.42, 0.42],
    'dltp': [-0.85, 0.85],
    'lth2': [78, 79.8],
    'lph2': [27.7, 28.9],
    'ltp2': [-0.45, 0.45],
}

PPD_2D_COMBS = (
    ('lth', 'lph'),
    ('lth', 'ltp'),
    ('lph', 'ltp'),

    ('lth2', 'lph2'),
    ('lth2', 'ltp2'),
    ('lph2', 'ltp2'),

    ('lth', 'dlth'),
    ('lph', 'dlph'),
    ('ltp', 'dltp'),
)

def get_ratio_graph(rfile, var='costh'):
    """
    Get the ratio graph from the file and re scale it such that the first point
    lies at 0.5. This allows for an easier importance sampling of the
    normalization, since the absolute value of the normalization doesn't really
    matter.
    """
    graph = rfile.Get(RATIO_NAME.format(var))
    return scale_graph(graph, 0.5 / graph.GetY()[0])


def generate_lambdas(n_points, condition_f, eff=None):
    """
    Generate lambdas such that they satisfy the condition posed by condition_f
    """
    logging.debug('Generating lambda values. condition_f = {}'
                  .format(condition_f.__name__))

    def gen_lambdas(n_ev):
        """Generate n_ev values in a cuboid"""
        # Generate the initial values in a cuboid that spans enough to allow
        # for all possible values for both chic1 and chic2. In this case these
        # are all defined by the chic2.
        lth = np.random.uniform(-0.6, 1, n_ev)
        lph = np.random.uniform(-np.sqrt(5) / 5, np.sqrt(5) / 5, n_ev)
        ltp = np.random.uniform(-np.sqrt(5) / 5, np.sqrt(5) / 5, n_ev)

        return lth, lph, ltp

    n_valid = 0

    n_gen_points = n_points
    # If an efficiency is provided generate more events according to
    # the efficiency
    if eff is not None and eff < 1:
        n_gen_points = int(n_points / eff)
        logging.debug('Using eff = {:.2f} and generating {} events per loop'
                      .format(eff, n_gen_points))

    valid_lth, valid_lph, valid_ltp = np.array([]), np.array([]), np.array([])
    # Keep looping until we have enough events
    while n_valid < n_points:
        lth, lph, ltp = gen_lambdas(n_gen_points)
        valid = condition_f(lth, lph, ltp)
        logging.debug('{} / {} generated tuples satisfy the condition function'
                      .format(np.sum(valid), n_gen_points))

        valid_lth = np.concatenate((valid_lth, lth[valid]))
        valid_lph = np.concatenate((valid_lph, lph[valid]))
        valid_ltp = np.concatenate((valid_ltp, ltp[valid]))

        n_valid += np.sum(valid)

    return valid_lth[:n_points], valid_lph[:n_points], valid_ltp[:n_points]


def cond_chi1_lambdas(lth, lph, ltp):
    """
    Check whether the values of lth, lph and ltp satisfy the condition of
    Eq. 27 of PRD 83, 096001 (2011)
    """
    cond1 = (lth >= -1./3.) & (lth <= 1)
    cond2 = np.abs(lph) <= ((1 - lth) * 0.25)
    cond3 = (2.25 * (lth - 1./3.)**2 + 6 * ltp**2) <= 1
    cond4 = np.abs(ltp) <= np.sqrt(3) * 0.5 * (lph + 1./3.)
    cond5 = ((lph > 1./9.) & (( (6 * lph - 1)**2 + 6 * ltp**2 ) <= 1)) | (lph <= 1./9.)

    return cond1 & cond2 & cond3 & cond4 & cond5


def cond_chi2_lambdas(lth, lph, ltp):
    """
    Check whether the values of lth, lph and ltp satisfy the condition of
    Eq. 28 of PRD 83, 096001 (2011)
    """
    return (0.3125 * (lth - 0.2)**2 + lph**2 + ltp**2) <= 0.2


def get_scan_points(n_points):
    """
    Get a DataFrame of scan values that should be used for scanning the PPD
    """
    logging.info('Generating random scan points: n_points = {}'.format(n_points))
    logging.debug('Generating norm values')
    norm_phi_vals = np.random.normal(NORM_MU, NORM_SIGMA, n_points)
    norm_cth_vals = np.random.normal(NORM_MU, NORM_SIGMA, n_points)

    lth_1, lph_1, ltp_1 = generate_lambdas(n_points, cond_chi1_lambdas, 0.16)
    lth_2, lph_2, ltp_2 = generate_lambdas(n_points, cond_chi2_lambdas, 0.52)

    return pd.DataFrame({
        'norm_phi': norm_phi_vals,
        'norm_cth': norm_cth_vals,
        'lth_1': lth_1, 'lph_1': lph_1, 'ltp_1': ltp_1,
        'lth_2': lth_2, 'lph_2': lph_2, 'ltp_2': ltp_2,
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


def calc_kappa(lth, lph, costh_range=0.625):
    """
    Calculate the kappa values in a given costh range from lth and lph
    """
    return (3 - costh_range**2) / (3 + lth * costh_range**2) * lph


def calc_ppd(cth_graph, phi_graph, dfr):
    """
    Calculate the ppd using the given graphs and all the passed values
    """
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

    return np.exp(-0.5 * (chi2_costh + chi2_phi))


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


def ppd_1d(data, var, vfunc, weights=None):
    """
    Get the 1d ppd histogram for a given variable
    """
    xran = XRANGES.get(var, None)
    bounds = {'min': xran[0], 'max': xran[1],
              'nbins': get_number_bins(xran, BIN_WIDTH_1D)}

    return hist1d(_get_var(data, vfunc), weights=weights,
                  **bounds)


def ppd_2d(data, var1, vfunc1, var2, vfunc2, weights=None):
    """
    Get the 2d ppd histogram for a given variable
    """
    xran = XRANGES.get(var1)
    yran = XRANGES.get(var2)

    bounds = {'minx': xran[0], 'maxx': xran[1], 'nbinsx': NBINS_2D,
              'miny': yran[0], 'maxy': yran[1], 'nbinsy': NBINS_2D}

    return hist2d(_get_var(data, vfunc1), _get_var(data, vfunc2),
                  weights=weights, **bounds)


def produce_ppd_hists(data):
    """
    Make the PPD histograms from the passed data
    """
    logging.info('Producing prior and posterior probability density histograms')
    variables = {
        'lth': 'lth_1', 'lph': 'lph_1', 'ltp': 'ltp_1',
        'lth2': 'lth_2', 'lph2': 'lph_2', 'ltp2': 'ltp_2', # debugging

        'norm_phi': 'norm_phi',
        'norm_costh': 'norm_cth',

        'dlth': lambda d: d.lth_2 - d.lth_1,
        'dlph': lambda d: d.lph_2 - d.lph_1,
        'dltp': lambda d: d.ltp_2 - d.ltp_1
    }

    hists = []

    # Calculate the weights only once
    w_norm = 1 / data.norm_w_cth / data.norm_w_phi
    w_ppd_norm = data.w_ppd * w_norm

    for var, vfunc in variables.iteritems():
        logging.debug('Filling 1d histogram for {}'.format(var))
        vprior = ppd_1d(data, var, vfunc, w_norm)
        vprior.SetName('prior_1d_{}'.format(var))
        hists.append(vprior)

        vppd = ppd_1d(data, var, vfunc, w_ppd_norm)
        vppd.SetName('ppd_1d_{}'.format(var))
        hists.append(vppd)

    for var1, var2 in PPD_2D_COMBS:
        logging.debug('Filling 2d histogram for {} v {}'.format(var1, var2))
        vfunc1, vfunc2 = variables[var1], variables[var2]
        vprior = ppd_2d(data, var1, vfunc1, var2, vfunc2, w_norm)
        vprior.SetName('prior_2d_{}_{}'.format(var1, var2))
        hists.append(vprior)

        vppd = ppd_2d(data, var1, vfunc1, var2, vfunc2, w_ppd_norm)
        vppd.SetName('ppd_2d_{}_{}'.format(var1, var2))
        hists.append(vppd)

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
    dfr['w_ppd'] = calc_ppd(costh_ratio, phi_ratio, dfr)

    # For now shift some values by a constant random number
    dfr['lth_2'] += RAND_DLTH_SHIFT
    dfr['lph_2'] += RAND_DLPH_SHIFT

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
