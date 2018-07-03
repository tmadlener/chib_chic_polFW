#!/usr/bin/env python
"""
Script that generates the photon efficiency curves and stores them in a root
file.

For the moment only the pT curves for the different eta bins are created
"""

import re
import json
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np
import sympy as sp

from utils.symbolic import func_cov
from utils.graph_utils import get_lower_band, get_upper_band

from common_func import get_name

# Covariance matrix from the fit integrated over the whole eta range, where
# alpha and beta were fixed. This will be used to calculate the correlation
# coefficients between the fitted parameters, which will then be used to get
# the uncertainty bands for the parametrization
COVARIANCE = np.array([
    [1.181e-06,  1.545e-06,  -4.328e-06, 4.156e-06],
    [1.545e-06,  7.215e-06,  -1.714e-05, 5.177e-06],
    [-4.328e-06, -1.714e-05, 4.228e-05,  -1.481e-05],
    [4.156e-06,  5.177e-06,  -1.481e-05, 1.506e-05],
])

# corr = diag(cov)^{-1/2} * cov * diag(cov)^{-1/2}
CORRELATIONS = np.matmul(
    np.matmul(
        np.diag(1/np.sqrt(np.diag(COVARIANCE))), COVARIANCE,
    ), np.diag(1/np.sqrt(np.diag(COVARIANCE)))
)


def eff_param_string():
    """
    The parametrization of the efficiencies from AN-2015-11 as a string that can
    be used in a TF1 constructor.

    p0 * (1 - p1 * (Erf(pT + p2) - p1 / alpha * (pT - p3 * (pT^2 - p3 / beta * pT^3))))
    """
    return '[0] * (1 - [1] * (TMath::Erf(x[0] + [2]) - [1] / [4] * (x[0] - [3] * (pow(x[0], 2) - [3] / [5] * pow(x[0], 3)))))'


def eff_param():
    """
    Get the parametrization as ROOT.TF1
    """
    return r.TF1('photon_eff_param', eff_param_string(), 0, 7)


def eff_param_sym():
    """
    Get the parametrization as sympy symbolic expression by doing some string
    manipulation on the parametrization and then using sympy.sympify
    """
    param_str = eff_param_string()
    # replace call to ROOTs erf and give x[0] a parseable name
    param_str = param_str.replace('TMath::Erf', 'erf').replace('x[0]', 'x')
    # convert parameters from [x] notation to px notation
    param_str = re.sub(r'\[([0-9])\]', r'p\1', param_str)
    # replace pow(x, y) with x**y (pythonic) syntax
    param_str = re.sub(r'pow\((.*?)\s*?,\s*?([0-9])\)', r'\1**\2', param_str)

    return sp.sympify(param_str)


def get_corr_subs_values(corr):
    """
    Get the dictionary of substitution values for the correlation matrix
    """
    subs_dict = {}
    n_dim = corr.shape[0]
    for irow in xrange(0, n_dim):
        for icol in xrange(irow + 1, n_dim):
            subs_dict['rho_p{}p{}'.format(irow, icol)] = corr[irow, icol]

    return subs_dict


def get_cov_func(params, corr):
    """
    Get the uncertainty function where only pT is left as a free parameter.

    This will return a python function that can be evaluated at any given point
    """
    eff = eff_param_sym()
    # get the list of free parameters
    free_params = []
    for sym in eff.free_symbols:
        if sym.name in params and params[sym.name][1] != 0:
            free_params.append(sym)

    # sort the parameters according to their name, such that the correlation
    # coefficients actually match
    free_params.sort(key=lambda x: int(x.name.replace('p', '')))

    cov_eff = func_cov(eff, free_params)

    # build up the dictionary of symbol -> value that will be substituted.
    # In the end the returned function will only have one free parameter left
    subst_vals = {
        p: v[0] for p, v in params.iteritems()
    }
    subst_vals.update({
        'sigma_' + p: v[1] for p, v in params.iteritems()
    })
    subst_vals.update(
        get_corr_subs_values(corr)
    )

    # NOTE: here it is assumed that 'x' is the only free parameter left
    return sp.lambdify(sp.symbols('x'), cov_eff.subs(subst_vals))


def get_graph_err(params, corr, n_sigma=1.0, n_points=100):
    """
    Get the function evaluated at n_points with uncertainties taking into
    account correlations between the parameters
    """
    # central function
    eff_f = eff_param_sym()
    eff_f = eff_f.subs({p: v[0] for p, v in params.iteritems()})
    # NOTE: assume that 'x' is the only free parameter left
    eff_f = sp.lambdify(sp.symbols('x'), eff_f)

    # uncertainty function (as function of pT)
    var_f = get_cov_func(params, corr)

    x_bins = np.linspace(0.4, 7, n_points + 1)
    x_cent = 0.5 * (x_bins[1:] + x_bins[:-1]) # bin centers
    x_err = np.diff(x_bins) # "uncertainties" in x

    y_cent = np.array([eff_f(x) for x in x_cent])
    y_err = np.sqrt(np.array([var_f(x) for x in x_cent])) * n_sigma

    return r.TGraphErrors(len(x_cent), x_cent, y_cent, x_err, y_err)


def set_params_errors(func, *params):
    """
    Set all the parameters as pairs of value and uncertainty (in the order they)
    are in the params list. If uncertainty = 0, the parameter is fixed
    """
    central = np.array([p[0] for p in params])
    uncer = np.array([p[1] for p in params])

    func.SetParameters(central)
    func.SetParErrors(uncer)

    for idx, err in enumerate(uncer):
        if err == 0:
            func.FixParameter(idx, func.GetParameter(idx))


def load_params(param_file):
    """
    Load the parameter file and return the list of dicts stored in it
    """
    with open(param_file, 'r') as pfile:
        eff_params = json.load(pfile)
        return eff_params


def create_param(params, sigma_shift, uncorrelated):
    """
    Create the function from the passed params and give it an appropriate name
    """
    # if the central results are desired. Use the exact parametrization as TF1
    if sigma_shift == 0:
        func = eff_param()
        set_params_errors(func, params["p0"], params["p1"], params["p2"],
                          params["p3"], params["alpha"], params["beta"])

        func.SetName(get_name(params["eta"], 'photon_eff_pt'))
        return func

    # else get an aproximation by evaluating the function at a given number of
    # points and determine the uncertainties at these points, then store the
    # points as a TGraph where the y-values are the central + uncertainty values
    # at each evaluation point

    # NOTE: Since eff_param_sym names alpha and beta p4 and p5 respectively
    # (can't use beta in an expression that goes through sympy.sympify), we have
    # to clone them here. We can leave the original values in, since they will
    # not be picked up by the substitution command
    params['p4'] = params['alpha']
    params['p5'] = params['beta']

    # use the global correlation matrix or an identity matrix if uncorrelated
    # parameters are desired
    corr = np.identity(4) if uncorrelated else CORRELATIONS
    graph = get_graph_err(params, corr, np.abs(sigma_shift), 200)

    if sigma_shift < 0:
        graph = get_lower_band(graph)
    else:
        graph = get_upper_band(graph)

    graph.SetName(get_name(params['eta'], 'photon_eff_pt'))

    return graph


def main(args):
    """Main"""
    file_option = 'update' if args.update else 'recreate'
    outfile = r.TFile.Open(args.outfile, file_option)

    all_params = load_params(args.paramfile)
    for params in all_params:
        eff = create_param(params, args.sigma, args.uncorrelated)
        eff.Write()

    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to generate TF1 '
                                     'photon efficiency parametrizations from '
                                     'json file holding the fit parameters')
    parser.add_argument('paramfile', help='json file containing the fitted '
                        'parameters')
    parser.add_argument('-o', '--outfile', help='root file into which the TF1 '
                        'should be stored', default='photon_effs_param.root')
    parser.add_argument('-u', '--update', help='update the output file instead '
                        'of recreating it', default=False, action='store_true')
    parser.add_argument('-s', '--sigma', help='Use the central value + [sigma] '
                        '* uncertainty for each parameter', type=float,
                        default=0)
    parser.add_argument('--uncorrelated', default=False, action='store_true',
                        help='Assume that the free parameters are uncorrelated '
                        'instead of using correlation parameters from a global '
                        'fit')

    clargs = parser.parse_args()
    main(clargs)
