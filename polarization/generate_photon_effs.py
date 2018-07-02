#!/usr/bin/env python
"""
Script that generates the photon efficiency curves and stores them in a root
file.

For the moment only the pT curves for the different eta bins are created
"""

import json
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np

from common_func import get_name


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


def set_params_errors(func, *params, **kwargs):
    """
    Set all the parameters as pairs of value and uncertainty (in the order they)
    are in the params list. If uncertainty = 0, the parameter is fixed
    """
    sigma = kwargs.pop('sigma', 0)
    central = np.array([p[0] for p in params])
    uncer = np.array([p[1] for p in params])

    params = central + sigma * uncer
    func.SetParameters(params)
    if sigma == 0: # NOTE: only setting the uncertainties for the central result
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


def create_param(params, sigma_shift):
    """
    Create the function from the passed params and give it an appropriate name
    """
    func = eff_param()

    set_params_errors(func,
                      params["p0"], params["p1"], params["p2"], params["p3"],
                      params["alpha"], params["beta"],
                      sigma=sigma_shift)

    func.SetName(get_name(params["eta"], 'photon_eff_pt'))
    return func


def main(args):
    """Main"""
    file_option = 'update' if args.update else 'recreate'
    outfile = r.TFile.Open(args.outfile, file_option)

    all_params = load_params(args.paramfile)
    for params in all_params:
        eff = create_param(params, args.sigma)
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

    clargs = parser.parse_args()
    main(clargs)
