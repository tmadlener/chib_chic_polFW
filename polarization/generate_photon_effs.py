#!/usr/bin/env python
"""
Script that generates the photon efficiency curves and stores them in a root
file.

For the moment only the pT curves for the different eta bins are created
"""

import json
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

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


def set_params_errors(func, p0, p1, p2, p3, alpha, beta):
    """
    All the parameters as pairs of value and uncertainty.
    If uncertainty = 0, the parameter is fixed
    """
    def cond_set_fix(func, pidx, pval, perr):
        """Conditionally set the value and the error or fix the value"""
        func.SetParameter(pidx, pval)
        func.SetParError(pidx, perr)
        if perr == 0:
            func.FixParameter(pidx, pval)

    for idx, val_err in enumerate((p0, p1, p2, p3, alpha, beta)):
        cond_set_fix(func, idx, *val_err)


def load_params(param_file):
    """
    Load the parameter file and return the list of dicts stored in it
    """
    with open(param_file, 'r') as pfile:
        eff_params = json.load(pfile)
        return eff_params


def create_param(params):
    """
    Create the function from the passed params and give it an appropriate name
    """
    func = eff_param()

    set_params_errors(func,
                      params["p0"], params["p1"], params["p2"], params["p3"],
                      params["alpha"], params["beta"])

    func.SetName(get_name(params["eta"], 'photon_eff_pt'))
    return func


def main(args):
    """Main"""
    file_option = 'update' if args.update else 'recreate'
    outfile = r.TFile.Open(args.outfile, file_option)

    all_params = load_params(args.paramfile)
    for params in all_params:
        eff = create_param(params)
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

    clargs = parser.parse_args()
    main(clargs)
