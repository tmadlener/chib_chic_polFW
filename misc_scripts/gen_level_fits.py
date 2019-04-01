#!/usr/bin/env python
"""
Script to do a gen-level fit
"""

import re
import pickle

from collections import OrderedDict

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.hist_utils import project
from utils.pol_utils import w_costh_phi

GENNAME = 'fold_costh_phi_JpsiPt_JpsiRap_gen_HX'

def deduce_scenario(filename):
    """
    Deduce the scenario to which the input file belongs
    """
    frame, scen = None, None
    if 'HX' in filename and not 'CS' in filename:
        frame = 'HX'
    if 'CS' in filename and not 'HX' in filename:
        frame = 'CS'

    rrg = r'([0-9]+(o[0-9]+)?)'
    scen_rgx = re.compile(r'chic(1_R_{}|2_R1_{}_R2_{})'.format(rrg, rrg, rrg))
    match = scen_rgx.search(filename)
    if match:
        scen = match.group(0)

    if frame is None:
        frame = raw_input('Could not deduce frame from filename. Please enter '
                          'frame: ')
    if scen is None:
        scen = raw_input('Could not deduce scen from filename. Please enter '
                         'polarization scenario: ')

    return frame, scen


def fit_to_dist(dist, func, verbose=False):
    """
    Fit func to dist
    """
    opts = 'S0'
    if not verbose:
        opts += 'q'
    fit_res = dist.Fit(func, opts)
    return fit_res.Chi2(), fit_res.Ndf()


def update_fit_file(filename, fitresults, frame, scen):
    """
    Update the file containing fit results from other fits or create it if it
    is not yet present
    """
    try:
        with open(filename, 'r') as fitfile:
            data = pickle.load(fitfile)
    except IOError:
        data = OrderedDict()

    if frame not in data:
        data[frame] = OrderedDict()

    if scen not in data[frame]:
        data[frame][scen] = OrderedDict()

    data[frame][scen].update(fitresults)

    with open(filename, 'w') as fitfile:
        pickle.dump(data, fitfile)


def fit_dict(func, chi2, ndf):
    """
    Create an OrderedDict containing the par values and the ndf and chi2 vals
    """
    vals = OrderedDict()
    for i in [1, 2, 3]:
        vals[func.GetParName(i)] = func.GetParameter(i)
        vals[func.GetParName(i) + '_err'] = func.GetParError(i)

    vals['chi2'] = chi2
    vals['ndf'] = ndf

    return vals


def main(args):
    """Main"""
    frame, scen = deduce_scenario(args.genlevelfile)

    genfile = r.TFile.Open(args.genlevelfile)
    genhist = project(genfile.Get(GENNAME), [1, 0])

    func = w_costh_phi(set_vals={'norm': 1e5})
    chi2, ndf = fit_to_dist(genhist, func, args.verbose)

    fit_res = fit_dict(func, chi2, ndf)
    update_fit_file(args.outfile, fit_res, frame, scen)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Scrip to do gen-level fits')
    parser.add_argument('genlevelfile', help='Generation level histogram file')
    parser.add_argument('outfile', help='output pickle file into which the fit '
                        'results are stored')
    parser.add_argument('-v', '--verbose', help='Print verbose fit',
                        action='store_true', default=False)

    clargs = parser.parse_args()
    main(clargs)
