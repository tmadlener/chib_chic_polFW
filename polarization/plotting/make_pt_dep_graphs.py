#!/usr/bin/env python
"""
Script to make the pT dependent graphs from the individual ppds.
"""

from collections import OrderedDict

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.data_base import JsonDataBase
from utils.hist_utils import get_quantiles

from common_func import get_scaled_ppd


DATABASE = JsonDataBase()

SIGMA_BANDS = {
    1: [0.16, 0.84],
    2: [0.025, 0.975],
    3: [0.005, 0.995]
}

CL_BANDS = {
    90: 0.1,
    95: 0.05,
    99: 0.01
}

def open_files(inputfiles):
    """
    Open the input files and return them ordered by pt
    """
    pt_files = []
    for pt_file in inputfiles:
        pt_part, filen = pt_file.split(':')
        min_pt, max_pt = [int(v) for v in pt_part.split(',')]
        pt_files.append(((min_pt, max_pt), filen))

    files = OrderedDict()
    # Assume that sorting the bins by the lower boundary is enough
    for pt_bin, filen in sorted(pt_files, key=lambda x: x[0][0]):
        files[pt_bin] = r.TFile.Open(filen)

    return files


def get_ppds(infiles, variable):
    """
    Get the ppds for the passed variable
    """
    ppds = OrderedDict()
    for pt, rfile in infiles.iteritems():
        ppds[pt] = get_scaled_ppd(rfile, variable)
    return ppds


def get_uncertainty_graph(ppds, n_sigma, center_at_0=False):
    """
    Get the graph with the uncertainties corresponding to n_sigma
    """
    pt_vals = np.array([DATABASE.get_mean_pt(2012, p) for p in ppds.keys()])
    pt_lo = pt_vals - np.array([p[0] for p in ppds.keys()])
    pt_hi = np.array([p[1] for p in ppds.keys()]) - pt_vals

    med = np.array([get_quantiles(p, 0.5) for p in ppds.values()])
    hi_lo = np.array(
        [get_quantiles(p, SIGMA_BANDS[n_sigma]) for p in ppds.values()]
    )
    val_lo = med - hi_lo[:, 0]
    val_hi = hi_lo[:, 1] - med

    if center_at_0:
        med = np.zeros_like(med)

    return r.TGraphAsymmErrors(len(ppds), pt_vals, med, pt_lo, pt_hi,
                               val_lo, val_hi)


def get_exclusion_graph(ppds, limit):
    """
    Get the graph (TGraph) corresponding to the lower limit of the passed PPDs
    """
    pt_vals = np.array([DATABASE.get_mean_pt(2012, p) for p in ppds.keys()])
    pt_lo = pt_vals - np.array([p[0] for p in ppds.keys()])
    pt_hi = np.array([p[1] for p in ppds.keys()]) - pt_vals

    lim = np.array([get_quantiles(p, CL_BANDS[limit]) for p in ppds.values()])


    return r.TGraphAsymmErrors(len(ppds), pt_vals, lim, pt_lo, pt_hi)


def main(args):
    """Main"""
    infiles = open_files(args.infiles)
    ppds = get_ppds(infiles, args.variable)

    graphfile = r.TFile(args.graphfile, 'recreate')
    graphfile.cd()

    if args.variable in ['dlth', 'dlph']:
        for n_sig in [1, 2, 3]:
            graph = get_uncertainty_graph(ppds, n_sig, args.variable == 'dlth')
            graph.SetName('{}_v_pt_n_sig_{}'.format(args.variable, n_sig))
            graph.Write()
    else:
        for clb in [90, 95, 99]:
            graph = get_exclusion_graph(ppds, clb)
            graph.SetName('{}_v_pt_cl_{}'.format(args.variable, clb))
            graph.Write()

    graphfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make pT dependent '
                                     'graphs from the individual ppds')
    parser.add_argument('graphfile', help='Output file containing 1, 2 and 3 '
                        'sigma graphs')
    parser.add_argument('-i', '--input', action='append', required=True,
                        help='Add an input file. Format is: min_pt,max_pt:file',
                        dest='infiles')
    parser.add_argument('-v', '--variable', default='dlth', help='Variable name '
                        'for which the plot should be produced (dlth, dlph or '
                        'lth)')



    clargs = parser.parse_args()
    main(clargs)
