#!/usr/bin/env python
"""
Script to calculate the chi2 probability of data to Toy MC predictions
"""
import json
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from itertools import product

from utils.symbolic import lth_1, lth_2
from utils.roofit_utils import get_var_err
from utils.hist_utils import get_binning, divide
from utils.graph_utils import get_errors

# data event numbers (from costh integrated fit)
N_CHIC1 = 63139.0
N_CHIC2 = 28114.0

def get_fit_graph(wsp, costh_bins, costh_means):
    """
    Get the graph of the ratio from the workspace
    """
    ratio_central = []
    ratio_err = []

    n_bins = len(costh_bins)

    for i in xrange(n_bins):
        wsp.loadSnapshot('snap_costh_bin_{}'.format(i))
        central, err = get_var_err(wsp, 'r_chic2_chic1',
                                   wsp.genobj('fit_res_costh_bin_{}'.format(i)))
        ratio_central.append(central)
        ratio_err.append(err)

    ratio = np.array(ratio_central)
    err = np.array(ratio_err)

    c_lo = np.array([costh_means[i] - costh_bins[i][0] for i in xrange(n_bins)])
    c_hi = np.array([costh_bins[i][1] - costh_means[i] for i in xrange(n_bins)])

    return r.TGraphAsymmErrors(n_bins, np.array(costh_means), ratio,
                               c_lo, c_hi, err, err)


def calc_chi2(pred_hist, data_graph):
    """
    Calculate the chi2 between the predictor histogram and the data graph
    """
    # Need to select the bins without over- and underflow bin here
    pred = np.array([b for b in pred_hist][1:pred_hist.GetNbinsX() + 1])
    data_central = np.array(data_graph.GetY())
    _, _, data_err, _ = get_errors(data_graph)

    return np.sum((pred - data_central) / data_err)**2


def get_ratio_combs(chi1_hists, chi2_hists):
    """
    Get all the possible ratio combinations
    """
    return product(chi1_hists, chi2_hists)


def main(args):
    """Main"""

    with open(args.bininfo, 'r') as finfo:
        costh_info = json.load(finfo)

    fitfile = r.TFile.Open(args.fitres)
    wsp = fitfile.Get('ws_mass_fit')

    fit_graph = get_fit_graph(wsp, costh_info['costh_bins'],
                              costh_info['costh_means'])

    histfile = r.TFile.Open(args.histfile)
    histlist = list(set(b.GetName() for b in histfile.GetListOfKeys()))

    ratio_combs = get_ratio_combs([h for h in histlist if 'chic1' in h],
                                  [h for h in histlist if 'chic2' in h])

    for chi1_hist, chi2_hist in ratio_combs:
        # Scale individual histograms to data numbers
        h_chi1 = histfile.Get(chi1_hist)
        h_chi1.Scale(N_CHIC1 / h_chi1.Integral())
        h_chi2 = histfile.Get(chi2_hist)
        h_chi2.Scale(N_CHIC2 / h_chi2.Integral())
        ratio = divide(h_chi2, h_chi1)

        chi2 = calc_chi2(ratio, fit_graph)
        print chi1_hist, chi2_hist, chi2, r.TMath.Prob(chi2, 4)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to calculate the chi2 '
                                     'probability of data fit results to Toy MC'
                                     ' predictions')
    parser.add_argument('fitres', help='File containing the results of costh '
                        'binned mass fits')
    parser.add_argument('histfile', help='File containing the selected and '
                        'efficiency weighted costh histograms from Toy MC')
    parser.add_argument('bininfo', help='File containing the costh binning info '
                        '(as produced e.g. by the costh binnned mass fits)')

    clargs = parser.parse_args()
    main(clargs)
