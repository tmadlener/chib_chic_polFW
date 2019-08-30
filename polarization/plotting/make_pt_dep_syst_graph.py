#!/usr/bin/env python
"""
Script to generate and store the graph with the systematic uncertainties
as a function of pT
"""

import json
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True


def create_graph(systfile, var):
    """
    Read the json file and create the graph with systematic uncertainties.
    The y values will be the systematic uncertainties, while the x-axis will
    have the bin-centers as values and the the bin widths covering the full bin.
    """
    with open(systfile, 'r') as sysf:
        data = json.load(sysf)

        n_bins = len(data[var])
        bin_low, bin_high, errs = [], [], []

        for (low, high), err in data[var]:
            bin_low.append(low)
            bin_high.append(high)
            errs.append(err)

        bin_low, bin_high = np.array(bin_low), np.array(bin_high)
        errs = np.array(errs)
        centers = 0.5 * (bin_low + bin_high)
        x_errs = 0.5 * (bin_high - bin_low)

        return r.TGraphErrors(n_bins, centers, errs, x_errs)


def main(args):
    """Main"""
    graphfile = r.TFile(args.outfile, 'recreate')
    graph = create_graph(args.systfile, args.variable)
    graph.SetName('{}_v_pt_syst'.format(args.variable))
    graph.Write()
    graphfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to generate systematic '
                                     'graphs and store them in a file')
    parser.add_argument('systfile', help='json file containing the systematic '
                        'uncertainties')
    parser.add_argument('outfile', help='File that will be created with the '
                        'systematic uncertainty graph in it')

    var_sel = parser.add_mutually_exclusive_group()
    var_sel.add_argument('--dlth', action='store_const', dest='variable',
                         const='dlth')
    var_sel.add_argument('--dlph', action='store_const', dest='variable',
                         const='dlph')
    parser.set_defaults(variable='dlth')

    clargs = parser.parse_args()
    main(clargs)
