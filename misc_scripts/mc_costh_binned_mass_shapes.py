#!/usr/bin/env python
"""
Script to produce some plots for the study showing the mass shape vs costh on MC
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.data_handling import get_dataframe
from utils.misc_helpers import get_costh_binning, get_bin_cut_df
from utils.plot_helpers import mkplot, default_attributes, setup_legend
from utils.hist_utils import create_histogram

def create_mass_hist(dfr, costh_bin=None):
    """
    Create a mass histogram in a given costh_bin (or integrated if None)
    """
    if costh_bin is not None:
        plot_data = dfr[get_bin_cut_df(dfr, lambda d: d.costh_HX.abs(),
                                       *costh_bin)]
    else:
        plot_data = dfr

    return create_histogram(plot_data.chicMass, (75, 3.325, 3.725))


def main(args):
    """Main"""
    dfr = get_dataframe(args.mcfile)
    costh_bins = get_costh_binning(dfr, args.nbins)

    h_int = create_mass_hist(dfr)
    h_int.Scale(1.0 / args.nbins)
    h_bins = [create_mass_hist(dfr, cbin) for cbin in costh_bins]

    attr = [{'color': 1, 'marker': 20, 'size': 1.0}]
    attr += default_attributes(size=1.0, open_markers=True) * 3

    legEntries = ['costh integrated']
    legEntries += ['{:.2f} < |cos#vartheta^{{HX}}| < {:.2f}'.format(*cbin) for cbin in costh_bins]

    leg = setup_legend(0.675, 0.6, 0.95, 0.94)
    leg.SetTextSize(0.025)

    can = r.TCanvas('plot_can', '', 600, 600)
    can.SetRightMargin(0.02)
    can.SetTopMargin(0.05)
    can.SetLeftMargin(0.12)

    can = mkplot([h_int] + h_bins, can=can, xRange=[3.325, 3.725],
                        yRange=[0.1, None], attr=attr, logy=args.logy,
                 xLabel='M^{#chi}', yLabel='# events', legOpt='PLE',
                 leg=leg, legEntries=legEntries)

    can.SaveAs(args.output)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to produce some plots '
                                     'for studying the mass shape vs costh')
    parser.add_argument('mcfile', help='Input MC file')
    parser.add_argument('-n', '--nbins', default=6, help='Number of costh '
                        'bins', type=int)
    parser.add_argument('-o', '--output', default='costh_binned_MC_mass_dists.pdf',
                        help='Name of the created output plot')
    parser.add_argument('--logy', action='store_true', default=False,
                        help='Use log y axis')

    clargs = parser.parse_args()
    main(clargs)
