#!/usr/bin/env python
"""
Script for checking the additional effect of reconstruction efficiencies on
top of the acceptance cuts using MC (gen + reco)
"""

# import numpy as np

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from gen_mc_costh_singlemu_acc import (
    get_costh_hist, get_phi_hist, get_selections
)

from utils.data_handling import get_dataframe
# from utils.plot_helpers import mkplot, default_colors
from utils.PlotServer import PlotServer


test_selections = [
    'no_sel',
    'jpsi_kin_sel', 'jpsi_kin_sel_loose', 'jpsi_kin_sel_fiducial',
    'jpsi_kin_sel_flat'
]



def create_histos(df, selection, gen, var, frame):
    """
    Get the chi1, chi2 and ratio histograms for a given variable and selection.
    If gen is True, the gen level variables are used
    """
    if var == 'phi':
        hist_func = get_phi_hist
    if var == 'costh':
        hist_func = get_costh_hist

    chi1_data = df[df.wChic1 == 1]
    chi2_data = df[df.wChic2 == 1]

    chi1_hist = hist_func(chi1_data, frame, selection, gen)
    chi2_hist = hist_func(chi2_data, frame, selection, gen)
    ratio_hist = chi2_hist.Clone()
    ratio_hist.SetYTitle("#chi_{c2} / #chi_{c1}")
    ratio_hist.Divide(chi1_hist)

    return (chi1_hist, chi2_hist, ratio_hist)


def store_histos(pserver, histos, top_dir, sel_str, var, frame, gen):
    """Store on set of histograms"""
    tdir = '_'.join([top_dir, sel_str])
    if gen:
        tdir += '_gen'

    pserver.store_hist(histos[0], tdir, 2012, '8_Jpsi', 0, var, frame, 'chic1')
    pserver.store_hist(histos[1], tdir, 2012, '8_Jpsi', 0, var, frame, 'chic2')
    pserver.store_hist(histos[2], tdir, 2012, '8_Jpsi', 0, var, frame, 'ratio')


def main(args):
    """Main"""
    data = get_dataframe(args.infile, args.treename)
    pserver = PlotServer(args.outfile, 'update')

    selections = get_selections((args.ptmin, args.ptmax))

    for frame in args.frames.split(','):
        for var in ['phi', 'costh']:
            for sel in test_selections:
                selection = selections[sel]
                histset = create_histos(data, selection, args.genlevel, var,
                                        frame)
                store_histos(pserver, histset, args.prefix, sel, var, frame,
                             args.genlevel)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for checking the '
                                     'additional effects of reconstruction '
                                     'efficiency on top of single muon '
                                     'acceptances (cuts) using gen and reco MC')
    parser.add_argument('infile', help='root file containing the mc data (reco '
                        'or gen).')
    parser.add_argument('outfile', help='root file into which the histograms '
                        'will be stored')
    parser.add_argument('--ptmin', type=float, default=8, help='minimum jpsi pt')
    parser.add_argument('--ptmax', type=float, default=20, help='maximum jpsi pt')
    parser.add_argument('-t', '--treename', default='chic_mc_tuple', type=str,
                        help='name of the input tree')
    parser.add_argument('-f', '--frames', help='Comma separated list of frames',
                        default='HX,PX,CS')
    parser.add_argument('-p', '--prefix', help='prefix to the top level directory')
    parser.add_argument('-g', '--genlevel', default=False, action='store_true',
                        help='Use gen level variables for reco MC')


    clargs = parser.parse_args()
    main(clargs)
