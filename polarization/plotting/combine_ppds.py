#!/usr/bin/env python
"""
Script to combine PPDs into a final PPD
"""

from os.path import dirname, basename
from itertools import combinations

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.hist_utils import get_binning, get_array, from_array, rebin
from utils.setup_plot_style import set_basic_style
from utils.plot_helpers import set_color, mkplot, default_attributes

from common_func import (
    get_scaled_ppd, get_scaled_ppd_2d, plot_lth, plot_dlth, plot_norm
)


def get_combined_ppd(inputfiles, var):
    """
    Get the combined ppd from all inputfiles
    """
    ppds = [get_scaled_ppd(f, var) for f in inputfiles]
    # PPDs all have the same binning
    ppd_binning = get_binning(ppds[0])

    ppd_vals = np.array([get_array(p) for p in ppds])
    ppd_errs = np.array([get_array(p, errors=True) for p in ppds])

    # Get the maximum value in each gin and its uncertainty
    max_idx = np.argmax(ppd_vals, axis=0)
    # Necessary for 2d indexing. There might be an easier way for this
    idx = np.arange(0, len(ppd_vals[0]))
    max_ppd = ppd_vals[max_idx, idx]
    max_err = ppd_errs[max_idx, idx]

    return from_array(max_ppd, ppd_binning, errors=max_err)


def get_combined_ppd_2d(inputfiles, var1, var2):
    """
    Get the combined 2d ppd from all inputfiles
    """
    ppds = [get_scaled_ppd_2d(f, var1, var2, 100, 100) for f in inputfiles]
    ppd_binning = np.array([get_binning(ppds[0], 0), get_binning(ppds[0], 1)])

    ppd_vals = np.array([get_array(p) for p in ppds])
    # TODO: at some point find out how argmax works in multiple dimensions

    return from_array(np.max(ppd_vals, axis=0), ppd_binning)


PLOT_FUNCTIONS = {
    'lth': plot_lth, 'dlth': plot_dlth, 'norm': plot_norm
}

def make_debug_plot(comb_ppd, inputfiles, var):
    """
    Make a plot overlaying the combined ppd and the input ppds
    """
    can = PLOT_FUNCTIONS[var](rebin(comb_ppd, [(0, 200)]))
    set_color(can.pltables[1], 1) # make the combined ppd black

    ppds = [get_scaled_ppd(f, var, 200) for f in inputfiles]
    mkplot(ppds, can=can, drawOpt='samehist',
           attr=default_attributes(linewidth=1))

    return can


def main(args):
    """Main"""
    set_basic_style()
    infiles = [r.TFile.Open(f) for f in args.inputfiles]

    outfile = r.TFile(args.outfile, 'recreate')
    outfile.cd()

    for var in ['lth', 'dlth', 'norm']:
        comb_ppd = get_combined_ppd(infiles, var)
        comb_ppd.SetName('ppd_1d_{}'.format(var))
        comb_ppd.Write()

        if args.debug_plots:
            outdir = dirname(args.outfile)
            filebase = basename(args.outfile).replace('.root', '')
            can = make_debug_plot(comb_ppd, infiles, var)
            can.SaveAs('{}/{}_comb_ppd_{}.pdf'.format(outdir, filebase, var))

    for var1, var2 in combinations(['lth', 'dlth', 'norm'], 2):
        comb_ppd = get_combined_ppd_2d(infiles, var1, var2)
        comb_ppd.SetName('ppd_2d_{}_{}'.format(var1, var2))
        comb_ppd.Write()

    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to combine PPDs')
    parser.add_argument('outfile', help='File into which the combined PPD is '
                        'stored')
    parser.add_argument('inputfiles', nargs='+', help='PPD scan files containing'
                        ' histograms that are used for combining')
    parser.add_argument('--debug-plots', action='store_true', default=False,
                        help='Create some debug overview plots of the 1d ppds')


    clargs = parser.parse_args()
    main(clargs)
