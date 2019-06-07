#!/usr/bin/env python
"""
Script to make plots from the ppds
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, get_y_max
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info

from common_func import (
    get_scaled_ppd, plot_lth, plot_dlth, plot_norm, plot_dlph, plot_lph,
    plot_ltilde, plot_dltilde
)


PLOT_FUNC = {
    'dlth': plot_dlth,
    'dlph': plot_dlph,
    'lph': plot_lph,
    'norm_phi': plot_norm,
    'norm_costh': plot_norm,
    'ltilde': plot_ltilde,
    'dltilde': plot_dltilde
}

def make_nice(can):
    """
    Add auxiliary info and adapt y-axis labeling
    """
    add_auxiliary_info(can, 2012, prelim=True)
    can.pltables[0].GetYaxis().SetMaxDigits(3)


def make_lth_plot(hfile):
    """
    Make the lambda1 plot
    """
    ppd = get_scaled_ppd(hfile, 'lth', 100)
    ppdmax = get_y_max(ppd)
    can = plot_lth(ppd)

    mkplot(r.TLine(-0.3333, 0, -0.3333, ppdmax * 1.1, ), can=can, drawOpt='same',
           attr=[{'color': 12, 'line': 7, 'width': 2}])

    make_nice(can)
    return can


def make_lph_plot(hfile):
    """
    Make the lambda_phi 1 plot
    """
    ppd = get_scaled_ppd(hfile, 'lph', 100)
    ppdmax = get_y_max(ppd)
    can = plot_lph(ppd)

    mkplot([r.TLine(v, 0, v, ppdmax * 1.1) for v in [-1./3, 1./3]], can=can,
           drawOpt='same', attr=[{'color': 12, 'line': 7, 'width': 2}])

    make_nice(can)
    return can


def make_simple_plot(var, n_bins=200):
    """
    Make an ordinary 1d plot, without anything fancy going on
    """
    def _make_plot(hfile):
        """
        Make the plot
        """
        ppd = get_scaled_ppd(hfile, var, n_bins)
        can = PLOT_FUNC[var](ppd)
        make_nice(can)

        return can

    return _make_plot


PLOT_FUNCTIONS = {
    'lth': make_lth_plot,
    'lph': make_lph_plot,
    'dlth': make_simple_plot('dlth'),
    'dlph': make_simple_plot('dlph', 400),
    'norm_costh': make_simple_plot('norm_costh'),
    'norm_phi': make_simple_plot('norm_phi'),
    'ltilde': make_simple_plot('ltilde'),
    'dltilde': make_simple_plot('dltilde'),
}


def main(args):
    """Main"""
    set_TDR_style()
    ppd_file = r.TFile.Open(args.ppdfile)

    for var, func in PLOT_FUNCTIONS.iteritems():
        can = func(ppd_file)
        can.SaveAs('ppd_{}.pdf'.format(var))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the ppds '
                                     'calculated in calc_ppd.py')
    parser.add_argument('ppdfile', help='Root file containing the TTree from the'
                        ' scanning of the parameter space')

    clargs = parser.parse_args()
    main(clargs)
