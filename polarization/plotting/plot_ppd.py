#!/usr/bin/env python
"""
Script to make plots from the ppds
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
r.TGaxis.SetMaxDigits(3)

from utils.plot_helpers import mkplot, get_y_max, setup_latex, put_on_latex
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.hist_utils import (
    get_binning, get_array, from_array, get_quantiles, rebin
)
from utils.plot_decoration import YLABELS


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
    # can.pltables[0].GetYaxis().SetMaxDigits(3)


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
    ppd = get_scaled_ppd(hfile, 'lph', 50)
    ppdmax = get_y_max(ppd)
    can = plot_lph(ppd)

    mkplot([r.TLine(v, 0, v, ppdmax * 1.1) for v in [-1./3, 1./3]], can=can,
           drawOpt='same', attr=[{'color': 12, 'line': 7, 'width': 2}])

    make_nice(can)
    return can


def shift_by_median(ppd):
    """
    Shift the ppd by the median to center it around 0
    """
    med = get_quantiles(ppd, 0.5)
    binning = get_binning(ppd)
    return from_array(get_array(ppd), binning - med, errors=get_array(ppd, errors=True))


def make_dlth_plot(hfile):
    """
    Make the dlth plot
    """
    ppd = get_scaled_ppd(hfile, 'dlth')
    ppd = rebin(shift_by_median(ppd), [(0, 200)])
    ppdmax = get_y_max(ppd)
    can = mkplot(ppd,
                 xLabel= '{0}#minus#bar{{{0}}}'.format(YLABELS['dlth']),
                 xRange=[-2, 2],
                 drawOpt='hist', yLabel='PPD [a.u.]')

    mkplot([r.TLine(v, 0, v, ppdmax * 1.1) for v in [-1.6, 1.3333]], can=can,
           drawOpt='same', attr=[{'color': 12, 'line': 7, 'width': 2}])

    ltx = setup_latex()
    put_on_latex(ltx, [
        (0.65, 0.8, '+{:.2f}'.format(np.diff(get_quantiles(ppd, [0.5, 0.84]))[0])),
        (0.40, 0.8, '#minus{:.2f}'.format(np.diff(get_quantiles(ppd, [0.16, 0.5]))[0]))
    ], ndc=True)

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
    'dlth': make_dlth_plot,
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

    plot_vars = []
    if 'costh' in args.variables:
        plot_vars.extend(['lth', 'dlth', 'norm_costh'])
    if 'phi' in args.variables:
        plot_vars.extend(['lph', 'dlph', 'norm_phi'])
    if '_' in args.variables:
        plot_vars.extend(['ltilde', 'dltilde'])

    for var in plot_vars:
        func = PLOT_FUNCTIONS[var]
        can = func(ppd_file)
        can.SaveAs('ppd_{}.pdf'.format(var))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the ppds '
                                     'calculated in calc_ppd.py')
    parser.add_argument('ppdfile', help='Root file containing the TTree from the'
                        ' scanning of the parameter space')

    var_sel = parser.add_mutually_exclusive_group()
    var_sel.add_argument('--costh', action='store_const', dest='variables',
                         const='costh', help='Make only the plots in costh resp.'
                         'make the plots for a scan done only in costh')
    var_sel.add_argument('--phi', action='store_const', dest='variables',
                         const='phi', help='Make only the plots in phi, resp. '
                         'make the plots for a scan done only in phi')
    var_sel.add_argument('--twodim', action='store_const', dest='variables',
                         const='costh_phi', help='Make all possible plots that '
                         'can be drawn from a 2d scan')
    parser.set_defaults(variables='costh')

    clargs = parser.parse_args()
    main(clargs)
