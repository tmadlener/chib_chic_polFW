#!/usr/bin/env python
"""
Script to make plots from the ppds
"""

import numpy as np

from utils.data_handling import get_dataframe
from utils.plot_helpers import mkplot, setup_latex, put_on_latex
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.hist_utils import hist1d, hist2d
from utils.plot_decoration import YLABELS
from utils.misc_helpers import quantile

YLABELS.update({'norm': 'N'})

XRANGES = {
    'dlth': [-4.7, -1.75],
    'lth': [-0.33333, 1],
    'norm': [0.5 - 5*0.025, 0.5 + 5*0.025]
}
NBINS_1D = 100


def ppd_1d(data, var):
    """
    Get the 1d ppd histogram for a given variable
    """
    xran = XRANGES.get(var, None)
    bounds = {'min': xran[0] if xran is not None else None,
              'max': xran[1] if xran is not None else None,
              'nbins': NBINS_1D}

    return hist1d(data.loc[:, var], weights=data.ppd / data.norm_weight,
                  **bounds)


def add_median_bands(can, data, var):
    """
    Add information about the observed median and a 68 % CI onto the
    canvas
    """
    loci, med, hici = quantile(data.loc[:, var].values, [0.16, 0.5, 0.84],
                               (data.ppd / data.norm_weight).values)

    lerr, herr = med - loci, hici - med
    val_text = '{:.4f}^{{+{:.4f}}}_{{#minus{:.4f}}}'.format(med, lerr, herr)

    var_text = YLABELS.get(var, var) + ' = ' + val_text

    ltx = setup_latex()
    put_on_latex(ltx, [(0.725, 0.65, var_text)], True)


def add_low_limit(can, data, var):
    """
    Add information about the observed lower limit (90 %) onto the canvas
    """
    lolim = quantile(data.loc[:, var].values, 0.1,
                     (data.ppd / data.norm_weight).values)

    var_text = '{} > {:.4f} (90 % C.I.)'.format(YLABELS.get(var, var), lolim)

    ltx = setup_latex()
    put_on_latex(ltx, [(0.675, 0.45, var_text)], True)


def make_1d_plot(data, var, add_info_f=None):
    """
    Make a 1d plot from the data and return the canvas
    """
    ppd = ppd_1d(data, var)



    can = mkplot(ppd, xLabel=YLABELS.get(var, var), xRange=XRANGES.get(var, None),
                 yLabel='PPD value [a.u.]')
    add_auxiliary_info(can, 2012, prelim=True)

    if add_info_f is not None:
        add_info_f(can, data, var)

    return can


def main(args):
    """Main"""
    set_TDR_style()
    ppd_data = get_dataframe(args.ppdfile)

    for var in ['dlth', 'norm']:
        can = make_1d_plot(ppd_data, var, add_median_bands)
        can.SaveAs('ppd_{}.pdf'.format(var))

    for var in ['lth']:
        can = make_1d_plot(ppd_data, var, add_low_limit)
        can.SaveAs('ppd_{}.pdf'.format(var))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the ppds '
                                     'calculated in calc_ppd.py')
    parser.add_argument('ppdfile', help='Root file containing the TTree from the'
                        ' scanning of the parameter space')

    clargs = parser.parse_args()
    main(clargs)
