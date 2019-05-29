#!/usr/bin/env python
"""
Script to make plots from the ppds
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, setup_latex, put_on_latex
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.plot_decoration import YLABELS

from calc_ppd import XRANGES

YLABELS.update({'norm': 'N'})

def make_1d_plot(hfile, var):
    """
    Make a 1d plot from the data and return the canvas
    """
    ppd = hfile.Get('ppd_1d_{}'.format(var))

    can = mkplot(ppd, xLabel=YLABELS.get(var, var), xRange=XRANGES.get(var, None),
                 yLabel='PPD value [a.u.]')
    add_auxiliary_info(can, 2012, prelim=True)

    return can


def main(args):
    """Main"""
    set_TDR_style()
    ppd_file = r.TFile.Open(args.ppdfile)

    for var in ['dlth', 'norm']:
        can = make_1d_plot(ppd_file, var)
        can.SaveAs('ppd_{}.pdf'.format(var))

    for var in ['lth']:
        can = make_1d_plot(ppd_file, var)
        can.SaveAs('ppd_{}.pdf'.format(var))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the ppds '
                                     'calculated in calc_ppd.py')
    parser.add_argument('ppdfile', help='Root file containing the TTree from the'
                        ' scanning of the parameter space')

    clargs = parser.parse_args()
    main(clargs)
