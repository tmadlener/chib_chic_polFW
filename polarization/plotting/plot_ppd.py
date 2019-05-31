#!/usr/bin/env python
"""
Script to make plots from the ppds
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, default_attributes
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.plot_decoration import YLABELS
from utils.hist_utils import rebin, _get_y_max_hist


YLABELS.update({'norm': 'N'})
ATTR = default_attributes(size=0)

def get_scaled_ppd(hfile, var, nbins=100):
    """
    Get the ppd scaled to unity
    """
    ppd = hfile.Get('ppd_1d_{}'.format(var))
    ppd = rebin(ppd, [(0, nbins)])
    ppd.Scale(1 / ppd.Integral())
    return ppd


def make_lth_plot(hfile):
    """
    Make the lambda1 plot
    """
    ppd = get_scaled_ppd(hfile, 'lth', 50)
    ppdmax = _get_y_max_hist(ppd)

    can = mkplot(ppd, xLabel=YLABELS['lth'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[-1, 1])
    add_auxiliary_info(can, 2012, prelim=True)
    mkplot(r.TLine(-0.3333, 0, -0.3333, ppdmax * 1.1, ), can=can, drawOpt='same',
           attr=[{'color': 12, 'line': 7, 'width': 2}])

    phist =  can.pltables[0]
    phist.GetYaxis().SetMaxDigits(3)

    return can


def make_dlth_plot(hfile):
    """
    Make the delta_lambda plot
    """
    ppd = get_scaled_ppd(hfile, 'dlth')
    ppdmax = ppd.GetBinCenter(ppd.GetMaximumBin())

    can = mkplot(ppd, xLabel=YLABELS['dlth'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[ppdmax - 2, ppdmax + 2])
    add_auxiliary_info(can, 2012, prelim=True)

    return can


def make_norm_plot(hfile):
    """
    Make the norm plot
    """
    ppd = get_scaled_ppd(hfile, 'norm')

    can = mkplot(ppd, xLabel=YLABELS['norm'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[0.375, 0.625])
    add_auxiliary_info(can, 2012, prelim=True)

    return can


PLOT_FUNCTIONS = {
    'lth': make_lth_plot,
    'dlth': make_dlth_plot,
    'norm': make_norm_plot
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
