#!/usr/bin/env python
"""
Script to make plots of the PPDs combining all three pT ranges into one plot
"""

import os
from collections import OrderedDict
from numpy import sqrt

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
r.TGaxis.SetMaxDigits(3)

from utils.hist_utils import rebin, get_quantiles
from utils.plot_helpers import (
    mkplot, get_y_max, default_attributes, setup_legend
)
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info

RDIR = os.path.join(os.environ['MY_DATA_DIR'], 'ChicPol', 'Chic2012', 'Results')
# path of the ppd files relative to the nominal_fit dir
PPD_FILES = {
    'costh': {
        (8, 12): 'fit_nominal_model/ppd.root',
        (12, 18): 'ppd.root',
        (18, 30): 'fit_nominal_model_5bins/ppd.root'
    }, 'phi': {
        (8, 12): 'ppd.root',
        (12, 18): 'ppd.root',
        (18, 30): 'fit_nominal_model/ppd.root'
    }
}

OUTDIR = os.path.join(RDIR, 'thesis_plots')

# attributes for the vertical lines indicating physical boundaries
BOUND_ATTR = [{'color': 12, 'line': 7, 'width': 2}]
MAX_DLPH = 1./3. + 1. / sqrt(5)

def open_ppd_files(var):
    """Open the ppd files for the given variable"""
    files = OrderedDict()
    for pt in [(8, 12), (12, 18), (18, 30)]:
        files[pt] = r.TFile.Open(
            os.path.join(RDIR, 'fits_jpsipt_{}_{}'.format(*pt), var, 'nominal_fit',
                         PPD_FILES[var][pt])
        )

    return files


def add_info(can, pos='right'):
    add_auxiliary_info(can, 2012, prelim=True, pos=pos)
    return can


def get_scaled_ppd(hfile, var, nbins=None):
    """
    Get the ppd scaled to unity
    """
    ppd = hfile.Get('ppd_1d_{}'.format(var))
    if nbins is not None:
        ppd = rebin(ppd, [(0, nbins)])
    ppd.Scale(1 / ppd.Integral())
    return ppd


def make_lth_plot(ppd_files):
    """Make the plots for the lth1 ppds"""
    ppds = [get_scaled_ppd(f, 'lth') for f in ppd_files.values()]
    quantiles = [get_quantiles(p, 0.1) for p in ppds]
    ppds = [rebin(p, [(0, 100)]) for p in ppds] # for plotting

    leg = setup_legend(0.18, 0.7, 0.48, 0.84)
    leg.SetTextSize(0.0425)
    leg.SetFillColor(0)
    leg.SetFillStyle(1001)

    # Draw boundaries first to not draw it over the legend in the end
    can = mkplot(r.TLine(-1./3., 0, -1./3., get_y_max(ppds) * 1.1),
                 attr=BOUND_ATTR,
                 xRange=[-1, 1], xLabel='#lambda_{#vartheta}^{#chi_{c1}}',
                 yRange=[0, get_y_max(ppds) * 1.1], yLabel='PPD [a.u.]')

    can = mkplot(ppds, drawOpt='hist same', can=can,
                 leg=leg, legOpt='L',
                 legEntries=['{} < p_{{T}} < {} GeV'.format(*p) for p in ppd_files.keys()])
    # can.pltables[0].GetYaxis().SetNdivisions(505)
    can.pltables[0].GetXaxis().SetNdivisions(505)

    mkplot([r.TLine(q, 0, q, 0.01) for q in quantiles], can=can, drawOpt='same')

    arrow_attr = default_attributes(open_markers=False)
    for att in arrow_attr:
        att['fill'] = att['color']
        att['width'] = 2

    mkplot([
        r.TArrow(q, v, q + 0.125, v, 0.025, '<|') for q, v in zip(quantiles, [0.0015, 0.003, 0.0045])
    ],
           can=can, drawOpt='same <|', attr=arrow_attr)

    return add_info(can, 'left')


def make_lph_plot(ppd_files):
    """Make the plots for lph1 ppds"""
    ppds = [get_scaled_ppd(f, 'lph', 50) for f in ppd_files.values()]

    can = mkplot([r.TLine(v, 0, v, get_y_max(ppds) * 1.1) for v in [-1/3., 1/3.]],
                 xRange=[-1, 1], xLabel='#lambda_{#varphi}^{#chi_{c1}}',
                 yRange=[0, get_y_max(ppds) * 1.1], yLabel='PPD [a.u.]',
                 attr=BOUND_ATTR)

    leg = setup_legend(0.39, 0.25, 0.75, 0.39)
    leg.SetTextSize(0.0425)
    leg.SetFillStyle(1001)

    can = mkplot(ppds, drawOpt='hist same', can=can,
                 leg=leg, legOpt='L',
                 legEntries=['{} < p_{{T}} < {} GeV'.format(*p) for p in ppd_files.keys()])
    # can.pltables[0].GetYaxis().SetNdivisions(505)
    can.pltables[0].GetXaxis().SetNdivisions(505)

    return add_info(can)


def make_dlth_plot(ppd_files):
    """Make the plots for Delta lambda theta"""
    ppds = [get_scaled_ppd(f, 'dlth', 200) for f in ppd_files.values()]

    can = mkplot([r.TLine(v, 0, v, get_y_max(ppds) * 1.1) for v in [-1.6, 1.333]],
                 attr=BOUND_ATTR,
                 xRange=[-2, 2], xLabel='#Delta#lambda_{#vartheta}',
                 yRange=[0, get_y_max(ppds) * 1.1], yLabel='PPD [a.u.]')


    leg = setup_legend(0.56, 0.7, 0.88, 0.84)
    leg.SetTextSize(0.0425)
    leg.SetFillStyle(1001)

    can = mkplot(ppds, drawOpt='hist same', can=can,
                 leg=leg, legOpt='L',
                 legEntries=['{} < p_{{T}} < {} GeV'.format(*p) for p in ppd_files.keys()])

    return add_info(can)


def make_dlph_plot(ppd_files):
    """Make the plots for Delta lambda phi"""
    ppds = [get_scaled_ppd(f, 'dlph', 400) for f in ppd_files.values()]

    can = mkplot([r.TLine(v, 0, v, get_y_max(ppds) * 1.1) for v in [-MAX_DLPH, MAX_DLPH]],
                 attr=BOUND_ATTR,
                 xRange=[-0.85, 0.85], xLabel='#Delta#lambda_{#varphi}',
                 yRange=[0, get_y_max(ppds) * 1.1], yLabel='PPD [a.u.]')

    leg = setup_legend(0.18, 0.7, 0.48, 0.84)
    leg.SetTextSize(0.0425)
    leg.SetFillStyle(1001)

    can = mkplot(ppds, drawOpt='hist same', can=can,
                 leg=leg, legOpt='L',
                 legEntries=['{} < p_{{T}} < {} GeV'.format(*p) for p in ppd_files.keys()])
    can.pltables[0].GetXaxis().SetNdivisions(507)

    return add_info(can, 'left')


def main():
    """Main"""
    set_TDR_style()
    r.gStyle.SetPadTopMargin(r.gStyle.GetPadTopMargin() + 0.005)

    ppd_files = open_ppd_files('costh')
    can = make_lth_plot(ppd_files)
    can.SaveAs(os.path.join(OUTDIR, 'ppd_lth_w_lower_lim_90_combined.pdf'))

    can = make_dlth_plot(ppd_files)
    can.SaveAs(os.path.join(OUTDIR, 'ppd_dlth_combined.pdf'))


    ppd_files = open_ppd_files('phi')
    can = make_lph_plot(ppd_files)
    can.SaveAs(os.path.join(OUTDIR, 'ppd_lph_combined.pdf'))

    can = make_dlph_plot(ppd_files)
    can.SaveAs(os.path.join(OUTDIR, 'ppd_dlph_combined.pdf'))


if __name__ == '__main__':
    main()
