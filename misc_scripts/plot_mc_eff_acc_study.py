#!/usr/bin/env python
"""
Script for creating the plots of the acceptance - efficiency studies
"""

import itertools

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.PlotServer import PlotServer
from utils.plot_helpers import mkplot, default_attributes, get_y_max
from utils.misc_helpers import create_random_str, cond_mkdir
from utils.setup_plot_style import set_TDR_style

# plot gen curves as filled areas
gen_attributes = [
    {'color': 1, 'size': 0, 'fillalpha': (1, 0.5), 'width': 0}
]
comp_attributes = default_attributes(size=1, linewidth=1)

# provide consistent styles for the different


def get_plot_histos(pserver, input_sets, plot):
    """
    inputsets (list of tuples): which mc, reco, etc. and which selection
    selections (list), plot (tuple with var, frame, plot)

    Get all the plots for all selections
    """
    hists = {}
    for iset, sel in input_sets:
        input_id = '_'.join([iset, sel])
        hists[input_id] = pserver.get_hist(input_id, 2012, '8_Jpsi', 0,
                                           *plot)

    return hists


def setup_legend():
    """
    Setup the legend
    """
    leg = r.TLegend(0.2, 0.15, 0.3, 0.3)
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetEntrySeparation(0.01)
    leg.SetBorderSize(0)
    return leg


def do_ratio(ref, comp):
    """Ratio comp to ref"""
    ratio = comp.Clone('ratio')
    ratio.Divide(ref)
    ratio.GetYaxis().SetTitle('reco / gen')
    ratio.GetYaxis().SetTitleSize(0.08)
    ratio.GetYaxis().SetLabelSize(0.08)
    ratio.GetXaxis().SetTitleSize(0.08)
    ratio.GetXaxis().SetLabelSize(0.08)

    return ratio


def make_ratio_plot(can, hists, leg, legentries, **kwargs):
    """
    Wrapper around
    """
    if legentries:
        mkplot(hists, can=can, leg=leg, legEntries=legentries, **kwargs)
    else:
        mkplot(hists, can=can, **kwargs)


def make_plot(hists, savename, reference, doleg=False):
    """
    Make plot
    """
    set_TDR_style()

    can = r.TCanvas(create_random_str(16), '', 50, 50, 600, 600)
    can.cd()

    # preliminary collect them and delete entries if necessary
    legentries = []

    with_ref = False
    if reference is not None:
        ref_hist = hists[reference]
        with_ref = True
        rest_hists = [hists[k] for k in hists if k != reference]
        pull_hists = [do_ratio(ref_hist, h) for h in rest_hists]
        legentries = [reference] + [n for n in hists if n != reference]
    else:
        rest_hists = hists.values()
        legentries = hists.keys()

    if doleg:
        for i, entry in enumerate(legentries):
            legentries[i] = entry.replace('jpsi_kin_sel', '').replace('__', '_')
    else:
        legentries = []

    y_max = get_y_max(hists.values()) * 1.1
    leg=setup_legend()

    if with_ref:
        pad = r.TPad('ratio_pad', 'ratio_pad', 0, 0.3, 1, 1)
        r.SetOwnership(pad, False)
        pad.Draw()
        pad.cd()
        # passing all legentries here, since indexing doesn't work with an empty
        # list and mkplot picks only the first in this case
        make_ratio_plot(pad, [ref_hist], leg, legentries,
                        yRange=[0, y_max], attr=gen_attributes, drawOpt='E2',
                        legOpt='F')
        make_ratio_plot(pad, rest_hists, leg, legentries[1:],
                        yRange=[0, y_max], attr=comp_attributes, drawOpt='sameE1',
                        legOpt='PLE')

        can.cd()
        pull_pad = r.TPad('pull_pad', 'pull_pad', 0, 0, 1, 0.3)
        r.SetOwnership(pull_pad, False)
        pull_pad.Draw()
        pull_pad.cd()
        make_ratio_plot(pull_pad, pull_hists, None, [], yRange=[None, None],
                        attr=comp_attributes)

    else:
        make_ratio_plot(can, rest_hists, leg, legentries,
                        yRange=[0, y_max], attr=comp_attributes)

    can.SaveAs(savename)


def main(args):
    """Main"""
    pserver = PlotServer(args.inputfile)
    cond_mkdir(args.outdir)

    input_set_sel = list(itertools.product(args.inputsets.split(','),
                                           args.selections.split(',')))

    for plot in itertools.product(args.variables.split(','),
                                  args.frames.split(','),
                                  args.plots.split(',')):
        hists = get_plot_histos(pserver, input_set_sel, plot)

        plot_savename = '_'.join([args.plotbase] + list(plot)) + args.extension
        savename = '/'.join([args.outdir, plot_savename])
        make_plot(hists, savename, args.reference, args.legend)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to make plots')

    parser.add_argument('inputfile', help='file containing the histograms')
    parser.add_argument('-f', '--frames', help='frames for which to create the '
                        'plots (one per frame)', default='HX,CS')
    parser.add_argument('-p', '--plots', help='plots to produce', default='ratio')
    parser.add_argument('-v', '--variables', help='which variable to plot',
                        default='costh')
    parser.add_argument('-o', '--outdir', help='output directory for the plots',
                        default='.')
    parser.add_argument('-s', '--selections', help='which selections to plot',
                        default='jpsi_kin_sel_loose')
    parser.add_argument('-i', '--inputsets', help='which input sets to plot',
                        default='gen')
    parser.add_argument('-r', '--reference', help='inputset_selection to be used'
                        ' as a reference in the plots', default=None)
    parser.add_argument('-b', '--plotbase', help='basename of plots',
                        default='comp_gen_reco_mc')
    parser.add_argument('-l', '--legend', help='Do a legend', default=False,
                        action='store_true')
    parser.add_argument('-x', '--extension', help='file extension/plot format '
                        'to use', default='.pdf')
    
    clargs = parser.parse_args()
    main(clargs)
