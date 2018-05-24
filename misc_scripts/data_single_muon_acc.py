#!/usr/bin/env python
"""
Script for making single muon acceptance cut plots
"""

# from decorator import decorator
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gStyle.SetPadRightMargin(0.15)
r.gStyle.SetPadLeftMargin(0.12)
# r.gROOT.SetBatch() # comment for ipython


from utils.data_handling import get_dataframe, create_histogram
from utils.plot_helpers import plot_on_canvas
from utils.misc_helpers import get_bin_cut_df, cond_mkdir, create_random_str


def basic_sel(df, mc=False, gen=False):
    selection = np.ones(df.shape[0], dtype=bool)
    if not gen:
        selection &= get_bin_cut_df(df, 'chicMass', 3.325, 3.725)
    if not mc and not gen:
        selection &= df.vtxProb > 0.01

    # selection &= (df.costh_HX.abs() < 0.5)
    # selection &= (df.wChic1 == 1)
    return selection


def jpsi_kin_sel(df, jpsiPt, noPt=False, noRap=False):
    ptname = 'JpsiPt'
    rapname = 'JpsiRap'
    selection = np.ones(df.shape[0], dtype=bool)
    if not noPt:
        selection &= get_bin_cut_df(df, ptname, *jpsiPt)
    if not noRap:
        selection &= df[rapname].abs() < 1.2
    return selection

def fiducial_cuts():
    """Get the coordinates representing the standard fiducial cuts"""
    coords = (
        (20, 1.6),
        (3.0, 1.6), (3.0, 1.4), (3.5, 1.4), (3.5, 1.2), (4.5, 1.2),
        (4.5, 0)
    )

    return coords


def flat_pt(pt=3.5, eta=1.6):
    """
    Get coordinates representing flat pt cut
    """
    return ((20, eta), (pt, eta), (pt, 0))


def loose_cuts():
    """Loosest cut that does not have any acceptance holes in pt-eta"""
    coords = (
        (20, 1.6),
        (2, 1.6), (3.5, 1.2), (3.5, 0)
    )

    return coords


def draw_line(coords):
    """
    Draw a line through all coordinates using TLines
    """
    lines = []
    for start, end in zip(coords[:-1], coords[1:]):
        lines.append(r.TLine(start[0], start[1], *end))
    return lines




def make_pt_eta_plot(plotvar, cuts, cut_plot_set, savename):
    """
    Make pt-eta plot (2D) for passed variables and overlay the passed cuts
    """
    pt_eta_hist = (50, 0, 20, 24, 0, 2.4)
    pt_eta = create_histogram(plotvar, pt_eta_hist)

    pt_eta.SetXTitle('p_{T}')
    pt_eta.SetYTitle('|#eta|')

    can = r.TCanvas(create_random_str(), 'c', 600, 600)
    can.cd()
    pt_eta.Draw('colz')


    # lines = {c: [] for c in cuts} # keep the TLines alive long enough for saving
    lines = []
    for cut in cuts.values():
        lines.append(draw_line(cut()))

    for i, cut in enumerate(lines):
        plot_on_canvas(can, cut, attr=(cut_plot_set[i],))

    leg = r.TLegend(0.5, 0.3, 0.8, 0.4)
    leg.SetTextColor(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    for i, cut in enumerate(cuts.keys()):
        leg.AddEntry(lines[i][0], cut, 'l')
    leg.Draw()

    can.Update()
    # can.SaveAs(savename)

    return can, pt_eta, lines # for ipython so that they survive


def make_overview_plots(data, outdir, mc=False, gen=False):
    """
    Make a set of 2D overview plots for different selected cuts
    """
    cuts = {
        'loose': loose_cuts,
        'fiducial': fiducial_cuts
    }
    plot_sett = ({'color': 0, 'width': 2, 'line': 1},
                 {'color': 0, 'width': 2, 'line': 2},
                 {'color': 0, 'width': 2, 'line': 4},
                 {'color': 0, 'width': 2, 'line': 6})

    # applying abs to pt doesn't change anything since we have negative values only
    make_pt_eta_plot(np.abs(data[['muP_pt', 'muP_eta']]),
                     cuts, plot_sett, '{}/muP_pt_eta_no_sel.pdf'.format(outdir))
    make_pt_eta_plot(np.abs(data[['muN_pt', 'muN_eta']]),
                     cuts, plot_sett, '{}/muN_pt_eta_no_sel.pdf'.format(outdir))

    make_pt_eta_plot(np.abs(data.loc[basic_sel(data, mc, gen), ['muP_pt', 'muP_eta']]),
                     cuts, plot_sett, '{}/muP_pt_eta_basic_sel.pdf'.format(outdir))
    make_pt_eta_plot(np.abs(data.loc[basic_sel(data, mc, gen), ['muN_pt', 'muN_eta']]),
                     cuts, plot_sett, '{}/muN_pt_eta_basic_sel.pdf'.format(outdir))

    make_pt_eta_plot(np.abs(data.loc[jpsi_kin_sel(data, (8, 20)), ['muP_pt', 'muP_eta']]),
                     cuts, plot_sett, '{}/muP_pt_eta_jpsi_kin_sel.pdf'.format(outdir))
    make_pt_eta_plot(np.abs(data.loc[jpsi_kin_sel(data, (8, 20)), ['muN_pt', 'muN_eta']]),
                     cuts, plot_sett, '{}/muN_pt_eta_jpsi_kin_sel.pdf'.format(outdir))

    make_pt_eta_plot(np.abs(data.loc[basic_sel(data, mc, gen) & jpsi_kin_sel(data, (8, 20)), ['muP_pt', 'muP_eta']]),
                     cuts, plot_sett, '{}/muP_pt_eta_basic_jpsi_kin_sel.pdf'.format(outdir))
    make_pt_eta_plot(np.abs(data.loc[basic_sel(data, mc, gen) & jpsi_kin_sel(data, (8, 20)), ['muN_pt', 'muN_eta']]),
                     cuts, plot_sett, '{}/muN_pt_eta_basic_jpsi_kin_sel.pdf'.format(outdir))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for checking effects '
                                     'of single muon acceptance cuts')
    parser.add_argument('datafile', help='File containing the flat data tuple')
    # parser.add_argument('--ptmin', type=float, default=8, help='minimum jpsi pt')
    # parser.add_argument('--ptmax', type=float, default=50, help='maximum jpsi pt')
    parser.add_argument('-t', '--treename', default='chic_tuple', type=str,
                        help='name of the tree in which the original variables are '
                        '(used for storing output file).')
    parser.add_argument('-o', '--outdir', help='output directory', default='.')
    parser.add_argument('-g', '--genlevel', help='data is generator level', default=False,
                        action='store_true')
    parser.add_argument('-mc', help='data is mc', default=False, action='store_true')

    args = parser.parse_args()

    data = get_dataframe(args.datafile, args.treename)
    # right now defining a lambda for this, maybe there is a way to apply the abs function only partially
    # pt_eta = lambda df, c: np.array([df['mu{}_pt'.format(c)], df['mu{}_eta'.format(c)].abs()])

    cond_mkdir(args.outdir)
    make_overview_plots(data, args.outdir, args.mc or args.genlevel, args.genlevel)
