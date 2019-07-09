#!/usr/bin/env python
"""
Script to make plots comparing the PPDs for different variations
"""

import json
from collections import OrderedDict

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.hist_utils import (
    get_quantiles, rebin, from_array, get_array, get_binning
)
from utils.plot_helpers import (
    mkplot, setup_latex, put_on_latex, default_colors, default_attributes,
    setup_legend
)
from utils.setup_plot_style import set_basic_style
from utils.misc_helpers import cond_mkdir
from utils.plot_decoration import YLABELS


# Quantiles defining the lower bound of the 1 sigma uncertainty band,
# the central value and the upper bound of the 1 sigma uncertainty band
QUANTILES = [0.16, 0.5, 0.84]

VAR_ATTR = [{'color': 1, 'size': 0.75, 'marker': 20}]

def get_scaled_ppd(hfile, var, nbins=None):
    """
    Get the ppd scaled to unity
    """
    ppd = hfile.Get('ppd_1d_{}'.format(var))
    if nbins is not None:
        ppd = rebin(ppd, [(0, nbins)])
    ppd.Scale(1 / ppd.Integral())
    return ppd


def open_files(files_list):
    """
    Open all the files
    """
    files = OrderedDict()
    for filen, lab in files_list:
        files[lab] = r.TFile.Open(filen)

    return files


def get_nom_graph(num_variations, err_lo, err_hi):
    """
    Get the nominal graph (centered around 0)
    """
    y_val = np.array([0.5 * num_variations - 0.5])
    return r.TGraphAsymmErrors(1, np.array([0]), y_val,
                               np.array([err_lo]), np.array([err_hi]),
                               y_val + 0.5, y_val + 0.5)


def calc_err(nom_err, var_err):
    """Calculate approximately correlation corrected uncertanities"""
    return np.sqrt(np.abs(nom_err**2 - var_err**2))


def var_graph_wrt_nominal(nom_quants, var_ppd, y_val):
    """
    Get the graph of the variation wrt nominal values
    """
    var_mean, var_rms = var_ppd.GetMean(), var_ppd.GetRMS()
    err = np.array([calc_err(var_rms, nom_quants[1])])

    return r.TGraphErrors(1, np.array([var_mean - nom_quants[0]]),
                          np.array([y_val], dtype=float), err)

    # var_quants = get_quantiles(var_ppd, QUANTILES)
    # var_lo, var_hi = np.diff(var_quants)
    # nom_lo, nom_hi = np.diff(nom_quants)

    # errlo = np.array([calc_err(var_lo, nom_lo)], dtype=float)
    # errhi = np.array([calc_err(var_hi, nom_hi)], dtype=float)

    # return r.TGraphAsymmErrors(1, np.array([var_quants[1] - nom_quants[1]]),
    #                            np.array(y_val, dtype=float),
    #                            errlo, errhi)

def get_var_graphs_wrt_nominal(var_files, nom_quants, var):
    """
    Get all the variation graphs
    """
    var_ppds = [get_scaled_ppd(f, var) for f in var_files.values()]
    return [
        var_graph_wrt_nominal(nom_quants, p, i) for i, p in enumerate(var_ppds)
    ]


def get_label(var):
    """
    Get the x-axis label
    """
    var_str = YLABELS.get(var, var)
    return '{0} - {0}_{{nominal}}'.format(var_str)


def get_labels(labels):
    """
    Produce a list of labels that can be used for latex labels
    """
    lab_list = []
    lpos = 0.05
    deltay = 0.82 / (len(labels) - 0.5)
    offy = 0.13 + deltay * 0.25

    for ilab, lab in enumerate(labels):
        lab_list.append((lpos, ilab * deltay + offy, lab))

    return lab_list


def make_var_plot(nom_file, var_files, var):
    """
    Make a variation plot comparing the nominal PPD to the variations PPDs
    """
    nom_ppd = get_scaled_ppd(nom_file, var)
    nom_quants = get_quantiles(nom_ppd, QUANTILES)
    nom_lo, nom_hi = np.diff(nom_quants)
    nom_graph = get_nom_graph(len(var_files), nom_lo, nom_hi)

    # var_graphs = get_var_graphs_wrt_nominal(var_files, nom_quants, var)
    var_graphs = get_var_graphs_wrt_nominal(var_files, [nom_ppd.GetMean(), nom_ppd.GetRMS()], var)

    xran = np.max([nom_lo, nom_hi]) * 1.75
    yran = [-0.25, len(var_files) - 0.75]

    can = mkplot(nom_graph, attr=[{'fillalpha': (default_colors()[1], 0.5)}],
                 drawOpt='E2',
                 xLabel=get_label(var), xRange=[-xran, xran], yRange=yran)
    mkplot(r.TLine(0, yran[0], 0, yran[1]), can=can, drawOpt='same',
           attr=default_attributes(widht=2)[1:])
    mkplot(var_graphs, drawOpt='samePE', can=can, attr=VAR_ATTR)

    phist = can.pltables[0]
    phist.GetYaxis().SetLabelSize(0)
    phist.GetYaxis().SetTickSize(0)

    ltx = setup_latex()
    put_on_latex(ltx, get_labels(var_files.keys()), ndc=True)
    ltx.Draw()
    can.Draw()

    return can


def shift_by_median(ppd, median):
    """Shift the ppd by the passed median"""
    return from_array(get_array(ppd),
                      get_binning(ppd) - median,
                      errors=get_array(ppd, errors=True))


def create_legend(xlo, xhi, yhi, n_plots):
    """
    Create a legend that fits all of the keys
    """
    leg = setup_legend(xlo, yhi - n_plots * 0.045, xhi, yhi)
    leg.SetTextSize(0.03)
    return leg


def make_ppd_comp_plot(nom_file, var_files, var):
    """
    Make a plot comparing the full ppds directly (shifting all of them by the
    median of the nominal ppd).
    """
    n_bins = {'dlth': 100, 'dlph': 400, 'dkappa': 400}

    nom_ppd = get_scaled_ppd(nom_file, var)
    nom_med = get_quantiles(nom_ppd, 0.5)
    nom_ppd = shift_by_median(rebin(nom_ppd, [(0, n_bins[var])]), nom_med)
    nom_ppd.Scale(1.0 / nom_ppd.Integral())

    var_ppds = [get_scaled_ppd(f, var, n_bins[var]) for f in var_files.values()]
    var_ppds = [shift_by_median(p, nom_med) for p in var_ppds]
    [p.Scale(1.0 / p.Integral()) for p in var_ppds]

    xran = {'dlth': [-2, 2], 'dlph': [-0.5, 0.5], 'dkappa': [-0.5, 0.5]}

    leg = create_legend(0.69, 0.88, 0.94, len(var_ppds) + 1)

    label = '{0} - #bar{0}_{{nominal}}'.format(YLABELS.get(var))

    can = mkplot(nom_ppd, attr=[{'color': 1, 'fillalpha': (1, 0.25)}],
                 drawOpt='hist', yLabel='PPD [a.u.]',
                 xLabel=label, xRange=xran[var], yRange=[0, None],
                 leg=leg, legEntries=['nominal'], legOpt='F')

    mkplot(var_ppds, can=can, drawOpt='histsame',
           attr=default_attributes(linewidth=2),
           leg=leg, legEntries=var_files.keys(), legOpt='L')
    return can


def main(args):
    """Main"""
    set_basic_style()
    r.gStyle.SetPadLeftMargin(0.4)


    with open(args.plotconfig, 'r') as conff:
        plot_config = json.load(conff)

    var_files = open_files(plot_config['variations'])
    nom_file = r.TFile.Open(plot_config['nominal'])

    cond_mkdir(args.outdir)
    can = make_var_plot(nom_file, var_files, args.variable)
    can.SaveAs('{}/comp_ppd_variations_nominal.pdf'.format(args.outdir))

    if args.plot_ppd:
        set_basic_style()
        can = make_ppd_comp_plot(nom_file, var_files, args.variable)
        can.SaveAs('{}/comp_ppds_full_ppds_variations_nominal.pdf'
                   .format(args.outdir))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make plots comparing'
                                     ' the ppds from different variations')
    parser.add_argument('plotconfig', help='JSON file containing the plot '
                        'configuration')
    parser.add_argument('-o', '--outdir', help='Output directory into which the '
                        'produced plot should be put', default='.')
    parser.add_argument('--plot-ppd', help='Make a plot where the PPDs are '
                        'overlaid against each other', action='store_true',
                        default=False)

    var_sel = parser.add_mutually_exclusive_group()
    var_sel.add_argument('--dlth', action='store_const', dest='variable',
                         const='dlth', help='PPDs are for delta_lambda theta')
    var_sel.add_argument('--dlph', action='store_const', dest='variable',
                         const='dlph', help='PPDs are for delta_lambda phi')
    var_sel.add_argument('--dkappa', action='store_const', dest='variable',
                         const='dkappa', help='PPDs are for delta_kappa')

    clargs = parser.parse_args()
    main(clargs)
