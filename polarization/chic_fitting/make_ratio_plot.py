#!/usr/bin/env python
"""
Make the ratio plots combining the fit results from different costh bins
"""

import numpy as np
import re

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from os.path import dirname

from utils.misc_helpers import cond_mkdir, get_storable_name, get_bin_cut_df
from utils.plot_helpers import plot_on_canvas, default_colors
from utils.roofit_utils import get_var_err
from utils.data_handling import get_dataframe
from utils.setup_plot_style import set_TDR_style
from utils.graph_utils import scale_graph

from common_func import get_bin_sel_info


mc_attributes = [
    {'color': default_colors()[0], 'marker': 25, 'size': 1.5,
     'fillalpha': (default_colors()[0], 0.5)},
    {'color': default_colors()[1], 'marker': 26, 'size': 1.5,
     'fillalpha': (default_colors()[1], 0.5)},
    {'color': default_colors()[2], 'marker': 27, 'size': 1.5,
     'fillalpha': (default_colors()[2], 0.5)}
]

data_attributes = [
    {'color': 1, 'marker': 20, 'size': 1.5}
]


def create_graph(ratios, costh_bins, costh_means):
    """
    Create graph
    """
    n_bins = len(costh_bins)
    if n_bins != len(ratios[0]) or n_bins != len(costh_means):
        logging.error('Ratios, costh_bins and costh_means must have same length'
                      ' are: {}, {}, {}'.format(len(ratios), len(costh_means),
                                                n_bins))
        return r.TGraph()

    # y_vals = np.array([v[0] for v in ratios])
    # y_err = np.array([v[1] for v in ratios])
    y_vals = ratios[0]
    y_err = ratios[1]

    x_lo = np.array([costh_means[i] - costh_bins[i][0] for i in xrange(n_bins)])
    x_hi = np.array([costh_bins[i][1] - costh_means[i] for i in xrange(n_bins)])

    graph = r.TGraphAsymmErrors(n_bins, np.array(costh_means), y_vals,
                                x_lo, x_hi, y_err, y_err)

    return graph


def get_var_in_bin(wsp, costh_bin, varname):
    """
    Get the variable from the fit results in a costh bin
    """
    wsp.loadSnapshot('snap_costh_bin_{}'.format(costh_bin))
    return get_var_err(wsp, varname) # NOTE: breaks if var is dependent


def get_ratio(num_vals, denom_vals):
    """
    Get the ratio (including uncertainties) from the numerator and denominator
    values

    Inputs have to have the same length (not checked) and be list of tuples
    with the values as first and the errors as second element
    """
    num_c = np.array([v[0] for v in num_vals])
    num_e = np.array([v[1] for v in num_vals])
    denom_c = np.array([v[0] for v in denom_vals])
    denom_e = np.array([v[1] for v in denom_vals])

    ratio = num_c / denom_c
    ratio_err = ratio * (np.sqrt(denom_e**2 / denom_c**2 + num_e**2 / num_c**2))

    return ratio, ratio_err


def get_ratio_graph(num_ws, num_var, denom_ws, denom_var, costh_bins,
                    costh_means):
    """
    Get the ratio graph using the two different variables from two (possibly)
    different workspaces
    """
    num_vals = []
    denom_vals = []
    for ibin, _ in enumerate(costh_bins):
        num_vals.append(get_var_in_bin(num_ws, ibin, num_var))
        denom_vals.append(get_var_in_bin(denom_ws, ibin, denom_var))

    ratio = get_ratio(num_vals, denom_vals)

    return create_graph(ratio, costh_bins, costh_means)


def get_outdir(outdir, fitfile):
    """
    Get the output directory (handle default case storing it to fitfile dir)
    Also creates it if not already present
    """
    if not outdir:
        outdir = dirname(fitfile)
    cond_mkdir(outdir)
    return outdir


def get_pt_range(sel_string):
    """
    Get the Jpsi pt range from the selection string
    """
    # print sel_string
    # NOTE: currently only handles integer bin boundaries
    pt_rgx = r'JpsiPt\s?>\s?((\d*\.?)\d*)\s?\&{2}\s?JpsiPt\s?<\s?((\d*\.?)\d*)'
    match = re.search(pt_rgx, sel_string)
    if match:
        return float(match.group(1)), float(match.group(3))

    logging.error('Could not get pt bins from {}'.format(sel_string))
    return -1, -1


def get_mc_vals(dfr, state, costh_bins, pol=None):
    """
    Get the number of events in each costh_bin for a given pol scenario
    """
    if pol is None:
        weight = state
    else:
        weight = get_storable_name('wPol_HX_lth_{:.2f}'.format(pol))

    state_sel = dfr[state] == 1 # select the events for a given state

    vals = []
    for ctbin in costh_bins:
        costh_sel = get_bin_cut_df(dfr, lambda d: np.abs(d.costh_HX), *ctbin)
        bin_events = dfr[costh_sel & state_sel][weight].sum()
        vals.append((bin_events, np.sqrt(bin_events),))

    return vals


def get_mc_ratio_graph(num_df, denom_df, num_sel, denom_sel,
                       costh_bins, costh_means, num_pol=None, denom_pol=None):
    """
    Get the MC ratio graph using data from dataframes

    *_sel is the name of the weight that identifies the state (i.e. only 1 and 0)
    *_pol will get converted to the correct pol weight internally
    """
    num_vals = get_mc_vals(num_df, num_sel, costh_bins, num_pol)
    denom_vals = get_mc_vals(denom_df, denom_sel, costh_bins, denom_pol)

    ratio = get_ratio(num_vals, denom_vals)

    return create_graph(ratio, costh_bins, costh_means)


def setup_legend():
    """
    Setup the legend
    """
    leg = r.TLegend(0.2, 0.15, 0.65, 0.3)
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetEntrySeparation(0.01)
    leg.SetBorderSize(0)

    return leg


def fit_const(graph):
    """
    fit a constant to the graph and return it (incl. uncertainties)
    """
    const_func = r.TF1('const_func', '[0]', 0, 1)
    fit_res = graph.Fit(const_func, 'EX00S')

    if int(fit_res) != 0:
        logging.error('Problem in fitting constant to graph')
        return -1, -1

    return fit_res.Parameter(0), fit_res.Error(0)


def rescale_MC_graphs(data_graph, mc_graphs):
    """
    rescale the MC graphs for an easier comparison of shapes

    Returns rescaled graphs
    """
    data_c = fit_const(data_graph)[0]
    r_graphs = []
    for graph in mc_graphs:
        mc_const = fit_const(graph)[0]
        logging.debug('Rescaling graph by {}'.format(data_c / mc_const))
        r_graphs.append(scale_graph(graph, data_c / mc_const))

    return r_graphs


def make_plot(data_graph, mc_graphs, pdfname, canvas_sett):
    """
    Make plot and save
    """
    mc_attributes = [
        {'color': default_colors()[0], 'marker': 25, 'size': 1.5,
         'fillalpha': (default_colors()[0], 0.5)},
        {'color': default_colors()[1], 'marker': 26, 'size': 1.5,
         'fillalpha': (default_colors()[1], 0.5)},
        {'color': default_colors()[2], 'marker': 27, 'size': 1.5,
         'fillalpha': (default_colors()[2], 0.5)}
    ]

    data_attributes = [
        {'color': 1, 'marker': 20, 'size': 1.5}
    ]

    set_TDR_style()
    can = r.TCanvas('rcan', 'rcan', 50, 50, 600, 600)
    frame = can.DrawFrame(*canvas_sett['range'])
    frame.SetXTitle(canvas_sett['xtitle'])
    frame.SetYTitle(canvas_sett['ytitle'])

    leg = setup_legend()

    mc_rescaled = rescale_MC_graphs(data_graph, mc_graphs.values())

    plot_on_canvas(can, mc_rescaled, drawOpt='P2', leg=leg, legOpt='PF',
                   legEntries=mc_graphs.keys(), attr=mc_attributes)
    plot_on_canvas(can, [data_graph], drawOpt='samePE', leg=leg,
                   legEntries=['data'], attr=data_attributes)

    can.SaveAs(pdfname)



def do_chic_ratio(args):
    """
    Make the chic ratio plot
    """
    ffile = r.TFile.Open(args.fitfile)
    ws = ffile.Get('ws_mass_fit')
    outdir = get_outdir(args.outdir, args.fitfile)

    bin_sel_info = get_bin_sel_info(args.pklfile, args.fitfile)
    costh_bins = bin_sel_info['costh_bins']
    costh_means = bin_sel_info['costh_means']

    graph = get_ratio_graph(ws, 'Nchic2', ws, 'Nchic1', costh_bins, costh_means)

    ptmin, ptmax = get_pt_range(bin_sel_info['basic_sel'])

    mcdf = get_dataframe(args.mcfile, 'chic_mc_tuple')
    mc_sel = (np.abs(mcdf.JpsiRap) < 1.2) & ((mcdf.trigger & 2) == 2) & \
             (mcdf.vtxProb > 0.01) & \
             (get_bin_cut_df(mcdf, 'chicPt', 0, 990)) & \
             (get_bin_cut_df(mcdf, 'JpsiPt', ptmin, ptmax))

    unpol = get_mc_ratio_graph(mcdf[mc_sel], mcdf[mc_sel], 'wChic2', 'wChic1',
                               costh_bins, costh_means, None, None)
    nrqcd = get_mc_ratio_graph(mcdf[mc_sel], mcdf[mc_sel], 'wChic2', 'wChic1',
                               costh_bins, costh_means, -0.3, 0.5)
    extreme = get_mc_ratio_graph(mcdf[mc_sel], mcdf[mc_sel], 'wChic2', 'wChic1',
                                 costh_bins, costh_means, -0.6, 1)

    mc_graphs = {'unpol': unpol, 'nrqcd': nrqcd, 'extreme': extreme}

    plot_sett = {
        'range': [0, 0, 1, 0.75],
        'xtitle': '|cos#theta^{HX}|', 'ytitle': '#chi_{c2} / #chi_{c1}'
    }

    make_plot(graph, mc_graphs, 'test.pdf', plot_sett)


def do_jpsi_ratios(args):
    """
    Make the jpsi ratio plots
    """
    chic_ff = r.TFile.Open(args.chic_ff)
    chic_ws = chic_ff.Get('ws_mass_fit')
    jpsi_ff = r.TFile.Open(args.jpsi_ff)
    jpsi_ws = jpsi_ff.Get('ws_mass_fit')
    outdir = get_outdir(args.outdir, args.chic_ff)

    bin_info_chic = get_bin_sel_info(args.pklfile_chic, args.chic_ff)
    costh_bins = bin_info_chic['costh_bins']
    costh_means = bin_info_chic['costh_means']
    # Not yet needed
    # bin_info_jpsi = get_bin_sel_info(args.pklfile_jpsi, args.jpsi_ff)

    chic1_graph = get_ratio_graph(chic_ws, 'Nchic1', jpsi_ws, 'Njpsi',
                                  costh_bins, costh_means)
    chic2_graph = get_ratio_graph(chic_ws, 'Nchic2', jpsi_ws, 'Njpsi',
                                  costh_bins, costh_means)

    can = r.TCanvas('rcan', 'rcan', 50, 50, 600, 600)
    frame = can.DrawFrame(0, 0, 1, 0.1)


    mc_graphs = {
        # TODO
    }
    plot_sett = {
        'range': [0, 0, 1, 0.1],
        'xtitle': '|cos#theta^{HX}|', 'ytitle': '#chi_{c1} / J/#psi'
    }
    make_plot(chic1_graph, mc_graphs, 'test_chic1.pdf', plot_sett)

    plot_sett['range'] = [0, 0, 1, 0.05]
    plot_sett['ytitle'] = '#chi_{c2} / J/#psi'
    make_plot(chic2_graph, mc_graphs, 'test_chic2.pdf', plot_sett)


def add_chic_parser(parsers):
    """
    Add the chic ratio mode subparser to the subparsers
    """
    chic_r_parser = parsers.add_parser('chic', description='Produce '
                                          'chic2 / chic1 ratio plots')
    chic_r_parser.add_argument('fitfile', help='file containing the workspace '
                               'with the chic fit results')
    chic_r_parser.add_argument('mcfile', help='mc file containing flat tuple '
                               'and all weights for the desired polarization '
                               'scenarios')
    chic_r_parser.add_argument('-pf', '--pklfile', help='Pickle file containing'
                               ' the costh binning and selection information. '
                               'Use this to override the default which derives '
                               'the name from the fitfile', default='',
                               type=str)
    chic_r_parser.add_argument('-o', '--outdir', help='Directory to which the '
                               'plots get stored (defaults to directory of '
                               'fit file)',
                               type=str, default='')
    chic_r_parser.set_defaults(func=do_chic_ratio)


def add_jpsi_parser(parsers):
    """
    Add the jpsi ratio mode subparser to the subparsers
    """
    jpsi_r_parser = parsers.add_parser('jpsi', description='Produce '
                                          'chicJ / Jpsi ratio plots')
    jpsi_r_parser.add_argument('chic_ff', help='file containing the workspace '
                               'with the chic fit results')
    jpsi_r_parser.add_argument('jpsi_ff', help='file containing the workspace '
                               'with the jpsi fit results')
    # TODO: add MC files

    jpsi_r_parser.add_argument('-pfc', '--pklfile-chic', help='Pickle file '
                               'containing the costh binning and selection '
                               'information for the chic.')
    jpsi_r_parser.add_argument('-pfj', '--pklfile-jpsi', help='Pickle file '
                               'containing the costh binning and selection '
                               'information for the jpsi.')
    jpsi_r_parser.add_argument('-o', '--outdir', help='Directory to which the '
                               'plots get stored (defaults to directory of '
                               'the chic fit file)',
                               type=str, default='')
    jpsi_r_parser.set_defaults(func=do_jpsi_ratios)


if __name__ == '__main__':
    import argparse

    main_parser = argparse.ArgumentParser(description='Script to produce ratio '
                                          'plots')

    subparsers = main_parser.add_subparsers(help='Mode to run', dest='mode')

    add_chic_parser(subparsers)
    add_jpsi_parser(subparsers)

    clargs = main_parser.parse_args()
    clargs.func(clargs)
