#!/usr/bin/env python
"""
Make the plots of the fit results for each costh bin
"""

import re
import pickle
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
from os.path import dirname

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.chic_fitting import make_mass_fit_plot, print_fit_results
from utils.hist_utils import combine_cuts
from utils.misc_helpers import (
    get_bin_cut_root, get_storable_name, cond_mkdir, get_bin_cut_df
)
from utils.plot_helpers import plot_on_canvas, default_colors
from utils.data_handling import get_dataframe
from utils.roofit_utils import get_var
from utils.setup_plot_style import set_TDR_style
from utils.graph_utils import scale_graph


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

def get_bin_sel_info(pklfile, fitfile):
    """
    Get the binning and selection info
    """
    if not pklfile:
        pklfile = fitfile.replace('.root', '_bin_sel_info.pkl')

    with open(pklfile, 'r') as pklf:
        bin_sel_info = pickle.load(pklf)

    return bin_sel_info


def make_fit_res_plots(wsp, costh_bins, base_sel, outdir):
    """
    Make the plots with the fit results
    """
    for i, ctbin in enumerate(costh_bins):
        costh_cut = get_bin_cut_root('TMath::Abs(costh_HX)', *ctbin)
        snapname = 'snap_costh_bin_{}'.format(i)
        wsp.loadSnapshot(snapname)
        full_selection = combine_cuts([costh_cut, base_sel])

        pdfname = '/'.join([outdir, 'mass_fit_chic_costh_bin_{}.pdf'.format(i)])

        make_mass_fit_plot(wsp, pdfname, 'chicMass', full_selection)
        print_fit_results(wsp, pdfname.replace('.pdf', '_fit_res.pdf'))


def get_ratio_in_bin(wsp, costh_bin):
    """
    Get the fitted ratio of for the given bin
    """
    savename = 'costh_bin_{}'.format(costh_bin)
    wsp.loadSnapshot('snap_' + savename)

    ratio = get_var(wsp, 'r_chic2_chic1')
    return (ratio.getVal(), ratio.getErrorLo(), ratio.getErrorHi())


def create_graph(ratios, costh_bins, costh_means):
    """
    Create graph
    """
    n_bins = len(costh_bins)
    if n_bins != len(ratios) or n_bins != len(costh_means):
        logging.error('Ratios, costh_bins and costh_means must have same length'
                      ' are: {}, {}, {}'.format(len(ratios), len(costh_means),
                                                n_bins))
        return r.TGraph()

    y_vals = np.array([v[0] for v in ratios])
    y_lo = np.array([-v[1] for v in ratios])
    y_hi = np.array([v[2] for v in ratios])

    x_lo = np.array([costh_means[i] - costh_bins[i][0] for i in xrange(n_bins)])
    x_hi = np.array([costh_bins[i][1] - costh_means[i] for i in xrange(n_bins)])

    graph = r.TGraphAsymmErrors(len(ratios), np.array(costh_means), y_vals,
                                x_lo, x_hi, y_lo, y_hi)

    return graph


def get_data_graph(wsp, costh_bins, costh_means):
    """
    Get the data graph using the fit results
    """
    ratios = []
    for ibin, _ in enumerate(costh_bins):
        ratios.append(get_ratio_in_bin(wsp, ibin))

    return create_graph(ratios, costh_bins, costh_means)


def get_sum_var_bins(wsp, var, costh_bins):
    """
    Get the sum of a given variable in all costh_bins
    """
    var_sum = 0
    for ibin, _ in enumerate(costh_bins):
        wsp.loadSnapshot('snap_costh_bin_{}'.format(ibin))
        var_sum += get_var(wsp, var).getVal()

    return var_sum


def get_corr_weight(mcdf, wsp, costh_bins):
    """
    Get the data - mc correction weight to ensure the fraction of chic2 / chic1
    events is the same as on data
    """
    d_N1 = get_sum_var_bins(wsp, 'Nchic1', costh_bins)
    d_N2 = get_sum_var_bins(wsp, 'Nchic2', costh_bins)

    mc_N1 = mcdf.wChic1.sum()
    mc_N2 = mcdf.wChic2.sum()

    corr_w = d_N2 / d_N1 * mc_N1 / mc_N2

    print('Corr weight: N1_data, N2_data, r_data, N1_mc, N2_mc, r_mc, corr_w')
    print('{:.1f}, {:.1f}, {:.3f}, {}, {}, {:.3f}, {:.3f}'
          .format(d_N1, d_N2, d_N2 / d_N1, mc_N1, mc_N2, mc_N2 / mc_N1, corr_w))

    return corr_w


def get_mc_ratios(mcdf, pol_scen, costh_bins, corr_weight=1):
    """
    Get the MC ratios from the passed dataframe for a given polarization
    scenario

    corr_weight is used to ensure that the fraction of chic2 events in all
    signal events is the same as on data:

    --> N_1_mc / N_2_mc == N_1_data / N_2_data
    """
    if pol_scen is None:
        w_chic1, w_chic2 = 'wChic1', 'wChic2'
    else:
        # NOTE: assuming 1D reweighting in HX frame here
        w_chic1 = get_storable_name('wPol_HX_lth_{:.2f}'
                                                 .format(pol_scen[0]))
        w_chic2 = get_storable_name('wPol_HX_lth_{:.2f}'
                                                  .format(pol_scen[1]))

    ratios = []
    for ctbin in costh_bins:
        costh_sel = get_bin_cut_df(mcdf, lambda d: np.abs(d.costh_HX), *ctbin)
        chic1_sel = mcdf.wChic1 == 1
        chic2_sel = mcdf.wChic2 == 1

        n_chic1 = mcdf[costh_sel & chic1_sel][w_chic1].sum()
        n_chic2 = mcdf[costh_sel & chic2_sel][w_chic2].sum() * corr_weight

        ratio = n_chic2 / n_chic1
        ratio_err = ratio * np.sqrt(1.0 / n_chic1 + 1.0 / n_chic2)

        # print n_chic1, n_chic2, ratio, ratio_err
        ratios.append((ratio, -ratio_err, ratio_err))

    return ratios


def get_mc_graph(mcdf, pol_scen, costh_bins, costh_means, corr_weight):
    """
    Get the MC graph for a given polarization scenario
    """
    ratios = get_mc_ratios(mcdf, pol_scen, costh_bins, corr_weight)
    return create_graph(ratios, costh_bins, costh_means)


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


def analytical_ratio(name):
    """Get a TF1 with the analytical description"""
    ffunc = r.TF1(name,
                  '[0] * (3 + [1]) / (3 + [2]) * '
                  '(1 + [2] * x[0]*x[0]) / (1 + [1] * x[0]*x[0])',
                  0, 1)

    return ffunc


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


def main(args):
    """Main"""
    ffile = r.TFile.Open(args.fitfile)
    ws = ffile.Get('ws_mass_fit')

    bin_sel_info = get_bin_sel_info(args.pklfile, args.fitfile)
    costh_bins = bin_sel_info['costh_bins']
    costh_means = bin_sel_info['costh_means']

    outdir = args.outdir
    if not outdir:
        outdir = dirname(args.fitfile)
    cond_mkdir(outdir)

    if not args.no_fitresults:
        make_fit_res_plots(ws, costh_bins, bin_sel_info['basic_sel'], outdir)

    dgraph = get_data_graph(ws, costh_bins, costh_means)

    ptmin, ptmax = get_pt_range(bin_sel_info['basic_sel'])

    mcdf = get_dataframe(args.mcfile, 'chic_mc_tuple')
    mc_sel = (np.abs(mcdf.JpsiRap) < 1.2) & ((mcdf.trigger & 2) == 2) & \
             (mcdf.vtxProb > 0.01) & \
             (get_bin_cut_df(mcdf, 'chicPt', 0, 990)) & \
             (get_bin_cut_df(mcdf, 'JpsiPt', ptmin, ptmax))

    data_mc_cw = get_corr_weight(mcdf[mc_sel], ws, costh_bins)

    unpol = get_mc_graph(mcdf[mc_sel], None,
                         costh_bins, costh_means, data_mc_cw)
    nrqcd = get_mc_graph(mcdf[mc_sel], (0.5, -0.3),
                         costh_bins, costh_means, data_mc_cw)
    extreme = get_mc_graph(mcdf[mc_sel], (1, -0.6),
                           costh_bins, costh_means, data_mc_cw)

    scaled_graphs = rescale_MC_graphs(dgraph, [unpol, nrqcd, extreme])

    set_TDR_style()
    can = r.TCanvas('rcan', 'rcan', 50, 50, 600, 600)
    frame = can.DrawFrame(0, 0, 1, 0.65)
    frame.SetXTitle('|cos#theta^{HX}|')
    frame.SetYTitle('#chi_{c2} / #chi_{c1}')

    leg = setup_legend()

    plot_on_canvas(can, scaled_graphs, drawOpt='sameP2',
                   attr=mc_attributes, leg=leg, legOpt='PF',
                   legEntries=['unpolarized', 'nrqcd', 'extreme'])
    plot_on_canvas(can, [dgraph], drawOpt='PE', attr=data_attributes, leg=leg,
                   legEntries=['data'])

    can.SaveAs('{}/chic2_chic1_pt{}_{}_nbins{}_costh_HX.pdf'
               .format(outdir, ptmin, ptmax, len(costh_bins)))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to produce the fit '
                                     'of the chic mass')
    parser.add_argument('fitfile', help='file containing the workspace with '
                        'the fit results and the dataset used for fitting')
    parser.add_argument('mcfile', help='mc file containing flat tuple and '
                        'all weights for the desired polarization scenarios')
    parser.add_argument('-pf', '--pklfile', help='Pickle file containing the '
                        'costh binning and selection informations. Use this to '
                        'override the default which derives the name from the '
                        'fitfile', default='', type=str)
    parser.add_argument('-o', '--outdir', help='Directory to which the plots '
                        'get stored (defaults to same directory as fitfile)',
                        default='', type=str)
    parser.add_argument('-nf', '--no-fitresults', default=False,
                        action='store_true', help='Prevent plotting of fit '
                        'result plots')

    clargs = parser.parse_args()
    main(clargs)
