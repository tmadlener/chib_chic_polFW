#!/usr/bin/env python
"""
Script to make the plots showing the delta parameters vs. pt
"""

from collections import OrderedDict

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import (
    mkplot, setup_legend, default_colors, setup_latex, put_on_latex
)
from utils.plot_decoration import YLABELS
from utils.data_handling import list_obj
from utils.graph_utils import get_errors
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.misc_helpers import cond_mkdir, fmt_float
from utils.constants import m_psiPDG

PCOLOR = default_colors()[0]
SCOLOR = default_colors()[2]
FCOLOR = default_colors()[1]


ATTR = {
    'central': [{'color': PCOLOR, 'marker': 20, 'size': 1.0, 'width': 2}],
    '1': [{'fillalpha': (PCOLOR, 0.65), 'size': 0, 'width': 0}],
    '2': [{'fillalpha': (PCOLOR, 0.35), 'size': 0, 'width': 0}],
    '3': [{'fillalpha': (PCOLOR, 0.25), 'size': 0, 'width': 0}],
    'combined': [{'color': PCOLOR, 'marker': 20, 'size': 1.0, 'width': 0,
                  'fillalpha': (PCOLOR, 1.0)}],
    'stat only': [{'color': SCOLOR, 'marker': 20, 'size': 0, 'width': 0,
                   'fillalpha': (SCOLOR, 1.0)}],
    'fit': [{'color': FCOLOR, 'width': 2}]
}

YRANGES = {
    'dlth': [-2.0, 2.0],
    'dlph': [-0.5, 0.5]
}
PTLABEL = 'p_{T}^{J/#psi} / GeV'
PTMLABEL = '(p_{T}/M)^{J/#psi}'

LEG_KEYS = {'1': '68 % CI', '2': '95 % CI', '3': '99 % CI'}

def get_graphs(gfile):
    """Get the graphs and the variable for which the graphs are plotted"""
    graphnames = list_obj(gfile, 'TGraphAsymmErrors', r'v_pt_(n_sig|cl)_')
    graphs = OrderedDict()
    for name in graphnames:
        n_sig = name.split('_')[-1]
        graphs[n_sig] = gfile.Get(name)

    var = graphs.values()[0].GetName().split('_')[0]

    return var, graphs


def remove_y_errors(graph):
    """
    Get a (copy) of the passed graph with removed y errors
    """
    x_vals, y_vals = np.array(graph.GetX()), np.array(graph.GetY())
    x_lo, x_hi, _, _ = get_errors(graph)

    return r.TGraphAsymmErrors(len(x_vals), x_vals, y_vals, x_lo, x_hi)


def remove_x_errors(graph):
    """
    Get a (copy) of the passed graph with removed x errors
    """
    x_vals, y_vals = np.array(graph.GetX()), np.array(graph.GetY())
    _, _, y_lo, y_hi = get_errors(graph)

    nbins = len(x_vals)
    x_err = np.zeros(nbins, dtype=float)

    return r.TGraphAsymmErrors(nbins, x_vals, y_vals, x_err, x_err, y_lo, y_hi)


def make_uncer_plot(graphs, variable, levels):
    """
    Make the plot plotting the uncertainties and return the canvas
    """
    # All graphs have the same central values so it doesn't matter which one we
    # take
    cent_graph = remove_y_errors(graphs.values()[0])

    leg = setup_legend(0.675, 0.74, 0.88, 0.88)

    can = mkplot(cent_graph, drawOpt='PE', attr=ATTR['central'],
                 xRange=[8, 30], xLabel=PTLABEL,
                 yRange=YRANGES[variable], yLabel=YLABELS[variable],
                 leg=leg, legEntries=['central'], legOpt='P')

    for n_sig in levels:
        graph = graphs[n_sig]
        mkplot(graph, can=can, drawOpt='E2same', legOpt='F', attr=ATTR[n_sig],
               leg=leg, legEntries=[LEG_KEYS[n_sig]])

    mkplot(cent_graph, can=can, drawOpt='samePE', attr=ATTR['central'])

    add_auxiliary_info(can, 2012, prelim=True, pos='left')

    return can


def make_limit_plot(graphs, variable):
    """
    Make a limit plot
    """
    leg = setup_legend(0.675, 0.74, 0.88, 0.88)
    can = mkplot(graphs.values(), xRange=[8, 30], xLabel=PTLABEL, drawOpt='PE',
                 yRange=[-0.5, 1.1], yLabel=YLABELS[variable],
                 leg=leg, legOpt='L',
                 legEntries=['{} % C.L.'.format(v) for v in graphs.keys()])

    mkplot(r.TLine(8, -0.33333, 30, -0.33333), can=can, drawOpt='same',
           attr=[{'color': 12, 'line': 7, 'linewidth': 2}])

    add_auxiliary_info(can, 2012, prelim=True, pos='left')

    return can


def make_plot(graphs, variable, levels):
    """
    Make the graph and return the canvas
    """
    if variable in ['dlth', 'dlph']:
        return make_uncer_plot(graphs, variable, levels)
    else:
        return make_limit_plot(graphs, variable)


def get_combined_graph(stat, syst):
    """
    Get the combined graph of the statistical and systematical uncertainties
    (added in quadrature)
    """
    syst_uncer = np.array(syst.GetY()) # stored in y-values
    n_bins = len(syst_uncer)
    xlow, xhigh, ylow, yhigh = get_errors(stat)
    xcentral, ycentral = np.array(stat.GetX()), np.array(stat.GetY())

    ylow = np.sqrt(ylow**2 + syst_uncer**2)
    yhigh = np.sqrt(yhigh**2 + syst_uncer**2)

    return r.TGraphAsymmErrors(n_bins, xcentral, ycentral,
                               xlow, xhigh, ylow, yhigh)


def scale_graph_x(graph, scale):
    """
    Scale the graph in x-direction
    """
    x_vals, y_vals = np.array(graph.GetX()), np.array(graph.GetY())
    xlo, xhi, ylo, yhi = get_errors(graph)

    x_vals *= scale
    xlo *= scale
    xhi *= scale

    return r.TGraphAsymmErrors(len(xlo), x_vals, y_vals, xlo, xhi, ylo, yhi)


def make_plot_with_syst(graphs, variable, syst_file,
                        add_fit=False, vs_pt_m=False):
    """
    Make a plot including the systematic uncertainties.
    First combine the stat. and syst. errors by summing them in quadrature then
    plot the combined as well as the stat. only graphs onto the plot.
    Only uses the +/- 1 sigma stat. uncertainties
    """
    syst_graph = syst_file.Get('{}_v_pt_syst'.format(variable))
    stat_graph = graphs['1']
    comb_graph = get_combined_graph(stat_graph, syst_graph)

    xran = [8, 30]
    xlab = PTLABEL

    if vs_pt_m:
        xran = [v / m_psiPDG for v in xran]
        comb_graph = scale_graph_x(comb_graph, 1 / m_psiPDG)
        stat_graph = scale_graph_x(stat_graph, 1 / m_psiPDG)
        xlab = PTMLABEL

    leg = setup_legend(0.675, 0.74, 0.88, 0.88)
    can = mkplot(comb_graph, drawOpt='PE2', attr=ATTR['combined'],
                 xRange=xran, xLabel=xlab,
                 yRange=YRANGES[variable], yLabel=YLABELS[variable],
                 leg=leg, legEntries=['stat. + syst.'], legOpt='PF')

    # remove the x uncertainties to avoid having the vertical line in the plot
    # stat_graph = remove_x_errors(stat_graph)
    mkplot(stat_graph, can=can, drawOpt='PE2same', attr=ATTR['stat only'],
           leg=leg, legEntries=['stat. only'], legOpt='F')

    # stat only covers the point of the combined one. Plot them again
    mkplot(comb_graph, drawOpt='PEX0same', attr=ATTR['combined'], can=can)

    add_auxiliary_info(can, 2012, prelim=True, pos='left')

    if add_fit:
        lin_func = r.TF1('lin_func',
                         '[0] + [1] * x[0]', *xran)
        fit_res = comb_graph.Fit(lin_func, 'SEX00')
        chi2, ndf = fit_res.Chi2(), fit_res.Ndf()
        p0, p1 = [lin_func.GetParameter(i) for i in [0, 1]]

        const_func = r.TF1('const_func', '[0]', *xran)
        fit_res = comb_graph.Fit(const_func, 'SEX00')
        chi2_c, ndf_c = fit_res.Chi2(), fit_res.Ndf()
        p0c, p0ce = [const_func.GetParameter(0), const_func.GetParError(0)]

        mkplot(lin_func, drawOpt='Lsame', can=can, attr=ATTR['fit'])

        ltx = setup_latex()
        put_on_latex(ltx, [
            (0.5, 0.35, 'linear fit: {:.3f} + {:.3f} x p_{{T}}'.format(p0, p1)),
            (0.5, 0.3, '#chi^{{2}} / ndf = {:.2f} / {}'.format(chi2, ndf)),
            (0.5, 0.25, 'const. fit: {:.3f} #pm {:.3f}'.format(p0c, p0ce)),
            (0.5, 0.2, '#chi^{{2}} / ndf = {:.2f} / {}'.format(chi2_c, ndf_c))
        ], ndc=True)
        ltx.Draw()
        can.Update()


    return can


def main(args):
    """Main"""
    set_TDR_style()
    graphfile = r.TFile.Open(args.graphfile)
    var, graphs = get_graphs(graphfile)

    if args.systematics is None:
        can = make_plot(graphs, var, args.sigmas.split(','))
    else:
        syst_file = r.TFile.Open(args.systematics)
        can = make_plot_with_syst(graphs, var, syst_file, args.fit, args.pt_over_m)

    cond_mkdir(args.outdir)
    can.SaveAs('{}/{}_v_pt.pdf'.format(args.outdir, var))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make the plots of '
                                     'the graphs produced by make_pt_dep_graphs'
                                     '.py')
    parser.add_argument('graphfile', help='File containing the graphs with the '
                        '1, 2 and 3 sigma uncertainty bands')
    parser.add_argument('-o', '--outdir', help='Output directory into which the '
                        'plot will be stored', default='.')
    parser.add_argument('-n', '--sigmas', help='Comma separated list of '
                        'uncertainty levels that should be plotted. NOTE: '
                        'can only be 1, 2 or 3 (or any combination of them)',
                        default='1')
    parser.add_argument('-s', '--systematics', help='File containing a graph '
                        'with the systematic uncertainties that should be '
                        'combined with the statistical uncertainties',
                        default=None)
    parser.add_argument('-f', '--fit', help='Add a linear fit vs pt to the graph'
                        ' with the combined uncertainties', default=False,
                        action='store_true')
    parser.add_argument('-m', '--pt-over-m', help='Make the plot as a function of'
                        'pT/M instead of pT. (Using M(jpsi) for the mass)',
                        action='store_true', default=False)

    clargs = parser.parse_args()
    main(clargs)
