#!/usr/bin/env python
"""
Script to make the plots showing the delta parameters vs. pt
"""

from collections import OrderedDict

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, setup_legend, default_colors
from utils.plot_decoration import YLABELS
from utils.data_handling import list_obj
from utils.graph_utils import get_errors
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.misc_helpers import cond_mkdir

PCOLOR = default_colors()[0]

ATTR = {
    'central': [{'color': PCOLOR, 'marker': 20, 'size': 1.0, 'linewidth': 2}],
    '1': [{'fillalpha': (PCOLOR, 0.65), 'size': 0, 'linewidth': 0}],
    '2': [{'fillalpha': (PCOLOR, 0.35), 'size': 0, 'linewidth': 0}],
    '3': [{'fillalpha': (PCOLOR, 0.25), 'size': 0, 'linewidth': 0}],
}

YRANGES = {
    'dlth': [-2, 2],
    'dlph': [-0.5, 0.5]
}
PTLABEL = 'p_{T}^{J/#psi} / GeV'

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


def main(args):
    """Main"""
    set_TDR_style()
    graphfile = r.TFile.Open(args.graphfile)
    var, graphs = get_graphs(graphfile)
    can = make_plot(graphs, var, args.sigmas.split(','))
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

    clargs = parser.parse_args()
    main(clargs)

