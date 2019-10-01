#!/usr/bin/env python
"""
Script to make ratio plots with more than one variable on it
"""

from collections import OrderedDict

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, setup_legend
from utils.plot_decoration import YLABELS
from utils.graph_utils import scale_graph_x, shift_graph_horizontal
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.constants import m_psiPDG

from make_pt_dep_plot import get_combined_graph

YRANGE = [-2, 2]
XRANGE = [2, 12]
ATTR = {
    'dlth': {'color': 1, 'marker': 20, 'size': 1, 'line': 1, 'width': 2},
    'lth': {'color': 1, 'marker': 21, 'size': 1},
    'lth2': {'color': 1, 'marker': 24, 'size': 1, 'line': 1, 'width': 2}
}
PTMLABEL = '(p_{T}/M)^{J/#psi}'


def get_graph(graphfile, variable):
    """
    Get the graph from the file
    """
    if variable in ['dlth', 'lth2']:
        return graphfile.Get('{}_v_pt_n_sig_1'.format(variable))
    else:
        return graphfile.Get('{}_v_pt_central'.format(variable))


def make_plot(lth_g, dlth_g, lth2_g):
    """
    Make the plot as a function of pT/M according to Carlos specifications
    """
    lth_g = scale_graph_x(lth_g, 1 / m_psiPDG)
    dlth_g = scale_graph_x(dlth_g, 1 / m_psiPDG)
    lth2_g = scale_graph_x(lth2_g, 1 / m_psiPDG)
    lth2_g = shift_graph_horizontal(lth2_g, 0.1)

    leg = setup_legend(0.675, 0.74, 0.88, 0.88)
    can = mkplot([lth_g], drawOpt='P', attr=[ATTR['lth']],
                 xRange=XRANGE, xLabel=PTMLABEL,
                 yRange=YRANGE, yLabel='#lambda',
                 leg=leg, legEntries=[YLABELS['lth']])
    mkplot([dlth_g, lth2_g], drawOpt='samePE', can=can,
           attr=[ATTR['dlth'], ATTR['lth2']],
           leg=leg, legEntries=[YLABELS['dlth'], YLABELS['lth2']])

    add_auxiliary_info(can, 2012, prelim=True, pos='left')

    return can

def main(args):
    """Main"""
    set_TDR_style()
    graphfile = r.TFile.Open(args.graphfile)

    g_lth = get_graph(graphfile, 'lth')
    g_dlth = get_graph(graphfile, 'dlth')
    g_lth2 = get_graph(graphfile, 'lth2')

    if args.systematics is not None:
        syst_file = r.TFile.Open(args.systematics)
        syst_graph = syst_file.Get('dlth_v_pt_syst')
        g_dlth = get_combined_graph(g_dlth, syst_graph)
        g_lth2 = get_combined_graph(g_lth2, syst_graph)

    can = make_plot(g_lth, g_dlth, g_lth2)
    can.SaveAs(args.output)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make pt-dep plots '
                                     'with more than one variable int he plot')
    parser.add_argument('graphfile', help='File containing all necessary graphs')
    parser.add_argument('-o', '--output', help='Name of the output file that is '
                        'created.', default='lambdas_v_ptM.pdf')
    parser.add_argument('-s', '--systematics', help='File containing a graph '
                        'with the systematic uncertainties that should be '
                        'combined with the statistical uncertainties',
                        default=None)


    clargs = parser.parse_args()
    main(clargs)
