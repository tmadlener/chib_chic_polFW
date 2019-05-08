#!/usr/bin/env python
"""
Script to make comparison plots of costh integrated fit params and costh binned
free fit params
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
r.gROOT.ProcessLine('gErrorIgnoreLevel = 1001')

import numpy as np

from utils.plot_helpers import (
    mkplot, default_colors, setup_legend, setup_latex, put_on_latex
)
from utils.misc_helpers import cond_mkdir

# fix the plotting ranges for some of the parameters for more meaningful plots
FIX_RANGES = {
    'BK_p2': [-0.004, 0.004],
    'CBmass1': [3.505, 3.5075],
    'beta1': [-5, -3],
    'alpha1': [0.4, 1.2],
    'CBmass2': [3.546, 3.556],
}

def make_comp_plot(graph, intgraph, fit=False):
    """
    Make one comparison plot and return the TCanvasWrapper
    """
    plt_attr = [
        {'color': default_colors()[0], 'marker': 20, 'size': 1.5},
        {'color': default_colors()[1], 'marker': 24, 'size': 0.5,
         'fillalpha': (default_colors()[1], 0.5)},
        {'color': default_colors()[4], 'line': 2, 'width': 2}
    ]

    graph_name = graph.GetName().replace('_v_costh', '')
    if graph_name in FIX_RANGES:
        yRange = FIX_RANGES[graph_name]
    else:
        yRange=None

    can = mkplot(graph, xRange=[0, 1], drawOpt='PE', attr=plt_attr[0:],
                 yRange=yRange, yLabel=graph_name, xLabel='|cos#vartheta^{HX}|')

    # some cosmetics
    can.SetTopMargin(0.05)
    can.SetRightMargin(0.02)
    can.SetLeftMargin(0.15)
    can.pltables[0].SetTitleOffset(2.1, 'Y')

    if intgraph is not None or fit:
        leg = setup_legend(0.65, 0.85, 0.97, 0.94)

    if intgraph is not None:
        can = mkplot(intgraph, can=can, drawOpt='sameP2', attr=plt_attr[1:],
                     leg=leg, legEntries=['cos#vartheta^{HX} int. fit'],
                     legOpt='PLEF')
        can = mkplot(intgraph, can=can, drawOpt='samePE', attr=plt_attr[1:])
    if fit:
        const = r.TF1('const', '[0]', 0, 1)
        const.SetParameter(0, np.array(graph.GetX())[0])
        fit_res = graph.Fit(const, 'NSEX0')

        can = mkplot(const, can=can, drawOpt='sameL', attr=plt_attr[2:],
                     leg=leg, legEntries=['fit to constant'], legOpt='L')

        ltx = setup_latex()
        can.add_tobject(ltx)
        put_on_latex(ltx,
                     [(0.595, 0.2,
                       'c = {:.2e} +/- {:.2e}'.format(const.GetParameter(0),
                                                      const.GetParError(0))),
                      (0.595, 0.25,
                       '#chi^{{2}} / ndf = {:.2f} / {}'.format(fit_res.Chi2(),
                                                               fit_res.Ndf()))
                     ],
                     ndc=True
        )
    can.Update()

    return can


def main(args):
    """Main"""
    graphf = r.TFile.Open(args.graphfile)
    graphs = {
        g.GetName().replace('_v_costh', ''):
        graphf.Get(g.GetName()) for g in graphf.GetListOfKeys()
    }

    if args.no_ratio:
        ratio_name = 'r_chic2_chic1'
        if ratio_name in graphs:
            del graphs[ratio_name]

    if args.intgraphfile is not None:
        intf = r.TFile.Open(args.intgraphfile)
        intgraphs = {
            g.GetName().replace('_v_costh', ''):
            intf.Get(g.GetName()) for g in intf.GetListOfKeys()
        }
    else:
        intgraphs = {g: None for g in graphs}

    cond_mkdir(args.outdir)
    for param in graphs:
        can = make_comp_plot(graphs[param], intgraphs[param], args.dofit)
        can.SaveAs(args.outdir + '/' + param + '_v_costh.pdf')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to make fit params vs '
                                     'costh plots')
    parser.add_argument('graphfile', help='File containing the graphs with the '
                        'free parameters vs costh')
    parser.add_argument('intgraphfile', help='File containing the graphs with '
                        'the free parameters from a costh integrated fit',
                        default=None, nargs='?')
    parser.add_argument('-o', '--outdir', help='output directory for the '
                        'created plots', default='.')
    parser.add_argument('-f', '--dofit', help='for each parameter fit the trend'
                        ' against costh with a constant and add the fit to the '
                        'plot', action='store_true', default=False)
    parser.add_argument('--no-ratio', action='store_true', default=False,
                        help='Do not create the ratio plot')

    clargs = parser.parse_args()
    main(clargs)
