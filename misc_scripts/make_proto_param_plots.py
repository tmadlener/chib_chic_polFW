#!/usr/bin/env python
"""
Script to make the plots of the proto parameters of the simultaneous binned fits
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.ProcessLine('gErrorIgnoreLevel = 1001')
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, default_attributes
from utils.plot_decoration import YLABELS, FIX_RANGES, CTH_PLOT
from utils.misc_helpers import cond_mkdir
from utils.data_handling import list_obj
from utils.setup_plot_style import set_basic_style

# common name of all the plots against costh
# NOTE: Assuming that only one pT bin is done
BIN_BASE = '_v_costh_HX_fold'

def load_graphs(rfile, no_ratio=True):
    """
    Load all the graphs
    """
    # TODO: Make this work for more than one bin!
    graph_names = list_obj(rfile, 'TGraphAsymmErrors', BIN_BASE + '_bin_0')
    func_names = list_obj(rfile, 'TF1', BIN_BASE)

    graphs_funcs = {
        n.replace(BIN_BASE + '_bin_0', ''): rfile.Get(n) for n in graph_names
    }
    graphs_funcs.update({
        n.replace(BIN_BASE, ''): rfile.Get(n) for n in func_names
    })

    if no_ratio:
        graphs_funcs.pop('r_chic2_chic1', None)

    return graphs_funcs


def make_plot(pname, param_fg):
    """
    Make plot for one parameter graph or function and return the canvas
    """
    can = mkplot(param_fg, drawOpt='PE',
                 attr=default_attributes(open_markers=False),
                 yLabel=YLABELS.get(pname, pname),
                 yRange=FIX_RANGES.get(pname, None),
                 **CTH_PLOT)

    return can


def main(args):
    """Main"""
    graphf = r.TFile.Open(args.graphfile)
    set_basic_style()
    graphs = load_graphs(graphf, args.no_ratio)

    cond_mkdir(args.outdir)
    for name, graph in graphs.iteritems():
        can = make_plot(name, graph)
        can.SaveAs('{}/{}_v_costh.pdf'.format(args.outdir, name))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the proto '
                                     'parameters from the simultaneous fits')
    parser.add_argument('graphfile', help='File containing the proto parameter '
                        'graphs (created by simultaneous fit plotting script '
                        'with option --graphs')
    parser.add_argument('-o', '--outdir', help='output directory of the '
                        'created plots', default='.')
    parser.add_argument('--no-ratio', action='store_true', default=False,
                        help='Do not create the chic2 / chic1 ratio plot')


    clargs = parser.parse_args()
    main(clargs)
