#!/usr/bin/env python
"""
Script to determine the systematic uncertainties
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.graph_utils import (
    divide_graphs, subtract_graphs, fit_to_graph, scale_graph, get_errors
)
from utils.plot_helpers import mkplot
from utils.setup_plot_style import set_basic_style
from utils.plot_decoration import CTH_PLOT


def get_graph(rfile, param):
    """Get the ratio graph from the file"""
    return rfile.Get('{}_v_costh_HX_fold_bin_0'.format(param))


def get_scaled_graphs(var_graph, cent_graph):
    """
    Rescale the variation graphs to have the same normalization as the central
    graph, by doing a fit using the central graph as template
    """
    scales = [fit_to_graph(g, cent_graph) for g in var_graph]
    sgraphs = [scale_graph(g, 1/scales[i][0]) for i, g in enumerate(var_graph)]

    return sgraphs


def get_diff(var_graphs, cent_graph, rel=False):
    """
    Get the (relative) difference between the central graph and the variations
    """
    sub_graphs = [subtract_graphs(g, cent_graph, 1.0) for g in var_graphs]
    if rel:
        # Create a new graph without uncertainties to define the central graph
        # as "truth" graph. While the uncertainties along the x-direction would
        # in principle not change, they are anyhow dropped in the division below
        cent_graph = r.TGraphAsymmErrors(cent_graph.GetN(),
                                         np.array(cent_graph.GetX()),
                                         np.array(cent_graph.GetY()))

        sub_graphs = [divide_graphs(g, cent_graph) for g in sub_graphs]

    return sub_graphs


def get_largest_dev_graph(diff_graphs):
    """
    Get the graph containing in each point the largest deviation from 0 as
    uncertainties. The uncertainties of the input graphs are ignored for now.
    """
    y_vals = np.array([g.GetY() for g in diff_graphs])
    max_dev = np.max(y_vals, axis=0)
    min_dev = -np.min(y_vals, axis=0)

    # Assume that all of them have the same binning
    elo, ehi, _, _ = get_errors(diff_graphs[0])
    xval = np.array(diff_graphs[0].GetX())

    return r.TGraphAsymmErrors(len(elo), xval, np.zeros_like(elo), elo, ehi,
                               min_dev, max_dev)


def create_outfile(fname, param, centg, vargs, svargs, rsvargs, ldgraph):
    """
    Store graphs into one rootfile for easier retrieval later
    """
    outfile = r.TFile(fname, 'recreate')
    centg.SetName('_'.join(['central', param]))
    centg.Write()

    ldgraph.SetName('_'.join(['largest_dev', param]))
    ldgraph.Write()

    for igr, graph in enumerate(vargs):
        graph.SetName('_'.join(['var', str(igr), param]))
        graph.Write()

    for igr, graph in enumerate(svargs):
        graph.SetName('_'.join(['var', str(igr), param, 'scaled']))
        graph.Write()

    for igr, graph in enumerate(rsvargs):
        graph.SetName('_'.join(['var', str(igr), param, 'scaled_rel_diff']))
        graph.Write()

    outfile.Close()


def process_param(cfile, vfiles, param, outfile=None):
    """Process one parameter"""
    cgraph = get_graph(cfile, param)
    vgraphs = [get_graph(f, param) for f in vfiles]

    scaled_vgraphs = get_scaled_graphs(vgraphs, cgraph)
    rel_scaled_vgraphs = get_diff(scaled_vgraphs, cgraph, rel=True)

    largest_dev_graph = get_largest_dev_graph(rel_scaled_vgraphs)


    if outfile is not None:
        create_outfile(outfile, param, cgraph, vgraphs, scaled_vgraphs,
                       rel_scaled_vgraphs, largest_dev_graph)

    # for testing / development
    # can = mkplot(rel_scaled_vgraphs, drawOpt='PE', **CTH_PLOT)
    # mkplot(largest_dev_graph, drawOpt='sameE2', can=can,
    #        attr=[{'fillalpha': (1, 0.5)}])

    # can = mkplot(scaled_vgraphs, drawOpt='PE', **CTH_PLOT)
    # can = mkplot(rel_scaled_vgraphs, drawOpt='PE', **CTH_PLOT)
    # can = mkplot(get_diff(scaled_vgraphs, cgraph), drawOpt='PE', **CTH_PLOT)


    # diff_graphs = get_diff(vgraphs, cgraph)
    # reldiff_graphs = get_diff(vgraphs, cgraph, True)

    # can = mkplot(diff_graphs, drawOpt='P', **CTH_PLOT)

    # can2 = mkplot(reldiff_graphs, drawOpt='P', **CTH_PLOT)



def main(args):
    """Main"""
    set_basic_style()
    cgraph_file = r.TFile.Open(args.centralfile)
    var_files = [r.TFile.Open(f) for f in args.variationfiles if f != args.centralfile]

    for param in args.params.split(','):
        process_param(cgraph_file, var_files, param, 'all_systematic_graphs.root')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to determine systematic'
                                     'uncertainties')
    parser.add_argument('centralfile', help='File containing the central result'
                        'graph')
    parser.add_argument('variationfiles', nargs='+', help='Files containing the '
                        'graphs for the systematic variations. (If the central '
                        'results file is present in this list it will be sorted'
                        ' out)')
    parser.add_argument('-p', '--params', help='Comma separated list of '
                        'parameters for which the systematic uncertainties '
                        'should be produced', default='r_chic2_chic1')


    clargs = parser.parse_args()
    main(clargs)
