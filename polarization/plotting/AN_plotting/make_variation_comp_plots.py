#!/usr/bin/env python
"""
Script to make pedagogical plots of the variations of the ratio
"""

import json

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.setup_plot_style import set_basic_style
from utils.plot_decoration import YLABELS, VAR_PLOT
from utils.plot_helpers import (
    mkplot, default_attributes, setup_legend, increase_label_size
)
from utils.misc_helpers import create_random_str, cond_mkdir, flatten
from utils.graph_utils import shift_graph_horizontal, scale_graph

# Attributes for central graph
CENTRAL_ATTR = [{'color': 1, 'marker': 20, 'size': 1}]
# For largest deviation
LARGE_DEV_ATTR = [{'fillalpha': (1, 0.25)}]
# for horizontal line
LINE_ATTR = [{'color': 1, 'line': 1, 'width': 2}]

# absolute value of the horizontal shift for visibility reasons
HORIZONTAL_SHIFT = {'costh': 0.01, 'phi':  1}
# y-ranges for the ratio comparison plots
YRANGE = {
    'costh': [0.3, 0.6],
    'phi': [0.25, 0.55]
}
# y-ranges for the plots with differences w.r.t. nominal
YRANGE_DIFF = {
    'costh': [-0.05 , 0.05],
    'phi': [-0.05 , 0.05]
}

YRANGE_REL_DIFF = {
    'costh': [-5, 5],
    'phi': [-5, 5]
}

def get_variation_graphs(graph_file, base, plot_config):
    """
    Get the scaled or difference graphs from the passed graph file
    for the variations specified in the plot_config
    """
    gbasen = 'var_{{}}_r_chic2_chic1{}'.format(base)

    return [
        graph_file.Get(gbasen.format(v[0])) for v in plot_config['variations']
    ]


def create_legend(xlo, xhi, ylo, plot_config, with_nominal):
    """
    Create a legend that fits all of the variations in plot_config
    """
    n_vars = len(plot_config['variations'])
    n_vars += with_nominal

    leg = setup_legend(xlo, ylo, xhi, ylo + n_vars * 0.045)
    leg.SetTextSize(0.03)
    return leg


def shift_graphs(graphs, abs_shift, offset=1):
    """
    Shift all graphs slightly along the horizontal direction.
    If offset==0 the first graph will remain unchanged and the others will be
    shifted symmetrically around this one, otherwise all the graphs will just be
    distributed symmetrically around the 0
    """
    # Just define 17 shifts and hope that this will be enough for all variations
    shifts = list(flatten([0] + [[-abs_shift * i, abs_shift * i] for i in xrange(1, 8)]))
    return [shift_graph_horizontal(g, shifts[i + offset]) for i, g in enumerate(graphs)]


def make_ratio_comp_plot(graph_file, plot_config, variable):
    """
    Make the plot comparing the unscaled and scaled graphs in two different panels
    """
    cgraph = graph_file.Get('central_r_chic2_chic1')
    vgraphs = get_variation_graphs(graph_file, '', plot_config)
    vsgraphs = get_variation_graphs(graph_file, '_scaled', plot_config)

    can = r.TCanvas(create_random_str(), '', 600, 600)
    can.cd()
    pad1 = r.TPad('unscaled_pad', '', 0, 0.5, 1, 1)
    r.SetOwnership(pad1, False)
    pad1.Draw()
    pad1.cd()

    leg = create_legend(0.675, 0.88, 0.16, plot_config, True)
    leg.SetTextSize(leg.GetTextSize() * 2)

    pad1 = mkplot(cgraph, drawOpt='PE', can=pad1, attr=CENTRAL_ATTR,
                  yRange=YRANGE[variable], yLabel='unscaled ' + YLABELS['r_chic2_chic1'],
                  leg=leg, legEntries=['nominal'],
                  **VAR_PLOT[variable])
    mkplot(shift_graphs(vgraphs, HORIZONTAL_SHIFT[variable]),
           can=pad1, drawOpt='samePEX0',
           leg=leg, legEntries=[v[1] for v in plot_config['variations']])
    increase_label_size(pad1.pltables[0], 1.8)
    pad1.pltables[0].SetXTitle('')
    pad1.pltables[0].SetTitleOffset(0.75, 'Y')

    can.cd()
    pad2 = r.TPad('scaled_pad', '', 0, 0, 1, 0.535)
    pad2.SetBottomMargin(0.19)
    r.SetOwnership(pad2, False)
    pad2.Draw()
    pad2.cd()

    pad2 = mkplot(cgraph, drawOpt='PE', can=pad2, attr=CENTRAL_ATTR,
                  yRange=YRANGE[variable], yLabel='normalized ' + YLABELS['r_chic2_chic1'],
                  **VAR_PLOT[variable])
    mkplot(shift_graphs(vsgraphs, HORIZONTAL_SHIFT[variable]),
           can=pad2, drawOpt='samePEX0')
    increase_label_size(pad2.pltables[0], 1.8)
    pad2.pltables[0].SetTitleOffset(0.9, 'X')
    pad2.pltables[0].SetTitleOffset(0.75, 'Y')

    # Make the pads survive
    can._pads = [pad1, pad2]
    return can


def make_abs_diff_plot(graph_file, plot_config, variable):
    """
    Make the plot of the absolute differences
    """
    dgraphs = get_variation_graphs(graph_file, '_diff', plot_config)

    leg = create_legend(0.675, 0.88, 0.16, plot_config, False)

    can = mkplot(shift_graphs(dgraphs, HORIZONTAL_SHIFT[variable], 0),
                 drawOpt='PLEX0',
                 yRange=YRANGE_DIFF[variable], yLabel='difference w.r.t. nominal',
                 leg=leg, legEntries=[v[1] for v in plot_config['variations']],
                 **VAR_PLOT[variable])

    return can


def make_rel_diff_plot(graph_file, plot_config, variable):
    """
    Make the plot of the relative, scaled differences
    """
    rgraphs = get_variation_graphs(graph_file, '_scaled_rel_diff', plot_config)
    largest_d_graph = graph_file.Get('largest_dev_r_chic2_chic1')

    # To have the y-axis in percent
    rgraphs = [scale_graph(g, 100) for g in rgraphs]
    largest_d_graph = scale_graph(largest_d_graph, 100)

    leg = create_legend(0.675, 0.88, 0.16, plot_config, True)

    can = mkplot(shift_graphs(rgraphs, HORIZONTAL_SHIFT[variable], 0),
                 drawOpt='PLEX0', yRange=YRANGE_REL_DIFF[variable],
                 yLabel='scaled relative difference w.r.t. nominal [%]',
                 leg=leg, legEntries=[v[1] for v in plot_config['variations']],
                 **VAR_PLOT[variable])
    mkplot(largest_d_graph, can=can, drawOpt='sameE2',
           attr=LARGE_DEV_ATTR, leg=leg, legEntries=['largest deviation'],
           legOpt='F')
    x_range = VAR_PLOT[variable]['xRange']
    mkplot(r.TLine(0, 0, x_range[1], 0), can=can, drawOpt='same',
           attr=LINE_ATTR)

    return can


def main(args):
    """Main"""
    set_basic_style()
    vargraph_f = r.TFile.Open(args.variationfile)

    with open(args.plot_config, 'r') as pconff:
        plot_config = json.load(pconff)

    cond_mkdir(args.outdir)

    can = make_ratio_comp_plot(vargraph_f, plot_config, args.variable)
    can.SaveAs('{}/r_chic2_chic1_comp_vars_nominal.pdf'.format(args.outdir))

    can = make_abs_diff_plot(vargraph_f, plot_config, args.variable)
    can.SaveAs('{}/r_chic2_chic1_abs_diff_comp_vars.pdf'.format(args.outdir))

    can = make_rel_diff_plot(vargraph_f, plot_config, args.variable)
    can.SaveAs('{}/r_chic2_chic1_rel_diff_comp_vars.pdf'.format(args.outdir))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make (pedagogical) '
                                     'plots of the variations of the ratio')
    parser.add_argument('variationfile', help='File containing the all the '
                        'variations of a graph (as produced by'
                        '\'make_systematic_variation_graphs.py\')')
    parser.add_argument('plot_config', help='json config file containing '
                        'information about the variations')
    parser.add_argument('-o', '--outdir', help='Output directory into which the '
                        'plots will be stored', default='.')

    var_sel = parser.add_mutually_exclusive_group()
    var_sel.add_argument('--costh', action='store_const', dest='variable',
                         const='costh', help='Input ratios are vs costh')
    var_sel.add_argument('--phi', action='store_const', dest='variable',
                         const='phi', help='Input ratios are vs phi')
    parser.set_defaults(variable='costh')



    clargs = parser.parse_args()
    main(clargs)
