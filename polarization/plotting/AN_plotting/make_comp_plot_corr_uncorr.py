#!/usr/bin/env python
"""
Script to make a comparison plot of the corrected and uncorrected ratio graphs
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, default_attributes, setup_legend
from utils.setup_plot_style import add_auxiliary_info, set_TDR_style
from utils.graph_utils import shift_graph_horizontal
from utils.plot_decoration import YLABELS, VAR_PLOT
from utils.misc_helpers import cond_mkdir

from common_utils import get_graph, rescale_graph

# absolute value of the horizontal shift for visibility reasons
HORIZONTAL_SHIFT = {'costh': 0.015, 'phi':  2}
# Legend position
LEGEND_POS = {'costh': (0.675, 0.6, 0.88, 0.72),
              'phi': (0.675, 0.3, 0.88, 0.42)}

def make_ratio_plot(uncorr_ratio, corr_ratio, variable):
    """
    Make the comparison plot of corrected and uncorrected graphs
    """
    leg = setup_legend(*LEGEND_POS[variable])
    leg.SetTextSize(0.03)

    can = mkplot(corr_ratio, drawOpt='PE',
                 attr=[{'color': 1, 'size': 1, 'marker': 20}],
                 yLabel=YLABELS.get('r_chic2_chic1'), yRange=[0, 0.75],
                 leg=leg, legEntries=['corrected'],
                 **VAR_PLOT[variable])

    shift = HORIZONTAL_SHIFT[variable]
    mkplot(
        [shift_graph_horizontal(rescale_graph(uncorr_ratio, corr_ratio), shift),
         shift_graph_horizontal(uncorr_ratio, -shift)],
        can=can, drawOpt='samePE', attr=default_attributes(open_markers=True),
        leg=leg, legEntries=['fitted (scaled)', 'fitted'])

    add_auxiliary_info(can, 2012, prelim=True)

    return can


def main(args):
    """Main"""
    set_TDR_style()

    ucfile = r.TFile.Open(args.uncorrfile)
    uncorr_g = get_graph(ucfile, 'r_chic2_chic1', args.variable)
    cfile = r.TFile.Open(args.corrfile)
    corr_g = get_graph(cfile, 'r_chic2_chic1', args.variable)

    can = make_ratio_plot(uncorr_g, corr_g, args.variable)
    cond_mkdir(args.outdir)

    can.SaveAs('{}/r_chic2_chic1_v_{}_comp_corr_uncorr.pdf'
               .format(args.outdir, args.variable))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make a comparison '
                                     'plot of the corrected and uncorrected '
                                     'ratio graphs')
    parser.add_argument('uncorrfile', help='File containing the uncorrected '
                        'ratio graph')
    parser.add_argument('corrfile', help='File containing the corrected ratio '
                        'graph')
    parser.add_argument('-o', '--outdir', help='Directory into which the '
                        'produced plot is placed', default='.')

    var_sel = parser.add_mutually_exclusive_group()
    var_sel.add_argument('--costh', action='store_const', dest='variable',
                         const='costh', help='Input ratios are vs costh')
    var_sel.add_argument('--phi', action='store_const', dest='variable',
                         const='phi', help='Input ratios are vs phi')
    parser.set_defaults(variable='costh')


    clargs = parser.parse_args()
    main(clargs)
