#!/usr/bin/env python
"""
Script to print an almost valid .tex table from a TGraphAsymmErrors
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.graph_utils import get_errors

RATIO_NAME = 'r_chic2_chic1_v_{}_HX_fold_bin_0'

def get_graph(gfile, direction):
    """Get the graph from the passed file"""
    return gfile.Get(RATIO_NAME.format(direction))

def get_fmt_str(direction):
    """Get the fmt string for the direction"""
    if direction == 'costh':
        return r'{0:.3f}--{1:.3f} & {2:.4f} & {3:.4f} & ${4}$ \\'
    return r'{0:.0f}--{1:.0f} & {2:.2f} & {3:.4f} & ${4}$ \\'


def dump(graph, direction):
    """
    Dump the graph to a .tex table
    """
    x_vals = np.array(graph.GetX())
    y_vals = np.array(graph.GetY())
    xlo, xhi, ylo, yhi = get_errors(graph)

    x_lo, x_hi = x_vals - xlo, x_vals + xhi

    fmt_str = get_fmt_str(direction)

    for i, x in enumerate(x_vals):
        uncer_lo = '{:.4f}'.format(ylo[i])
        uncer_hi = '{:.4f}'.format(yhi[i])
        if uncer_lo == uncer_hi:
            uncer = '\pm {}'.format(uncer_lo)
        else:
            uncer = '{{}}^{{+{}}}_{{-{}}}'.format(uncer_hi, uncer_lo)

        print(fmt_str.format(x_lo[i], x_hi[i], x, y_vals[i], uncer))


def main(args):
    """Main"""
    rfile = r.TFile.Open(args.infile)
    graph = get_graph(rfile, args.direction)

    dump(graph, args.direction)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to generate .tex table '
                                     'stumps from TGraphAsymmErrors')
    parser.add_argument('infile', help='File containing the TGraph')

    dir_sel = parser.add_mutually_exclusive_group()
    dir_sel.add_argument('--costh', action='store_const', dest='direction',
                         const='costh',
                         help='Assume that the ratio is vs costh')
    dir_sel.add_argument('--phi', action='store_const', dest='direction',
                         const='phi', help='Assume that the ratio is vs phi')
    parser.set_defaults(direction='costh')


    clargs = parser.parse_args()
    main(clargs)
