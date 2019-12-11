#!/usr/bin/env python
"""
Script to make 2D contour plots of lambda2 vs lambda1
"""

import numpy as np
from scipy.spatial import ConvexHull

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import (
    mkplot, default_attributes, setup_legend, setup_latex, put_on_latex
)
from utils.hist_utils import get_array, get_binning
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info


def get_contour(hist):
    """
    Get the outer contour of all filled points in the histogram
    """
    vals = get_array(hist) > 0
    xbinning, ybinning = get_binning(hist, 0), get_binning(hist, 1)

    xvals = 0.5 * (xbinning[:-1] + xbinning[1:])
    yvals = 0.5 * (ybinning[:-1] + ybinning[1:])

    filled = []
    for ix, xv in enumerate(xvals):
        for iy, yv in enumerate(yvals):
            if vals[ix, iy]:
                filled.append([xv, yv])
    filled = np.array(filled)

    hull = ConvexHull(filled)

    # Append the first point again at the end to "close" the contour
    xcont = filled[hull.vertices, 0]
    xcont = np.append(xcont, np.array([xcont[0]]))

    ycont = filled[hull.vertices, 1]
    ycont = np.append(ycont, np.array(ycont[0]))

    return r.TGraph(len(hull.vertices) + 1, xcont, ycont)


def get_phys_region():
    """
    Get the lines for the physically allowed region
    """
    return [
        r.TLine(1, 1, 1, -0.6),
        r.TLine(1, 1, -1./3., 1),
        r.TLine(-1./3., -0.6, 1, -0.6),
        r.TLine(-1./3., -0.6, -1./3., 1)
    ]


def get_jz_graph():
    """
    Get a graph with the markers for the J_z components
    """
    return r.TGraph(6,
                    np.array([-1./3, -1./3, -1./3, 1, 1, 1]),
                    np.array([-0.6, -1./3, 1, -0.6, -1./3, 1]))


CONT_ATTR = default_attributes(linewidth=4)
CONT_ATTR[0]['line'] = 1
CONT_ATTR[1]['line'] = 9
CONT_ATTR[2]['line'] = 2
CONT_ATTR[2]['color'] = 1


def main(args):
    """Main"""
    set_TDR_style()
    fitfile = r.TFile.Open(args.fitfile)

    cont_1sigma = get_contour(fitfile.Get('chi2_scan_contour_1sigma'))
    cont_2sigma = get_contour(fitfile.Get('chi2_scan_contour_2sigma'))
    cont_3sigma = get_contour(fitfile.Get('chi2_scan_contour_3sigma'))

    leg = setup_legend(0.5, 0.78, 0.8, 0.94)
    leg.SetTextSize(0.035)
    leg.SetEntrySeparation(0.005)
    leg.SetFillStyle(1001)

    can = mkplot([cont_1sigma, cont_2sigma, cont_3sigma], drawOpt='L',
                 xRange=[-1, 1.6], xLabel='#lambda_{#vartheta}^{#chi_{c1}}',
                 yRange=[-1.2, 1.6], yLabel='#lambda_{#vartheta}^{#chi_{c2}}',
                 attr=CONT_ATTR,
                 leg=leg, legEntries=['{} % CL'.format(v) for v in [68, 95, 99]],
                 legOpt='L')

    mkplot(get_phys_region(), can=can, drawOpt='same',
           attr=[{'color': r.kRed, 'line': 7, 'width': 3}])

    mkplot([r.TLine(0, -1.2, 0, 1.6), r.TLine(-1, 0, 1.6, 0)], can=can, drawOpt='same',
           attr=[{'color': 12, 'line': 7, 'width': 2}])

    mkplot(r.TGraph(1, np.array([0]), np.array([0])), can=can, drawOpt='Psame',
           attr=[{'color': 1, 'marker': 24, 'size': 2.75}])

    mkplot(get_jz_graph(), can=can, drawOpt='Psame',
           attr=[{'color': 1, 'marker': 30, 'size': 3}])

    ltx = setup_latex()
    ltx.SetTextSize(0.04)
    put_on_latex(ltx, [
        (-0.4, 1.3, 'J_{z}'),
        (-0.4, 1.1, '#pm1'),
        (0.975, 1.1, '0'),
        (-0.54, 0.95, '#pm2'),
        (-0.54, -0.38, '#pm1'),
        (-0.47, -0.65, '0')
    ], ndc=False)

    add_auxiliary_info(can, 2012, 'left')
    can.SaveAs(args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make 2D contour plot'
                                     ' of lambda2 vs lambda1')
    parser.add_argument('fitfile', help='File containing 2D histograms from chi2'
                        ' scans')
    parser.add_argument('-o', '--outfile', help='Name of the outputfile',
                        default='lth2_v_lth1_contours.pdf')

    clargs = parser.parse_args()
    main(clargs)
