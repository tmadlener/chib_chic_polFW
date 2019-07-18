#!/usr/bin/env python
"""
Script to make plots overlaying the data and the correction maps
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.data_handling import get_dataframe, apply_selections
from utils.hist_utils import (
    project, get_array, from_array, get_binning, find_bin, divide, hist2d
)
from utils.plot_helpers import mkplot
from utils.misc_helpers import select_bin, cond_mkdir


def get_correction_map(cmfile, use_pt=True, use_acc=False):
    """
    Get the correction map
    """
    map_name = 'fold_costh_phi_JpsiPt_JpsiRap_{}_HX'
    reco = 'acc' if use_acc else 'reco'

    if use_pt:
        gen_dist = project(cmfile.Get(map_name.format('gen')), [0, 1, 2])
        reco_dist = project(cmfile.Get(map_name.format(reco)), [0, 1, 2])
    else:
        gen_dist = project(cmfile.Get(map_name.format('gen')), [1, 0])
        reco_dist = project(cmfile.Get(map_name.format(reco)), [1, 0])

    accmap = divide(reco_dist, gen_dist)

    return accmap


def get_pt_bin(amap, pt_val):
    """
    Get the pt bin costh-phi map from the passed (3d) acceptance map
    """
    pt_bin = find_bin(get_binning(amap, 2), np.array([pt_val]))[0]
    val, err = get_array(amap), get_array(amap, errors=True)
    ctp_binning = np.array([get_binning(amap, i) for i in [0, 1]])

    return from_array(val[:, :, pt_bin], ctp_binning, errors=err[:, :, pt_bin])


def get_mask_graph(xbinning, ybinning, mask):
    """
    Get a TGraphErrors that can be drawn with option 'E5' to mark the
    bins where the mask applies in a TH2D
    """
    # TODO: error checking (e.g. dimensions of mask and binnings)
    xidx, yidx = np.where(mask)
    xcenters = 0.5 * (xbinning[xidx] + xbinning[xidx + 1])
    ycenters = 0.5 * (ybinning[yidx] + ybinning[yidx + 1])

    # TODO: Check if this actually works with a non-uniform binning
    xerrs, yerrs = 0.5 * np.diff(xbinning)[xidx], 0.5 * np.diff(ybinning)[yidx]
    return r.TGraphErrors(np.sum(mask), xcenters, ycenters, xerrs, yerrs)


def make_overlay_plot(pt_map, pt_data, **kwargs):
    """
    Plot the coverage of the pt_data onto the
    """
    amap_x, amap_y = get_binning(pt_map, 0), get_binning(pt_map, 1)
    if np.min(amap_x) == 0:
        costh = pt_data.costh_HX_fold.abs()
    else:
        costh = pt_data.costh_HX_fold

    data_dist = hist2d(costh, pt_data.phi_HX_fold,
                       x_hist_sett=(len(amap_x) - 1, amap_x),
                       y_hist_sett=(len(amap_y) - 1, amap_y))

    coverage = get_array(data_dist) > 0
    cov_graph = get_mask_graph(amap_x, amap_y, coverage)

    can = mkplot(pt_map, drawOpt='colz', **kwargs)
    mkplot(cov_graph, can=can, drawOpt='sameE5',
           attr=[{'color': r.kRed, 'fillalpha': (r.kRed, 0), 'marker': 1}])

    return can


def main(args):
    """Main"""
    data = get_dataframe(args.datafile)
    cmfile = r.TFile.Open(args.corrmapfile)
    accmap = get_correction_map(cmfile, True, False)

    cond_mkdir(args.outdir)

    if isinstance(accmap, r.TH2):
        plot = make_overlay_plot(accmap, data)
        plot.SaveAs('{}/corrmap_data_overlay_2d.pdf'.format(args.outdir))
    else:
        pt_binning = get_binning(accmap, 2)
        pt_bins = zip(pt_binning[:-1], pt_binning[1:])

        for pt_bin in pt_bins:
            pdata = apply_selections(data, select_bin('JpsiPt', *pt_bin))
            pmap = get_pt_bin(accmap, 0.5 * np.sum(pt_bin))
            plot = make_overlay_plot(pmap, pdata)
            plot.SaveAs('{}/corrmap_data_overlay_2d_{}_{}.pdf'
                        .format(args.outdir, int(pt_bin[0]), int(pt_bin[1])))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make some debug '
                                     'plots overlaying the data and the '
                                     'correction maps')
    parser.add_argument('datafile', help='Data file containing the data')
    parser.add_argument('corrmapfile', help='File containing the correction map')
    parser.add_argument('-o', '--outdir', help='Directory into which the plots '
                        'are stored', default='.')

    clargs = parser.parse_args()
    main(clargs)
