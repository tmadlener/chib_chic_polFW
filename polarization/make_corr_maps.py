#!/usr/bin/env python
"""
Script to determine correction maps
"""

import sys

from collections import OrderedDict

import logging
logging.basicConfig(level=logging.WARNING,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import numpy as np
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.data_handling import get_dataframe
from utils.hist_utils import (
    divide, hist2d, hist3d, get_array, from_array, get_binning
)
from utils.misc_helpers import _get_var, parse_binning

WEIGHT_F = lambda d: d.gamma_eff_sm * 0.01 * d.lepP_eff_sm * d.lepN_eff_sm


def fill_histograms(data, frames, hist_func, effs=None):
    """Fill the histograms"""
    hists = OrderedDict()
    weights = [None]
    if effs is not None:
        weights.append(effs(data))

    for frame in frames:
        for weight in weights:
            name = frame + ('_eff' if weight is not None else '')
            logging.debug('Filling histogram for {}'.format(name))
            hists[name] = hist_func(data, frame, weights=weight)

    return hists


def get_costh_phi_in_bins(hist_3d):
    """Get all costh phi histograms in each bin of the 3rd variable"""
    arr = get_array(hist_3d)
    binning = np.array([get_binning(hist_3d, 'X'), get_binning(hist_3d, 'Y')])
    err = get_array(hist_3d, errors=True)
    return [from_array(arr[:,:,i], binning, errors=err[:,:,i])
            for i in xrange(arr.shape[2])]


def store_hists(outfile, hists, basename, binvar=None):
    """Store histograms"""
    outfile.cd()
    if binvar is not None:
        # store 2d projections onto costh-phi
        for name, hist in hists.iteritems():
            projections = get_costh_phi_in_bins(hist)
            var_binning = get_binning(hist, 'Z')
            for ibin, proj in enumerate(projections):
                bin_bord = '{:.2f}_{:.2f}'.format(var_binning[ibin],
                                                  var_binning[ibin + 1])
                bin_bord = bin_bord.replace('.', 'p').replace('-', 'm')

                proj.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
                proj.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
                proj.SetName('_'.join(['proj', basename, binvar, bin_bord, name]))
                proj.Write()

            # also store the 3d maps that are used in the lookup
            hist.SetName('_'.join([basename, binvar, name]))
            hist.Write()

    else:
        for name, hist in hists.iteritems():
            hist.SetName('_'.join([basename, name]))
            hist.Write()


def create_histograms(data, frames, n_costh, n_phi, binvar, binning, effs=None):
    """
    Create histograms for each bin in the bin var
    """
    if binvar is not None:
        costh_binning = np.linspace(-1, 1, n_costh + 1)
        phi_binning = np.linspace(-180, 180, n_phi + 1)
        hist_sett = (n_costh, costh_binning, n_phi, phi_binning,
                     len(binning) - 1, binning)
        def fill_3d(data, frame, weights):
            """3D histograms"""
            return hist3d(_get_var(data, 'costh_' + frame),
                          _get_var(data, 'phi_' + frame),
                          _get_var(data, binvar),
                          hist_sett=hist_sett,
                          weights=weights)
        return fill_histograms(data, frames, fill_3d, effs)

    else:
        def fill_2d(data, frame, weights):
            """2D histograms"""
            return hist2d(_get_var(data, 'costh_' + frame),
                          _get_var(data, 'phi_' + frame),
                          hist_sett=(n_costh, -1, 1, n_phi, -180, 180),
                          weights=weights)

        return fill_histograms(data, frames, fill_2d, effs)


def main(args):
    """Main"""
    frames = args.frames.split(',')
    if args.bin_variable is not None:
        bin_var = args.bin_variable
        if args.binning is None:
            logging.error('You have to define a binning for \'{}\' if you want '
                          'to use it as binning variable'.format(bin_var))
            sys.exit(1)

        binning = parse_binning(args.binning)
    else:
        binning = None
        bin_var = None

    logging.info('Processing gen level file')
    gen_hists = create_histograms(get_dataframe(args.genlevelfile), frames,
                                  args.ncosth, args.nphi, bin_var, binning)

    logging.info('Processing reco level file')
    reco_hists = create_histograms(get_dataframe(args.recolevelfile), frames,
                                   args.ncosth, args.nphi, bin_var, binning,
                                   WEIGHT_F)

    logging.debug('Scaling gen level hists by '.format(args.scale_gen))
    for hist in gen_hists.values():
        hist.Scale(args.scale_gen)

    logging.info('calculating acceptance maps')
    acc_maps = OrderedDict()
    for name, hist in reco_hists.iteritems():
        gen_map_n = [n for n in gen_hists if n in name]
        if len(gen_map_n) > 1:
            logging.warning('Found more than 1 gen level map for \'{}\': {}'
                            .format(name, gen_map_n))

        # Still use just the first, since we should always just have 1
        acc_maps[name] = divide(hist, gen_hists[gen_map_n[0]])


    logging.debug('storing to output file')
    outfile = r.TFile(args.outfile, 'recreate' if args.recreate else 'update')
    store_hists(outfile, gen_hists, 'gen_costh_phi', bin_var)
    store_hists(outfile, reco_hists, 'reco_costh_phi', bin_var)
    store_hists(outfile, acc_maps, 'acc_map_costh_phi', bin_var)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to create and store '
                                     'acceptance correction maps')
    parser.add_argument('genlevelfile', help='File containing the gen level '
                        'costh,phi values')
    parser.add_argument('recolevelfile', help='File containing the reco level '
                        'costh,phi values')
    parser.add_argument('outfile', help='Output file into which the maps will be'
                        ' stored', default='corrmaps.root')
    parser.add_argument('-r', '--recreate', default=False, action='store_true',
                        help='Recreate the output file instead of updating it')
    parser.add_argument('-nc', '--ncosth', help='Number of bins in costh',
                        default=32)
    parser.add_argument('-np', '--nphi', help='Number of bins in phi',
                        default=45)
    parser.add_argument('-f', '--frames', help='comma separated list of frames '
                        'to produce correction maps in', default='PX,HX,CS')
    parser.add_argument('-s', '--scale-gen', help='Scale the gen level maps by '
                        'this number before doing the division', default=1.0,
                        type=float)
    parser.add_argument('-v', '--bin-variable', help='Add another binning in '
                        'this variable', default=None)
    parser.add_argument('-b', '--binning', help='Use this binning for the '
                        'bin variable', default=None)

    clargs = parser.parse_args()
    main(clargs)
