#!/usr/bin/env python
"""
Script to determine correction maps
"""

from collections import OrderedDict

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True


from utils.data_handling import get_dataframe
from utils.hist_utils import divide, hist2d
from utils.misc_helpers import _get_var

WEIGHT_F = lambda d: d.gamma_eff_sm * 0.01 * d.lepP_eff_sm * d.lepN_eff_sm

def create_histograms(data, frames, n_costh, n_phi, effs=None):
    """
    Create all histograms
    """
    hists = OrderedDict()
    weights = [None]
    if effs is not None:
        weights.append(effs(data))

    for frame in frames:
        for weight in weights:
            name = frame + ('_eff' if weight is not None else '')
            hists[name] = hist2d(_get_var(data, 'costh_' + frame),
                                 _get_var(data, 'phi_' + frame),
                                 hist_sett=(n_costh, -1, 1, n_phi, -180, 180),
                                 weights=weight)

    return hists


def store_hists(outfile, hists, basename):
    """Store histograms"""
    outfile.cd()
    for name, hist in hists.iteritems():
        hist.SetName('_'.join([basename, name]))
        hist.Write()


def main(args):
    """Main"""
    frames = args.frames.split(',')
    gen_hists = create_histograms(get_dataframe(args.genlevelfile), frames,
                                  args.ncosth, args.nphi)

    reco_hists = create_histograms(get_dataframe(args.recolevelfile), frames,
                                   args.ncosth, args.nphi, WEIGHT_F)

    for hist in gen_hists.values():
        hist.Scale(args.scale_gen)

    acc_maps = OrderedDict()
    for name, hist in reco_hists.iteritems():
        acc_maps[name] = divide(hist, gen_hists[name.split('_')[0]])

    outfile = r.TFile(args.outfile, 'recreate' if args.recreate else 'update')
    store_hists(outfile, gen_hists, 'gen_costh_phi')
    store_hists(outfile, reco_hists, 'reco_costh_phi')
    store_hists(outfile, acc_maps, 'acc_map_costh_phi')


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

    clargs = parser.parse_args()
    main(clargs)
