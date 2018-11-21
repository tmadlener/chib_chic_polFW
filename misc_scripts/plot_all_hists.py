#!/usr/bin/env python
"""
Script that simply stores all histograms in a root file
"""

import sys

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.ERROR,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.plot_helpers import mkplot
from utils.misc_helpers import cond_mkdir
from utils.data_handling import list_obj

def store_hist(hist, out_dir, formats, draw_opt, **kwargs):
    """Store one hist"""
    for axis in 'xyz':
        ax_range = kwargs.pop(axis + 'range', None)
        if ax_range is not None:
            range_vals = [float(v) for v in ax_range.split(',')]
            haxis = getattr(hist, 'Get' + axis.upper() + 'axis')()
            haxis.SetRangeUser(*range_vals)


    hist_cmd = kwargs.pop('hist_cmd', None)
    if hist_cmd is not None:
        if hasattr(hist, hist_cmd):
            hist = getattr(hist, hist_cmd)()
            if not isinstance(hist, (r.TH1, r.TGraph, r.TF1)):
                logging.fatal('Result of \'{}\' is not a plotable object'
                              .format(hist_cmd))
                sys.exit(1)
        else:
            logging.error('Cannot call \'{}\' on histogram object.'
                          .format(hist_cmd))

    hist.SetStats(0)
    can = mkplot(hist, drawOpt=draw_opt)

    name = '/'.join([out_dir, hist.GetName()])
    for form in formats:
        can.SaveAs('.'.join([name, form]))


def main(args):
    """Main"""
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    else:
        r.gROOT.ProcessLine('gErrorIgnoreLevel = 1001')

    histfile = r.TFile.Open(args.inputfile)
    cond_mkdir(args.outdir)

    formats = args.extensions.split(',')
    for hist in list_obj(histfile, 'TH1', args.filter):
        store_hist(histfile.Get(hist), args.outdir, formats, args.draw_opt,
                   xrange=args.xrange, yrange=args.yrange, zrange=args.zrange,
                   hist_cmd=args.hist_cmd)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to simply plot and '
                                     'store all histograms of a root file')
    parser.add_argument('inputfile', help='File with the histograms')
    parser.add_argument('-o', '--outdir', help='Directory to which the files '
                        'aree stored', default='./')
    parser.add_argument('-x', '--extensions', help='comma separated list of '
                        'image formats that should be created', default='pdf')
    parser.add_argument('-d', '--draw-opt', help='Use a specific draw option',
                        default='')
    parser.add_argument('-f', '--filter', help='Only plot the histograms for '
                        'the name matches this string', default=None)
    parser.add_argument('-zr', '--zrange', help='Constrain the z-range of a plot'
                        ' (mostly for the color axis in colz plots)',
                        default=None)
    parser.add_argument('-yr', '--yrange', help='Constrain the y-range of a '
                        'plot', default=None)
    parser.add_argument('-xr', '--xrange', help='Constrain the x-range of a '
                        'plot', default=None)
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Enable printouts to screen')
    parser.add_argument('-c', '--hist-cmd', help='Call this function (input is '
                        'the name as string) on each histogram and plot the '
                        'result', default=None)
    clargs = parser.parse_args()
    main(clargs)
