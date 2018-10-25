#!/usr/bin/env python
"""
Script that simply stores all histograms in a root file
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

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

    can = mkplot(hist, drawOpt=draw_opt)

    name ='/'.join([out_dir, hist.GetName()])
    for form in formats:
        can.SaveAs('.'.join([name, form]))


def main(args):
    """Main"""
    histfile = r.TFile.Open(args.inputfile)
    cond_mkdir(args.outdir)

    formats = args.extensions.split(',')
    for hist in list_obj(histfile, 'TH1', args.filter):
        store_hist(histfile.Get(hist), args.outdir, formats, args.draw_opt,
                   xrange=args.xrange, yrange=args.yrange, zrange=args.zrange)


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


    clargs = parser.parse_args()
    main(clargs)
