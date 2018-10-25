#!/usr/bin/env python
"""
Script that divides all histograms from one file by the histograms with the same
name from another file
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.WARNING,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.hist_utils import divide
from utils.data_handling import list_obj


def main(args):
    """Main"""
    numfile = r.TFile.Open(args.numfile)
    denomfile = r.TFile.Open(args.denomfile)

    # use sets so that it is easier to find the histograms with the same names
    num_hists = set(list_obj(numfile, 'TH1', args.filter))
    denom_hists = set(list_obj(denomfile, 'TH1', args.filter))
    ratio_hists = num_hists.intersection(denom_hists)

    outfile = r.TFile(args.outfile, 'update' if args.update else 'recreate')
    for ratio_n in ratio_hists:
        logging.debug('Creating: {}'.format(ratio_n))
        ratio = divide(numfile.Get(ratio_n), denomfile.Get(ratio_n))
        outfile.cd()
        ratio.Write()

    outfile.Write('', r.TObject.kWriteDelete)
    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script that divides hists from'
                                     ' one file by hists from another file if '
                                     'they have the same name')
    parser.add_argument('numfile', help='File containing the numerator '
                        'histograms')
    parser.add_argument('denomfile', help='File containing the denominator '
                        'histograms')
    parser.add_argument('outfile', help='Output file that contains the ratio '
                        'histograms', default='ratios.root')
    parser.add_argument('-u', '--update', default=False, action='store_true',
                        help='Update the output file instead of recreating it')
    parser.add_argument('-f', '--filter', help='Only plot the histograms for '
                        'the name matches this string', default=None)


    clargs = parser.parse_args()
    main(clargs)
