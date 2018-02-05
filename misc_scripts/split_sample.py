#!/usr/bin/env python
"""
Split passed root file into subsamples randomly
"""

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from numpy.random import permutation
from numpy import arange

from utils.data_handling import get_dataframe, store_dataframe


def main(args):
    """Main"""
    orig_frame = get_dataframe(args.infile, args.treename)
    rand_indices = permutation(arange(orig_frame.shape[0]) % args.nsamples)

    base_output = args.infile.replace('.root', '__{}.root')

    for i in xrange(args.nsamples):
        idx_frame = orig_frame.iloc[rand_indices == i]
        out_name = base_output.format(i)
        store_dataframe(idx_frame, out_name, args.treename)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for splitting a passed'
                                     'root file into several subsamples')
    parser.add_argument('infile', help='input file that should be split')
    parser.add_argument('-n', '--nsamples', type=int, default=2,
                        help='the number of desired (mutually exclusive) '
                        'subsamples')
    parser.add_argument('-t', '--treename', type=str, default='chic_mc_tuple',
                        help='treename under which sample is stored in root '
                        'file')

    clargs = parser.parse_args()
    main(clargs)
