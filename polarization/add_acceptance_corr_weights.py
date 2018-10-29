#!/usr/bin/env python
"""
Script to add state dependent acceptance correction weights to the passed events
"""
import sys

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.EfficiencyProvider import AcceptanceCorrectionProvider
from utils.data_handling import (
    get_dataframe, add_branch, check_branch_available, get_treename
)
from utils.misc_helpers import _get_var

def branch_already_present(infile, treename, branchname):
    """Check if a branch is already present"""
    rfile = r.TFile.Open(infile)
    tree = rfile.Get(treename)
    return check_branch_available(tree, branchname, nowarn=True)


def main(args):
    """Main"""
    treename = get_treename(args.inputfile)
    branchname = args.branchname
    while branch_already_present(args.inputfile, treename, branchname):
        logging.warning('\'{}\' is already present in tree \'{}\' in file {}'
                        .format(branchname, treename, args.inputfile))
        branchname = raw_input('Please enter a different branch name: ')

    map_file = r.TFile.Open(args.corrmapfile)
    corr_map = map_file.Get(args.name)
    if not corr_map:
        logging.fatal('Cannot find acceptance map \'{}\' in file {}'
                      .format(args.name, args.corrmapfile))
        sys.exit(1)

    acc_prov = AcceptanceCorrectionProvider(corr_map)
    variables = ['{costh,phi}_PX']
    if acc_prov.dim == 3:
        variables.append(args.variable)

    data = get_dataframe(args.inputfile, treename, columns=variables)

    if acc_prov.dim == 2:
        corr_weights = acc_prov.eval(data.costh_PX, data.phi_PX)
    else:
        corr_weights = acc_prov.eval(data.costh_PX, data.phi_PX,
                                     _get_var(data, args.variable))
    add_branch(corr_weights, branchname, args.inputfile, treename)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to add state dependent '
                                     'acceptance correction weights')
    parser.add_argument('inputfile', help='File to process')
    parser.add_argument('corrmapfile', help='File with corrmap in it')
    parser.add_argument('-n', '--name', help='Name of the correction map in the '
                        'corrmap file', default='acc_map_costh_phi_JpsiPt_PX')
    parser.add_argument('-v', '--variable', help='Additional variable (apart '
                        'from costh and phi) in which the correction map is '
                        'binned (ignored if corrections are only in 2d)',
                        default='JpsiPt')
    parser.add_argument('-b', '--branchname', help='Name of the branch that is '
                        'created to store the correction weights in',
                        default='w_acc_corr')

    clargs = parser.parse_args()
    main(clargs)
