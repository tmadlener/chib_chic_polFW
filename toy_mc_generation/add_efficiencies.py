#!/usr/bin/env python
"""
Script to add the single muon and photon efficiencies to the passed data
"""

import sys

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.data_handling import get_dataframe, store_dataframe, apply_selections
from utils.EfficiencyProvider import MuonEfficiencies, PhotonEfficiencies

# names of the branches
MUON_NAMES = ['lepP_sm', 'lepN_sm']
PHOTON_NAME = 'gamma_sm'


def get_effs(eff_file, oride_eff, eff_proto):
    """
    Setup the efficiencies using either a the orid_eff file (if set) or the
    efficiency file.
    Returns an instance of the eff_proto
    """
    used_file = oride_eff if oride_eff is not None else eff_file
    if used_file is None:
        logging.error('Need an efficiency file to setup {}'.format(eff_proto))
        sys.exit(1)

    return eff_proto(used_file)


def calc_effs(data, effs, branch):
    """
    Calc the efficiencies for a given branch (i.e. particle) and return the
    results
    """
    pt_branch = 'pT_' + branch
    eta_branch = 'eta_' + branch

    logging.info('Calculating efficiencies for {}'.format(branch))
    return map(effs.eval, data.loc[:, pt_branch], data.loc[:, eta_branch])


def get_treename(filename):
    """
    Get the name of the ONLY tree in the rootfile

    NOTE: this does not do any checks at the moment!
    """
    rfile = r.TFile.Open(filename)
    keys = [k.GetName() for k in rfile.GetListOfKeys()]
    return keys[0]


def main(args):
    """Main"""
    muon_effs = get_effs(args.efficiencies, args.muoneffs, MuonEfficiencies)
    phot_effs = get_effs(args.efficiencies, args.photoneffs, PhotonEfficiencies)
    treename = get_treename(args.datafile) # do this here and fail early in case
    data = get_dataframe(args.datafile)

    # TODO:
    # Decide if we want selections here already, currently do not use any
    # selection. This could mean that some events are assigned efficiencies
    # that might be outside of the valid range of where they have been
    # determined
    # TODO: select only those events for which we have efficiencies and maybe
    # make that cl args
    data.loc[:, PHOTON_NAME + 'Eff'] = calc_effs(data, phot_effs, PHOTON_NAME)
    for muon in MUON_NAMES:
        data.loc[:, muon + 'Eff'] = calc_effs(data, muon_effs, muon)

    outfile = args.outfile if args.outfile else args.datafile
    store_dataframe(data, outfile, treename)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to add the single muon '
                                     'and photon efficiences')
    parser.add_argument('datafile', help='file containing all the data to which '
                        'the efficiencies should be added')
    parser.add_argument('-e', '--efficiencies', help='file containing photon and'
                        ' single muon efficiencies', default=None)
    parser.add_argument('-em', '--muoneffs', help='file containing single muon '
                        'efficiencies. overrides efficiencies argument',
                        default=None)
    parser.add_argument('-ep', '--photoneffs', help='file containing photon '
                        'efficiencies. Overrides efficiencies argument',
                        default=None)
    parser.add_argument('-o', '--outfile', help='Store output to a new file '
                        'instead of updating the input file', type=str,
                        default='')

    clargs = parser.parse_args()
    main(clargs)
