#!/usr/bin/env python
"""
Script to add the single muon and photon efficiencies to the passed data
"""

import sys

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.data_handling import get_dataframe, add_branch
from utils.EfficiencyProvider import MuonEfficiencies, PhotonEfficiencies

# names of the branches
MUON_NAMES = ['muP', 'muN']
PHOTON_NAME = 'photon'


def get_effs(eff_file, oride_eff, eff_proto):
    """
    Setup the efficiencies using either a the orid_eff file (if set) or the
    efficiency file.
    Returns an instance of the eff_proto
    """
    used_file = oride_eff if oride_eff is not None else eff_file
    if used_file is None:
        return None

    return eff_proto(used_file)


def calc_effs(data, effs, branch):
    """
    Calc the efficiencies for a given branch (i.e. particle) and return the
    results
    """
    pt_branch = branch + 'Pt'
    eta_branch = branch + 'Eta'

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

    # only read in the necessary columns and add newly created branches to file
    columns = [n + '{Pt,Eta}' for n in MUON_NAMES + [PHOTON_NAME]]
    data = get_dataframe(args.datafile, treename, columns=columns)

    if phot_effs is not None:
        photon_effs = np.array(calc_effs(data, phot_effs, PHOTON_NAME))
        add_branch(photon_effs, '_'.join(['gamma_sm', args.name]),
                   args.datafile, treename)

    if muon_effs is not None:
        for muon in MUON_NAMES:
            mu_effs = np.array(calc_effs(data, muon_effs, muon))
            lep = 'lepP' if 'P' in muon else 'lepN' # for consistency
            add_branch(mu_effs, '_'.join([lep, 'sm', args.name]),
                       args.datafile, treename)


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
    # parser.add_argument('-o', '--outfile', help='Store output to a new file '
    #                     'instead of updating the input file', type=str,
    #                     default='')
    parser.add_argument('-n', '--name', help='Name of the efficiencies that '
                        'will be used in the output branch to identify them. '
                        'The branch will be named (e.g.) gamma_sm_[name]. '
                        'NOTE: The default will overwrite possibly existing '
                        'efficiencies from the generation.', type=str,
                        default='eff')


    clargs = parser.parse_args()
    main(clargs)
