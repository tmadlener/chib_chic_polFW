#!/usr/bin/env python
"""
Convert the muon efficiencies from another root file and store them using the
same conventions as the photon efficiencies
"""

import re

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from common_func import get_name

# The abs(eta) bin borders in which the single muon effs have been determined
ETA_RANGE = (0., 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

def get_bin_idx(eff_name, bin_var):
    """
    Get the bin index from a given efficiency name and the variable that has
    been used for binning
    """
    bin_rgx = r''.join([bin_var, r'([0-9]+)'])
    match = re.search(bin_rgx, eff_name)
    if match:
        return int(match.group(1))


def get_eta_bin(eff_name):
    """
    Get the eta bin borders from the efficiency name
    """
    bin_idx = get_bin_idx(eff_name, 'AETA')
    return ETA_RANGE[bin_idx], ETA_RANGE[bin_idx + 1]


def get_effs(eff_file, eff_base):
    """
    Get all the efficiencies matching the eff_base
    """
    keys = [k.GetName() for k in eff_file.GetListOfKeys()]
    effs = [eff_file.Get(k) for k in keys if eff_base in k]
    return effs


def main(args):
    """Main"""
    infile = r.TFile.Open(args.infile)
    old_effs = get_effs(infile, args.efficiencies)

    file_option = 'update' if args.update else 'recreate'
    outfile = r.TFile(args.outfile, file_option)
    outfile.cd()
    for eff in old_effs:
        eta_bin = get_eta_bin(eff.GetName())
        eff.SetName(get_name(eta_bin, 'muon_eff_pt'))
        eff.Write()

    outfile.cd()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for converting single '
                                     'efficiencies from the old naming scheme '
                                     'to the one used in this framework')
    parser.add_argument('infile', help='Input file containing the efficiencies '
                        'with the old naming scheme')
    parser.add_argument('-o', '--outfile', help='Output file to which the '
                        'efficiencies will be stored with the new naming '
                        'scheme', default='single_muon_effs.root')
    parser.add_argument('-u', '--update', help='update the output file instead '
                        'of recreating it', default=False, action='store_true')
    parser.add_argument('-e', '--efficiencies', help='Base name of the '
                        'efficiencies to convert', default='gEff_MC_PT_')

    clargs = parser.parse_args()
    main(clargs)
