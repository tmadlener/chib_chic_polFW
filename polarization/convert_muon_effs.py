#!/usr/bin/env python
"""
Convert the muon efficiencies from another root file and store them using the
same conventions as the photon efficiencies
"""

import re

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.graph_utils import get_binning
from utils.misc_helpers import find_common_binning, get_bin, get_bin_edges

from common_func import get_name, get_eta_range

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
    if bin_idx is None: # if AETA is not in the name, try again with etaBin
        bin_idx = get_bin_idx(eff_name, 'etaBin')
    return ETA_RANGE[bin_idx], ETA_RANGE[bin_idx + 1]


def get_effs(eff_file, eff_base):
    """
    Get all the efficiencies matching the eff_base
    """
    keys = [k.GetName() for k in eff_file.GetListOfKeys()]
    effs = [eff_file.Get(k) for k in keys if eff_base in k]
    return effs


def make_2D_plot(efficiency_graphs):
    """
    Make a 2D plot from all efficiency graphs (assuming the new naming scheme)
    """
    bin_borders = [get_binning(g) for g in efficiency_graphs]
    eff_binning = find_common_binning(bin_borders)
    eff_binning = get_bin_edges(eff_binning)

    if eff_binning is None:
        logging.info('Cannot create a 2D efficiency plot because no common '
                     'binning could be found')
        return None

    eta_binning = np.array(ETA_RANGE)
    hist_eff = r.TH2D("efficiencies_2D", ";p_{T} / GeV;|#eta|",
                      len(eff_binning) - 1, eff_binning,
                      len(eta_binning) - 1, eta_binning)

    for graph in efficiency_graphs:
        eta_bin = get_eta_range(graph.GetName())
        # use the average of eta, since that is certainly in the bin
        ieta = get_bin(eta_binning, 0.5 * (eta_bin[0] + eta_bin[1]))

        effs = np.array(graph.GetY())

        for ieff, val in enumerate(np.array(graph.GetX())):
            ibin = get_bin(eff_binning, val)
            if ibin >= 0:
                hist_eff.SetBinContent(ibin + 1, ieta + 1, effs[ieff])


    return hist_eff


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

    if args.plot:
        hist_eff = make_2D_plot(old_effs)
        if hist_eff is not None:
            hist_eff.Write()

    outfile.Write('', r.TObject.kWriteDelete)
    outfile.Close()


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
    parser.add_argument('-p', '--plot', help='make a 2D plot of the resulting '
                        'efficiencies and also store that to the produced file',
                        action='store_true', default=False)

    clargs = parser.parse_args()
    main(clargs)
