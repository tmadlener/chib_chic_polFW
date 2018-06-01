#!/usr/bin/env python
"""
Module providing the functionality to easily provide single muon and photon
efficiencies binned in eta as a function of pt
"""

import re

import numpy as np
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.WARNING,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.misc_helpers import stringify, flatten


def get_eta(eff_name):
    """
    Get the eta bin borders from the efficiency name
    """
    float_rgx = r'(-?\d\.?\d*)'
    eta_rgx = r'_'.join([r'eta', float_rgx, float_rgx])
    match = re.search(eta_rgx, stringify(eff_name, reverse=True))
    if match:
        return float(match.group(1)), float(match.group(2))


def get_effs(rfile, eff_type):
    """
    Get all the efficiencies matching the eff_type
    """
    keys = [k.GetName() for k in rfile.GetListOfKeys()]
    return {get_eta(k): rfile.Get(k) for k in keys if eff_type in k}


def get_eta_binning(effs):
    """
    Get the eta binning from the efficiencies and also check that the eta range
    is contiguous throughout the whole range
    """
    # first get all unique bin borders than check if all adjacent pairs are
    # also present in the efficiencies
    # to facilitate the task sort the unique borders
    uniq_bin_bord = sorted({b for b in flatten(effs)})
    possible_bins = zip(uniq_bin_bord[:-1], uniq_bin_bord[1:])

    for pbin in possible_bins:
        if pbin not in effs:
            logging.error('The eta binning is not contiguous! The possible bin '
                          '{:.1f} - {:.1f} is not present'.format(*pbin))
            return None

    return uniq_bin_bord


class EfficiencyProvider(object):
    """
    Class providing interface to easily obtain single muon and photon
    efficiencies
    """
    def __init__(self, eff_file, eff_type):
        """
        Args:
            eff_file (str): name of the ROOT file where the efficiencies are
                stored
            eff_type (str): either 'photon' or 'muon' depending on which
                efficiencies are desired
        """
        self.eff_file = r.TFile.Open(eff_file)
        self.effs = get_effs(self.eff_file, eff_type)
        self.eta_binning = get_eta_binning(self.effs)


    def eval(self, pt, eta):
        """
        Evaluate the efficiency at a given pt and eta value
        """
        eta_bin = self._find_eta_bin(np.abs(eta))
        if eta_bin is not None:
            return self.effs[eta_bin].Eval(pt)

        return 0


    def _find_eta_bin(self, eta):
        """
        Find the eta bin
        """
        for i in xrange(len(self.eta_binning) - 1):
            if eta >= self.eta_binning[i] and eta < self.eta_binning[i + 1]:
                return self.eta_binning[i], self.eta_binning[i + 1]

        logging.warning('Could not find efficiency for eta value of {:.4f}'
                        .format(eta))

        return None


class MuonEfficiencies(EfficiencyProvider):
    """
    Single muon efficiencies binned in eta as function of pT
    """
    def __init__(self, eff_file):
        """
        Args:
            eff_file (str): name of the ROOT file where the efficiencies are
                stored
        """
        EfficiencyProvider.__init__(self, eff_file, 'muon')


class PhotonEfficiencies(EfficiencyProvider):
    """
    Photon efficiencies binned in eta as function of pT
    """
    def __init__(self, eff_file):
        """
        Args:
            eff_file (str): name of the ROOT file where the efficiencies are
                stored
        """
        EfficiencyProvider.__init__(self, eff_file, 'photon')
