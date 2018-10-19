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
from utils.hist_utils import get_array, get_binning, find_bin


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

        # no valid efficiency: mark by returning -1
        return -1


    def _find_eta_bin(self, eta):
        """
        Find the eta bin
        """
        for i in xrange(len(self.eta_binning) - 1):
            if eta >= self.eta_binning[i] and eta < self.eta_binning[i + 1]:
                return self.eta_binning[i], self.eta_binning[i + 1]

        logging.debug('Could not find efficiency for eta value of {:.4f}'
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


class AcceptanceCorrectionProvider(object):
    """
    Class providing interface to easily obtain acceptance (and efficiency)
    corrections depending on 2D (costh-phi) variables
    """
    def __init__(self, acc_map):
        """
        Args:
            acc_map (TH2D): costh-phi map obtained by applying all cuts and
                selections (and possibly efficiency weightings). For each bin
                1 / (bin content) will be the weight for the acceptance
                correction
        """
        self.hist = acc_map
        # Corrections are 1 / acceptance map
        acc_values = get_array(self.hist)
        # mask the values without acceptance in the acceptance map
        # this will also make them return -1 for the correction map
        acc_values -= 1 * (acc_values == 0)
        self.corr_map = 1.0 / acc_values
        self.costh_binning = get_binning(acc_map, 'X')
        self.phi_binning = get_binning(acc_map, 'Y')


    def eval(self, costh, phi):
        """
        Evaluate the correction map at all given costh and phi values

        Args:
            costh, phi (np.array): Arrays (of equal length) containing all pairs
                of costh and phi values

        Return:
            np.array: Array of the values in the correction map at the given
                 costh and phi coordinates. For coordinates where the acceptance
                 map contained 0 (i.e. infinite correction) -1 is returned
        """
        costh_bins = find_bin(self.costh_binning, costh)
        phi_bins = find_bin(self.phi_binning, phi)
        return self.corr_map[costh_bins, phi_bins]
