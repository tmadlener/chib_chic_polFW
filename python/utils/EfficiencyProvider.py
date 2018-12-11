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
    corrections depending on 2D (costh-phi) variables or 3D (costh-phi plus
    another variable)
    """
    def __init__(self, acc_map, min_acc=0, mask_prec=None, mask=None):
        """
        Args:
            acc_map (TH2D, TH3D or THnD): costh-phi map or costh-phi-var map
                obtained by applying all cuts and selections (and possibly
                efficiency weightings). For each bin 1 / (bin content) will be
                the weight for the acceptance correction
            min_acc (float, optional): Mask all bins with an acceptance below
                this value (default = 0)
            mask_prec (float, optional): If not None, mask all bins for which
                the relative error is larger than the passed value
            mask (np.array, optional): Array with the same dimensions as the
                acceptance map. All bins containing a non False value will be
                masked. Overrides the min_acc and mask_prec argument (i.e. they
                will be ignored) but still respects zero bin masking
        """
        self.hist = acc_map
        logging.debug('Using acceptance map \'{}\''.format(self.hist.GetName()))
        # Corrections are 1 / acceptance map
        acc_values = get_array(self.hist)

        if mask is not None:
            if mask.shape != acc_values.shape:
                logging.error('mask and acceptance map need to have the same '
                              'dimensions. mask: {}, map: {}'
                              .format(mask.shape, acc_values.shape))
            if min_acc != 0:
                logging.info('Ignoring min_acc={} because a mask is used'
                             .format(min_acc))
            if mask_prec is not None:
                logging.info('Ignoring mask_prec={} because a mask is used'
                             .format(mask_prec))
            # mask the values without acceptance in the acceptance map
            # this will also make them return -1 for the correction map
            logging.debug('Masking {} bins according to the mask'
                          .format(np.sum(mask)))
            empty_mask = (acc_values == 0)
            logging.debug('Masking {} empty bins'.format(np.sum(empty_mask)))
            masked_vals = empty_mask | mask
        else:
            if min_acc < 0 or min_acc > 1:
                logging.warning('The minimum acceptance should be a value '
                                'between 0 and 1, but is {}'.format(min_acc))
            masked_vals = (acc_values <= min_acc).astype(bool)
            logging.debug('Minimum acceptance = {}: Masking {} bins'
                          .format(min_acc, np.sum(masked_vals)))
            if mask_prec is not None:
                if isinstance(mask_prec, float):
                    acc_errs = get_array(self.hist, errors=True)
                    rel_uncer = np.divide(acc_errs, acc_values,
                                          where=acc_values!=0)
                    mask_uncer = (rel_uncer > mask_prec).astype(bool)
                    logging.debug('Minimum precision = {}: Masking {} bins'
                                  .format(mask_prec, np.sum(mask_uncer)))
                    masked_vals |= mask_uncer
                else:
                    logging.error('mask_prec has to be a float value. Not using'
                                  ' it to mask bins with too low precision.')

        acc_values = ~masked_vals * acc_values + -1 * masked_vals
        logging.debug('{} of {} bins are masked in the correction map'
                      .format(np.sum(masked_vals), acc_values.size))

        self.corr_map = 1.0 / acc_values
        self.var_binnings = []
        self.ndim = self.corr_map.ndim
        for i in xrange(self.ndim):
            self.var_binnings.append(get_binning(acc_map, i))


    def eval(self, *args):
        """
        Evaluate the correction map at all the passed values

        Args:
            np.arrays: As many arrays as there are dimensions in the correction
                map. All arrays must be one dimensional and have the same length

        Returns:
            np.array: Array of the values in the correction map at the given
                input coordinates. For coordinates where the acceptance map
                contained 0 (i.e. infinite correction) or where any masking
                condition was true -1 is returned
        """
        if len(args) != self.ndim:
            logging.error('Correction map has dimension {} but trying to '
                          'evaluate it with {}'.format(self.ndim, len(args)))
            return None

        shape = args[0].shape
        for arg in args:
            if arg.shape != shape:
                logging.error('All input values need to have the same shape')
                return None

        var_bins = []
        for i in xrange(self.ndim):
            var_bins.append(find_bin(self.var_binnings[i], args[i]))

        return self.corr_map[var_bins]
