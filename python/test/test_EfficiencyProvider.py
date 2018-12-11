#!/usr/bin/env python
"""
Tests for EfficiencyProvider classes
"""

import unittest
from mock import patch

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np
import numpy.testing as npt

from utils.hist_utils import create_histogram, get_array

from utils.EfficiencyProvider import AcceptanceCorrectionProvider


class TestAcceptanceCorrectionProvider(unittest.TestCase):
    def setUp(self):
        self.acc_map = create_histogram(np.random.uniform(0, 1, (10000, 2)),
                                        (20, 0, 1, 10, 0, 1))
        self.acc_map.Scale(1.0 / self.acc_map.Integral())

        self.acc_map3 = create_histogram(np.random.uniform(0, 1, (100000, 3)),
                                         (10, 0, 1, 10, 0, 1, 5, 0, 1))
        self.acc_map3.Scale(1.0 / self.acc_map3.Integral())


    def test_masking_zero_init(self):
        """
        Test if the masking in init works
        (i.e. if bins with 0 entries have -1 in it afterwards)
        """
        # Randomly set some bins to 0 in the acc map
        mask_bin_vals = np.random.uniform(0, 1, (10, 2))
        idcs = set()
        for v in mask_bin_vals:
            idxx = self.acc_map.GetXaxis().FindBin(v[0])
            idxy = self.acc_map.GetYaxis().FindBin(v[1])
            idcs.add((idxx, idxy))
            self.acc_map.SetBinContent(idxx, idxy, 0)

        corr_prov = AcceptanceCorrectionProvider(self.acc_map)
        self.assertEqual(np.sum(corr_prov.corr_map == -1), len(idcs))

        # When obtaining the masked values keep in mind that ROOT histograms start at 1
        masked_vals = np.array([corr_prov.corr_map[x-1, y-1] for (x, y) in idcs])
        npt.assert_equal(masked_vals, -1*np.ones(len(idcs)))


    def test_masking_low_prec_init(self):
        """
        Test if the masking in init works for bins with too low precision
        """
        # sprt(1e4 / (20*10)) / (1e4 / 20 * 10), average poisson uncertainty
        min_prec = 0.143 # tuned to be slightly below the expected uncertainty
        arr, err = get_array(self.acc_map), get_array(self.acc_map, errors=True)
        masked_bins = (err / arr  > min_prec).astype(bool)

        corr_prov = AcceptanceCorrectionProvider(self.acc_map, mask_prec=min_prec)
        # the correction map has to have the same bins masked as the acceptance map
        npt.assert_equal((corr_prov.corr_map == -1).astype(bool), masked_bins)


    def test_masking_min_acc_init(self):
        """Test if masking works for bins with too low acceptance (i.e. content)"""
        min_acc = 0.00475 # (slightly below average occupancy)
        arr = get_array(self.acc_map)
        masked_bins = (arr < min_acc).astype(bool)

        corr_prov = AcceptanceCorrectionProvider(self.acc_map, min_acc=min_acc)
        npt.assert_equal((corr_prov.corr_map == -1).astype(bool), masked_bins)


    @patch('utils.EfficiencyProvider.logging') # to catch info messages
    def test_masking_mask_init(self, mock_logger):
        """Test if using a mask works as intended"""
        mask = np.random.uniform(0, 1, (20, 10)) < 0.5 # mask half randomly
        corr_prov = AcceptanceCorrectionProvider(self.acc_map, mask=mask)
        npt.assert_equal((corr_prov.corr_map == -1).astype(bool), mask)

        # other arguments are ignored if a mask is passed
        # NOTE: using the hardest possible cuts for them
        corr_prov = AcceptanceCorrectionProvider(self.acc_map, mask=mask,
                                                 min_acc=1, mask_prec=0)
        npt.assert_equal((corr_prov.corr_map == -1).astype(bool), mask)

        # zero bins are still masked
        idcs = set()
        mask_bin_vals = np.random.uniform(0, 1, (10, 2))
        for v in mask_bin_vals:
            idxx = self.acc_map.GetXaxis().FindBin(v[0])
            idxy = self.acc_map.GetYaxis().FindBin(v[1])
            idcs.add((idxx, idxy))
            self.acc_map.SetBinContent(idxx, idxy, 0)

        zero_masked_bins = np.zeros((20, 10), dtype=bool)
        for x, y in idcs:
            zero_masked_bins[x-1, y-1] = True # ROOT binning starts at 1

        corr_prov = AcceptanceCorrectionProvider(self.acc_map, mask=mask)
        npt.assert_equal((corr_prov.corr_map == -1).astype(bool), mask | zero_masked_bins)


    def test_masking_init(self):
        """
        Test if masking also works when finding zero bins and using precision
        limit or acceptance limit
        """
        # Mask the precision bins first, to avoid division by zero
        min_prec = 0.143
        arr, err = get_array(self.acc_map), get_array(self.acc_map, errors=True)
        masked_bins = (err / arr > min_prec).astype(bool)

        # mask the min acceptance bins
        min_acc = 0.00425
        masked_bins |= (arr < min_acc).astype(bool)

        mask_bin_vals = np.random.uniform(0, 1, (10, 2))
        idcs = set()
        for v in mask_bin_vals:
            idxx = self.acc_map.GetXaxis().FindBin(v[0])
            idxy = self.acc_map.GetYaxis().FindBin(v[1])
            idcs.add((idxx, idxy))
            self.acc_map.SetBinContent(idxx, idxy, 0)

        # add the bins that are masked due to being zero
        for x, y in idcs:
            masked_bins[x-1, y-1] = True # ROOT binning starts at 1

        corr_prov = AcceptanceCorrectionProvider(self.acc_map, mask_prec=min_prec)
        npt.assert_equal((corr_prov.corr_map == -1).astype(bool), masked_bins)


    def test_eval(self):
        corr_map_prov = AcceptanceCorrectionProvider(self.acc_map)
        test_points = np.random.uniform(0, 1, (10000, 2))
        corrs = corr_map_prov.eval(test_points[:, 0], test_points[:, 1])

        exp_corrs = np.array(
            [1.0 / self.acc_map.GetBinContent(self.acc_map.FindBin(v[0], v[1]))
             for v in test_points]
        )
        npt.assert_allclose(corrs, exp_corrs)


    def test_eval_3d(self):
        map_3d = AcceptanceCorrectionProvider(self.acc_map3)
        test_ps = np.random.uniform(0, 1, (10000, 3))
        corrs = map_3d.eval(test_ps[:, 0], test_ps[:, 1], test_ps[:, 2])

        exp_corrs = np.array(
            [1.0 / self.acc_map3.GetBinContent(self.acc_map3.FindBin(v[0], v[1], v[2]))
             for v in test_ps]
        )
        npt.assert_allclose(corrs, exp_corrs)


    def test_eval_4d(self):
        acc_map4 = r.THnD('', '', 4, np.array([10, 10, 5, 5], dtype='i4'),
                          np.zeros(4), np.ones(4))
        fill_vals = np.random.uniform(0, 1, (100000, 4))
        for i in xrange(fill_vals.shape[0]):
            acc_map4.Fill(fill_vals[i, :])
        acc_map4.Scale(1.0 / acc_map4.ComputeIntegral())

        map_4d = AcceptanceCorrectionProvider(acc_map4)
        test_ps = np.random.uniform(0, 1, (10000, 4))
        corrs = map_4d.eval(test_ps[:, 0], test_ps[:, 1], test_ps[:, 2], test_ps[:, 3])

        exp_corrs = np.array(
            [1.0 / acc_map4.GetBinContent(acc_map4.GetBin(v))
             for v in test_ps]
        )
        npt.assert_allclose(corrs, exp_corrs)



    @patch('utils.EfficiencyProvider.logging')
    def test_eval_warnings(self, mock_logger):
        corr_map2d = AcceptanceCorrectionProvider(self.acc_map)
        test_ps = np.random.uniform(0, 1, (10, 3))
        corr_map2d.eval(test_ps[:, 0], test_ps[:, 1], test_ps[:, 2])
        mock_logger.error.assert_called_with('Correction map has dimension 2 but'
                                             ' trying to evaluate it with 3')

        corr_map2d.eval(test_ps[:, 0], test_ps[:-1, 1])
        mock_logger.error.assert_called_with('All input values need to have the'
                                             ' same shape')


if __name__ == '__main__':
    unittest.main()
