#!/usr/bin/env python
"""
Tests for EfficiencyProvider classes
"""

import unittest
from mock import patch

import numpy as np
import numpy.testing as npt

from utils.hist_utils import create_histogram, find_bin

from utils.EfficiencyProvider import AcceptanceCorrectionProvider


class TestAcceptanceCorrectionProvider(unittest.TestCase):
    def setUp(self):
        self.acc_map = create_histogram(np.random.uniform(0, 1, (10000, 2)),
                                        (20, 0, 1, 10, 0, 1))
        self.acc_map.Scale(1.0 / self.acc_map.Integral())

        self.acc_map3 = create_histogram(np.random.uniform(0, 1, (100000, 3)),
                                         (10, 0, 1, 10, 0, 1, 5, 0, 1))
        self.acc_map3.Scale(1.0 / self.acc_map3.Integral())



    def test_masking_init(self):
        """
        Test if the masking in init works
        (i.e. if bins with 0 entries have -1 in it afterwards)
        """
        # Randomly set some bins to 0 in the acc map
        mask_bin_vals = np.random.uniform(0, 1, (10, 2))
        idcs = []
        for v in mask_bin_vals:
            idcs.append((
                self.acc_map.GetXaxis().FindBin(v[0]),
                self.acc_map.GetYaxis().FindBin(v[1])
            ))
            self.acc_map.SetBinContent(idcs[-1][0], idcs[-1][1], 0)

        corr_prov = AcceptanceCorrectionProvider(self.acc_map)
        self.assertEqual(np.sum(corr_prov.corr_map == -1), len(idcs))

        # When obtaining the masked values keep in mind that ROOT histograms start at 1
        masked_vals = np.array([corr_prov.corr_map[x-1, y-1] for (x, y) in idcs])
        npt.assert_equal(masked_vals, -1*np.ones(len(idcs)))


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


    @patch('utils.EfficiencyProvider.logging')
    def test_eval_warnings(self, mock_logger):
        corr_map2d = AcceptanceCorrectionProvider(self.acc_map)
        test_ps = np.random.uniform(0, 1, (10, 3))
        corr_map2d.eval(test_ps[:, 0], test_ps[:, 1], test_ps[:, 2])
        mock_logger.warning.assert_called_with('Ignoring third variable for 2d corrections')

        corr_map2d.eval(test_ps[:, 0], test_ps[:-1, 1])
        mock_logger.error.assert_called_with('Need same number of costh and phi values')

        corr_map3d = AcceptanceCorrectionProvider(self.acc_map3)
        corr_map3d.eval(test_ps[:, 0], test_ps[:, 1])
        mock_logger.error.assert_called_with('Need third variable for 3d corrections')

        corr_map3d.eval(test_ps[:, 0], test_ps[:, 1], test_ps[:-1, 2])
        mock_logger.error.assert_called_with('Need the same number of var and costh and phi values')


if __name__ == '__main__':
    unittest.main()
