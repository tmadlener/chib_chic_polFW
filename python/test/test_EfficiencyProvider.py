#!/usr/bin/env python
"""
Tests for EfficiencyProvider classes
"""

import unittest
import numpy as np
import numpy.testing as npt

from utils.hist_utils import create_histogram, find_bin

from utils.EfficiencyProvider import AcceptanceCorrectionProvider


class TestAcceptanceCorrectionProvider(unittest.TestCase):
    def setUp(self):
        self.acc_map = create_histogram(np.random.uniform(0, 1, (10000, 2)),
                                        (20, 0, 1, 10, 0, 1))
        self.acc_map.Scale(1.0 / self.acc_map.Integral())


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


if __name__ == '__main__':
    unittest.main()
