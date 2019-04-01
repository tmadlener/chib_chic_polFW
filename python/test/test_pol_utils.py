#!/usr/bin/env python
"""
Tests for consistency between numpy and TF2 implementations of functions
"""

import unittest

import numpy as np
import numpy.testing as npt

import utils.pol_utils as pu


class Test2DAngularDistributions(unittest.TestCase):
    """Test for checking that the 2d distribution definitions are consistent"""
    def test_consistent(self):
        cth_tp = np.random.uniform(-1, 1, 1000)
        phi_tp = np.random.uniform(-180, 180, 1000)

        # test different sets of lambda just to be really sure
        for _ in xrange(5):
            lam = np.random.uniform(-1, 1, 3) # not necessarily physical lambdas
            tf2_lam = {'lth': lam[0], 'lph': lam[1], 'ltp': lam[2]}
            tf2 = pu.w_costh_phi(set_vals=tf2_lam) # no need to fix pars here
            tf2_vals = np.array([
                tf2.Eval(cth_tp[i], phi_tp[i]) for i in xrange(1000)
            ])
            np_vals = pu.ang_dist_2d(cth_tp, phi_tp, lam)

            # not using the default settings here to avoid false positives for
            # small return values
            npt.assert_allclose(np_vals, tf2_vals, rtol=0, atol=1e-7)
