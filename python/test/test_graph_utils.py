#!/usr/bin/env python
"""
Test the graph utils module
"""

import unittest

import numpy as np
import numpy.testing as npt

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.misc_helpers import create_random_str

import utils.graph_utils as gu

class TestGraphRange(unittest.TestCase):
    def test_tgraph(self):
        """Test that TGraph are handled correctly"""
        n_vals = 25
        x_vals = np.random.uniform(0, 10, n_vals)
        y_vals = np.random.uniform(0, 10, n_vals)


        in_graph = r.TGraph(n_vals, x_vals, y_vals)

        r_low, r_high = 3, 8
        in_range_idcs = (x_vals > r_low) & (x_vals < r_high)
        r_graph = gu.graph_in_range(in_graph, r_low, r_high)

        self.assertTrue(isinstance(r_graph, r.TGraph))
        npt.assert_allclose(np.array(r_graph.GetX()),
                            x_vals[in_range_idcs])
        npt.assert_allclose(np.array(r_graph.GetY()),
                            y_vals[in_range_idcs])


    def test_tgraph_errors(self):
        """Test that TGraphErrors are handled correctly"""
        n_vals = 25
        x_vals = np.random.uniform(0, 10, n_vals)
        y_vals = np.random.uniform(0, 10, n_vals)
        x_errs = np.random.uniform(0, 1, n_vals)
        y_errs = np.random.uniform(0, 1, n_vals)

        in_graph = r.TGraphErrors(n_vals, x_vals, y_vals, x_errs, y_errs)

        r_low, r_high = 3, 8
        in_range_idcs = (x_vals > r_low) & (x_vals < r_high)
        r_graph = gu.graph_in_range(in_graph, r_low, r_high)

        self.assertTrue(isinstance(r_graph, r.TGraphErrors))

        npt.assert_allclose(np.array(r_graph.GetX()),
                            x_vals[in_range_idcs])
        npt.assert_allclose(np.array(r_graph.GetY()),
                            y_vals[in_range_idcs])


        rx_errs, ry_errs = gu.get_errors(r_graph)
        npt.assert_allclose(rx_errs, x_errs[in_range_idcs])
        npt.assert_allclose(ry_errs, y_errs[in_range_idcs])


    def test_tgraph_asymm_errors(self):
        """Test that TGraphAsymmErrors are handled correctly"""
        n_vals = 25
        x_vals = np.random.uniform(0, 10, n_vals)
        y_vals = np.random.uniform(0, 10, n_vals)
        xl_errs = np.random.uniform(0, 1, n_vals)
        yl_errs = np.random.uniform(0, 1, n_vals)
        xh_errs = np.random.uniform(0, 1, n_vals)
        yh_errs = np.random.uniform(0, 1, n_vals)

        in_graph = r.TGraphAsymmErrors(n_vals, x_vals, y_vals,
                                       xl_errs, xh_errs, yl_errs, yh_errs)

        r_low, r_high = 3, 8
        in_range_idcs = (x_vals > r_low) & (x_vals < r_high)
        r_graph = gu.graph_in_range(in_graph, r_low, r_high)

        self.assertTrue(isinstance(r_graph, r.TGraphAsymmErrors))

        npt.assert_allclose(np.array(r_graph.GetX()),
                            x_vals[in_range_idcs])
        npt.assert_allclose(np.array(r_graph.GetY()),
                            y_vals[in_range_idcs])

        rlx_errs, rhx_errs, rly_errs, rhy_errs = gu.get_errors(r_graph)
        npt.assert_allclose(rlx_errs, xl_errs[in_range_idcs])
        npt.assert_allclose(rly_errs, yl_errs[in_range_idcs])
        npt.assert_allclose(rhx_errs, xh_errs[in_range_idcs])
        npt.assert_allclose(rhy_errs, yh_errs[in_range_idcs])


if __name__ == '__main__':
    unittest.main()
