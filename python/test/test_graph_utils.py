#!/usr/bin/env python
"""
Test the graph utils module
"""

import unittest

import numpy as np
import numpy.testing as npt

import sympy as sp

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.misc_helpers import create_random_str
from utils.symbolic import func_cov

import utils.graph_utils as gu

def get_random_graph(gtype, n_points=5):
    """Get a random graph"""
    x_vals = np.random.uniform(0, 10, n_points)
    y_vals = np.random.uniform(0, 10, n_points)

    try:
        errors = [np.random.uniform(0, 1, n_points) for _ in xrange(4)]
        return gtype(n_points, x_vals, y_vals, *errors)
    except TypeError:
        try:
            return gtype(n_points, x_vals, y_vals, *errors[:2])
        except TypeError:
            return gtype(n_points, x_vals, y_vals)

    return None


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


class TestDivideGraphs(unittest.TestCase):
    """Test that divide_graphs works as intended"""
    def __init__(self, *args, **kwargs):
        super(TestDivideGraphs, self).__init__(*args, **kwargs)
        x, y = sp.symbols('x, y')
        self.ratio_cov = func_cov(x / y, (x, y))[0]


    def _get_analytic_err(self, vx, vy, sx, sy, rho=0):
        """
        Calculate the uncertainty analytically
        """
        uncer = []
        # have to loop over the values manually and do the calculation for every
        # set of values
        for i, v in enumerate(vx):
            var = self.ratio_cov.subs({
                'x': v, 'y': vy[i],
                'sigma_x': sx[i], 'sigma_y': sy[i], 'rho_xy': rho
            })
            uncer.append(float(sp.sqrt(var)))

        return np.array(uncer)


    def test_tgraph(self):
        """Check that TGraphs are handled correctly"""
        n_graph = get_random_graph(r.TGraph, 5)
        d_graph = get_random_graph(r.TGraph, 5)

        r_graph = gu.divide_graphs(n_graph, d_graph)
        npt.assert_allclose(np.array(r_graph.GetY()),
                            np.array(n_graph.GetY()) / np.array(d_graph.GetY()))

        npt.assert_allclose(np.array(r_graph.GetX()), np.array(n_graph.GetX()))


    def test_tgraph_errors(self):
        """Check that the TGraphErrors are handled correctly"""
        n_graph = get_random_graph(r.TGraphErrors, 5)
        d_graph = get_random_graph(r.TGraphErrors, 5)
        r_graph = gu.divide_graphs(n_graph, d_graph)

        npt.assert_allclose(np.array(r_graph.GetY()),
                            np.array(n_graph.GetY()) / np.array(d_graph.GetY()))

        npt.assert_allclose(np.array(r_graph.GetX()), np.array(n_graph.GetX()))

        x_errs, y_errs = gu.get_errors(r_graph)
        xn_errs, yn_errs = gu.get_errors(n_graph)
        npt.assert_allclose(x_errs, xn_errs)

        _, yd_errs = gu.get_errors(d_graph)

        an_errs  = self._get_analytic_err(np.array(n_graph.GetY()),
                                          np.array(d_graph.GetY()),
                                          yn_errs, yd_errs)

        npt.assert_allclose(y_errs, an_errs)


    def test_tgraph_asymm_errors(self):
        """Check that TGraphAsymmErrors are handled correctly"""
        n_graph = get_random_graph(r.TGraphAsymmErrors, 5)
        d_graph = get_random_graph(r.TGraphAsymmErrors, 5)
        r_graph = gu.divide_graphs(n_graph, d_graph)

        xn_vals, yn_vals = np.array(n_graph.GetX()), np.array(n_graph.GetY())
        yd_vals = np.array(d_graph.GetY())

        x_vals, y_vals = np.array(r_graph.GetX()), np.array(r_graph.GetY())
        npt.assert_allclose(y_vals, yn_vals / yd_vals)
        npt.assert_allclose(x_vals, xn_vals)

        xle, xhe, yle, yhe = gu.get_errors(r_graph)
        xnle, xnhe, ynle, ynhe = gu.get_errors(n_graph)
        _, _, ydle, ydhe = gu.get_errors(d_graph)
        npt.assert_allclose(xle, xnle)
        npt.assert_allclose(xhe, xnhe)

        anl_errs = self._get_analytic_err(yn_vals, yd_vals, ynle, ydle)
        npt.assert_allclose(yle, anl_errs)
        anh_errs = self._get_analytic_err(yn_vals, yd_vals, ynhe, ydhe)
        npt.assert_allclose(yhe, anh_errs)


    def test_tgraph_asymm_errors_correlation(self):
        """
        Check that the uncertainties are calculated correctly with a correlation
        coefficient different from 0
        """
        n_graph = get_random_graph(r.TGraphAsymmErrors, 5)
        d_graph = get_random_graph(r.TGraphAsymmErrors, 5)
        rho = np.random.uniform(-1, 1, 1)
        r_graph = gu.divide_graphs(n_graph, d_graph, rho)

        # Just to be sure that there is no influence on the ratio check it
        # All checks concerning x-direction are not done here!
        y_vals = np.array(r_graph.GetY())
        yn_vals = np.array(n_graph.GetY())
        yd_vals = np.array(d_graph.GetY())
        npt.assert_allclose(y_vals, yn_vals / yd_vals)

        _, _, yle, yhe = gu.get_errors(r_graph)
        _, _, ydle, ydhe = gu.get_errors(d_graph)
        _, _, ynle, ynhe = gu.get_errors(n_graph)

        anl_errs = self._get_analytic_err(yn_vals, yd_vals, ynle, ydle, rho)
        npt.assert_allclose(yle, anl_errs)
        anh_errs = self._get_analytic_err(yn_vals, yd_vals, ynhe, ydhe, rho)
        npt.assert_allclose(yhe, anh_errs)


    def test_tgraph_errors_correlation(self):
        """
        Check that the uncertainties are calculated correctly with a correlation
        coefficient different from 0
        """
        n_graph = get_random_graph(r.TGraphErrors, 5)
        d_graph = get_random_graph(r.TGraphErrors, 5)
        rho = np.random.uniform(-1, 1, 1)
        r_graph = gu.divide_graphs(n_graph, d_graph, rho)

        # Just to be sure that there is no influence on the ratio check it
        # All checks concerning x-direction are not done here!
        y_vals = np.array(r_graph.GetY())
        yn_vals = np.array(n_graph.GetY())
        yd_vals = np.array(d_graph.GetY())
        npt.assert_allclose(y_vals, yn_vals / yd_vals)

        _, y_err = gu.get_errors(r_graph)
        _, yn_err = gu.get_errors(n_graph)
        _, yd_err = gu.get_errors(d_graph)

        an_err = self._get_analytic_err(yn_vals, yd_vals, yn_err, yd_err, rho)
        npt.assert_allclose(y_err, an_err)


if __name__ == '__main__':
    unittest.main()
