#!/usr/bin/env python
"""
Some tests for the plot_helpers package
"""

import unittest
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np

from utils.misc_helpers import create_random_str
from utils.graph_utils import get_errors

import utils.plot_helpers as ph

class TestGetXYMinMax(unittest.TestCase):
    def setUp(self):
        self.func = r.TF1('test_func',
                          '[0] * (1 + [1] * x[0] + [2] * x[0]*x[0])', -2, 2)
        self.func.SetParameters(1, 2.0, -2.0)

        self.func_x_max = 2
        self.func_x_min = -2
        self.func_y_min = -11
        self.func_y_max = 1.5

        # Use TGraphAsymmErrors with random errors
        self.graph = r.TGraphAsymmErrors(10, np.linspace(0, 4, 10),
                                         np.sqrt(np.linspace(4, 16, 10)),
                                         np.random.uniform(size=10),
                                         np.random.uniform(size=10),
                                         np.random.uniform(size=10),
                                         np.random.uniform(size=10))
        gxlo, gxhi, gylo, gyhi = get_errors(self.graph)

        graph_x = np.array(self.graph.GetX())
        graph_y = np.array(self.graph.GetY())

        self.graph_x_max = np.max(graph_x + gxhi)
        self.graph_x_min = np.min(graph_x - gxlo)
        self.graph_y_max = np.max(graph_y + gyhi)
        self.graph_y_min = np.min(graph_y - gylo)

        self.hist = r.TH1D(create_random_str(8), '', 10, -3, 1)
        self.hist.Fill(-1)
        # Fill one event into each bin to have non-trivial results
        for i in xrange(12):
            self.hist.Fill(-3.2 + 4.0 / 10 * i)

        self.hist_x_max = 1
        self.hist_x_min = -3
        self.hist_y_max = 2
        self.hist_y_min = 1

        # For 2D histograms we don't need to fill the histogram since the min
        # and max values do not depend on it
        self.hist2d = r.TH2D(create_random_str(8), '', 10, -3, 1, 10, 4, 5)
        self.hist2d_x_max = 1
        self.hist2d_x_min = -3
        self.hist2d_y_max = 5
        self.hist2d_y_min = 4


    def test_get_x_max(self):
        # First test the three types individually, then combined
        self.assertAlmostEqual(self.func_x_max, ph.get_x_max(self.func))
        self.assertAlmostEqual(self.hist_x_max, ph.get_x_max(self.hist))
        self.assertAlmostEqual(self.graph_x_max, ph.get_x_max(self.graph))
        self.assertAlmostEqual(self.hist2d_x_max, ph.get_x_max(self.hist2d))

        overall_max = max([self.graph_x_max, self.func_x_max, self.hist_x_max,
                           self.hist2d_x_max])

        self.assertAlmostEqual(overall_max, ph.get_x_max([self.hist,
                                                          self.func,
                                                          self.graph,
                                                          self.hist2d]))


    def test_get_x_min(self):
        # First test the three types individually, then combined
        self.assertAlmostEqual(self.func_x_min, ph.get_x_min(self.func))
        self.assertAlmostEqual(self.hist_x_min, ph.get_x_min(self.hist))
        self.assertAlmostEqual(self.graph_x_min, ph.get_x_min(self.graph))
        self.assertAlmostEqual(self.hist2d_x_min, ph.get_x_min(self.hist2d))

        overall_min = min([self.graph_x_min, self.func_x_min, self.hist_x_min,
                           self.hist2d_x_min])

        self.assertAlmostEqual(overall_min, ph.get_x_min([self.hist,
                                                          self.func,
                                                          self.graph,
                                                          self.hist2d]))


    def test_get_y_max(self):
        # First test the three types individually, then combined
        self.assertAlmostEqual(self.func_y_max, ph.get_y_max(self.func))
        self.assertAlmostEqual(self.hist_y_max, ph.get_y_max(self.hist))
        self.assertAlmostEqual(self.graph_y_max, ph.get_y_max(self.graph))
        self.assertAlmostEqual(self.hist2d_y_max, ph.get_y_max(self.hist2d))

        overall_max = max([self.graph_y_max, self.func_y_max, self.hist_y_max,
                           self.hist2d_y_max])

        self.assertAlmostEqual(overall_max, ph.get_y_max([self.hist,
                                                          self.func,
                                                          self.graph,
                                                          self.hist2d]))


    def test_get_y_min(self):
        # First test the three types individually, then combined
        self.assertAlmostEqual(self.func_y_min, ph.get_y_min(self.func))
        self.assertAlmostEqual(self.hist_y_min, ph.get_y_min(self.hist))
        self.assertAlmostEqual(self.graph_y_min, ph.get_y_min(self.graph))
        self.assertAlmostEqual(self.hist2d_y_min, ph.get_y_min(self.hist2d))

        overall_min = min([self.graph_y_min, self.func_y_min, self.hist_y_min,
                           self.hist2d_y_min])

        self.assertAlmostEqual(overall_min, ph.get_y_min([self.hist,
                                                          self.func,
                                                          self.graph,
                                                          self.hist2d]))


if __name__ == '__main__':
    unittest.main()
