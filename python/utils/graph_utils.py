#!/usr/bin/env python
"""
Module containing utility functionality for dealing with TGraph (and
derivatives)
"""

import numpy as np
import ROOT as r

from utils.misc_helpers import get_vals_from_rwbuffer


def scale_graph(graph, scale):
    """
    Scale the graph by a constant factor (only in y-direction)

    Args:
        graph (ROOT.TGraph or inheriting): graph to rescale
        scale (float): factor by which the graph should be rescaled

    Returns:
        ROOT.TGraph: graph of the same type as the passed in graph scaled by a
            factor
    """
    x_vals = np.array(graph.GetX())
    n_points = graph.GetN()

    y_vals = np.array(graph.GetY()) * scale

    # have to test from very specific to less specific here
    if graph.InheritsFrom('TGraphAsymmErrors'):
        x_lo_errs = get_vals_from_rwbuffer(graph.GetEXlow(), n_points)
        x_hi_errs = get_vals_from_rwbuffer(graph.GetEXhigh(), n_points)
        y_lo_errs = get_vals_from_rwbuffer(graph.GetEYlow(), n_points) * scale
        y_hi_errs = get_vals_from_rwbuffer(graph.GetEYhigh(), n_points) * scale
        return r.TGraphAsymmErrors(n_points, x_vals, y_vals,
                                   x_lo_errs, x_hi_errs, y_lo_errs, y_hi_errs)

    if graph.InheritsFrom('TGraphErrors'):
        x_errs = np.array(graph.GetEX())
        y_errs = np.array(graph.GetEY()) * scale
        return r.TGraphErrors(n_points, x_vals, y_vals, x_errs, y_errs)

    return r.TGraph(n_points, x_vals, y_vals)
