#!/usr/bin/env python
"""
Module containing utility functionality for dealing with TGraph (and
derivatives)
"""

import numpy as np
import ROOT as r

from decorator import decorator
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s: %(message)s')

from utils.misc_helpers import get_vals_from_rwbuffer

@decorator
def out_of_range_default(func, *args):
    """
    Decorator for accessing points of a TGraph 'safely' by emitting a warning
    and returning a nan value

    Returns:
        float: Return value of function taking graph as first argument and index
            as second or nan if out of range
    """
    graph, idx = args[0:2]
    n_points = graph.GetN()
    if n_points < idx:
        logging.error('Cannot access index {} in graph having only {} points'
                      .format(idx, n_points))
        return np.nan

    return func(*args)


def get_errors(graph):
    """
    Get the errors of the passed in graph

    Args:
        graph (ROOT.TGraphErrors or ROOT.TGraphAsymmErrors): graph for which
            errors should be obtained

    Returns:
        tuple: Tuple containing numpy.arrays with all the error values. If
            graph was a TGraphErors two elements are returned (x, y). If graph
            was a TGraphAsymmErrors four elements are returend
            (x_lo, x_hi, y_lo, y_hi)
    """
    n_points = graph.GetN()

    # have to test from specific to less specific here
    if graph.InheritsFrom('TGraphAsymmErrors'):
        x_lo_errs = get_vals_from_rwbuffer(graph.GetEXlow(), n_points)
        x_hi_errs = get_vals_from_rwbuffer(graph.GetEXhigh(), n_points)
        y_lo_errs = get_vals_from_rwbuffer(graph.GetEYlow(), n_points)
        y_hi_errs = get_vals_from_rwbuffer(graph.GetEYhigh(), n_points)
        return x_lo_errs, x_hi_errs, y_lo_errs, y_hi_errs
    if graph.InheritsFrom('TGraphErrors'):
        x_errs = np.array(graph.GetEX())
        y_errs = np.array(graph.GetEY())
        return x_errs, y_errs


def shift_graph(graph, shift):
    """
    Shift the central y-values of a graph

    Args:
        graph (ROOT.TGraph or inheriting) graph to shift
        shift (float): shift value by which the y-value of the graph should be
            shifted

    Returns:
        ROOT.TGraph: graph of the same type as the passed in graph shifted by
            a constant factor
    """
    x_vals = np.array(graph.GetX())
    n_points = graph.GetN()

    y_vals = np.array(graph.GetY()) + shift

    # have to test from specific to less specific here
    if graph.InheritsFrom('TGraphAsymmErrors'):
        return r.TGraphAsymmErrors(n_points, x_vals, y_vals, *get_errors(graph))
    if graph.InheritsFrom('TGraphErrors'):
        return r.TGraphErrors(n_points, x_vals, y_vals, *get_errors(graph))

    return r.TGraph(n_points, x_vals, y_vals)


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
        x_lo_errs, x_hi_errs, y_lo_errs, y_hi_errs = get_errors(graph)
        y_lo_errs *= scale
        y_hi_errs *= scale
        return r.TGraphAsymmErrors(n_points, x_vals, y_vals,
                                   x_lo_errs, x_hi_errs, y_lo_errs, y_hi_errs)

    if graph.InheritsFrom('TGraphErrors'):
        x_errs, y_errs = get_errors(graph)
        y_errs *= scale
        return r.TGraphErrors(n_points, x_vals, y_vals, x_errs, y_errs)

    return r.TGraph(n_points, x_vals, y_vals)


@out_of_range_default
def get_y(graph, point_idx):
    """
    Get the y-value of a point

    Args:
        graph (ROOT.TGraph or inheriting): graph to obtain point from
        point_idx (int): Index of the point in the graph (0 indexed)

    Returns:
        float: y-value of the graph at point_idx or nan (if out of range)
    """
    x, y = r.Double(0), r.Double(0)
    graph.GetPoint(point_idx, x, y)
    return y
