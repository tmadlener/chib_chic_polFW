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


def divide_graphs(ngraph, dgraph):
    """
    Divide the two graphs and return the ratio graph (only in y-direction)

    The x-values are taken from the numerator graph, the ones from the
    denominator graph is ignored.
    The graphs have to be of the same type

    Args:
        ngraph (ROOT.TGraph or inheriting): numerator graph
        dgraph (ROOT.TGraph or inheriting): denominator graph

    Returns:
        ROOT.TGraph: ratio of the numerator to denominator graph (in
            y-direction). If an error occurs, None is returned
    """
    n_points = ngraph.GetN()
    if n_points != dgraph.GetN():
        logging.error('Cannot divide graphs with different numbers of points! '
                      'Number of points: {} / {}'
                      .format(n_points, dgraph.GetN()))
        return None
    if ngraph.ClassName() != dgraph.ClassName():
        logging.error('Cannot divide graphs of different types! '
                      'Types: {} / {}'.format(ngraph.ClassName(),
                                              dgraph.ClassName()))
        return None

    nx_vals = np.array(ngraph.GetX())
    ny_vals = np.array(ngraph.GetY())
    dy_vals = np.array(dgraph.GetY())

    ry_vals = ny_vals / dy_vals

    if ngraph.InheritsFrom('TGraphAsymmErrors'):
        nx_lo_err, nx_hi_err, ny_lo_err, ny_hi_err = get_errors(ngraph)
        _, _, dy_lo_err, dy_hi_err = get_errors(dgraph)

        r_lo_err = ry_vals * \
                   np.sqrt((ny_lo_err / ny_vals)**2 + (dy_lo_err / dy_vals)**2)
        r_hi_err = ry_vals * \
                   np.sqrt((ny_hi_err / ny_vals)**2 + (dy_hi_err / dy_vals)**2)
        return r.TGraphAsymmErrors(n_points, nx_vals, ry_vals,
                                   nx_lo_err, nx_hi_err, r_lo_err, r_hi_err)

    if ngraph.InheritsFrom('TGraphErrors'):
        nx_errs, ny_errs = get_errors(ngraph)
        _, dy_errs = get_errors(dgraph)

        r_errs = np.sqrt((ny_errs / ny_vals)**2 + (dy_errs / dy_vals)**2)
        return r.TGraphErrors(n_points, nx_vals, ry_vals, nx_errs, r_errs)

    return r.TGraph(n_points, nx_vals, ry_vals)


def get_binning(graph):
    """
    Get the binning along the x direction by calculating the bins from the
    central values and the errors

    Args:
        graph (ROOT.TGraph[Asymm]Errors): graph for which the binning is
            desired

    Returns:
        np.array or None: 2D array with the bins (low, high) as columns for the
            different bins. None if not a TGraphErrors or a TGraphAsymmErrors
    """
    x_vals = np.array(graph.GetX())

    if graph.InheritsFrom('TGraphAsymmErrors'):
        x_lo_errs, x_hi_errs, _, _ = get_errors(graph)
        bins_lo = x_vals - x_lo_errs
        bins_hi = x_vals + x_hi_errs
    elif graph.InheritsFrom('TGraphErrors'):
        x_errs, _ = get_errors(graph)
        bins_lo = x_vals - x_errs
        bins_hi = x_vals + x_errs
    else:
        logging.warning('Trying to obtain binning from something else than a '
                        'TGraphErrors or a TGraphAsymmErrors')
        return None # without uncertainties there are no bins

    return np.array(zip(bins_lo, bins_hi))
