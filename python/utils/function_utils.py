#!/usr/bin/env python
"""
Module containing some TF1, etc. helper functions
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.misc_helpers import make_iterable

def _get_y_min_func(func):
    """
    Get the minimum y value of the passed function(s)
    """
    return min(f.GetMinimum() for f in make_iterable(func))


def _get_y_max_func(func):
    """
    Get the maximum y value of the passed function(s)
    """
    return max(f.GetMaximum() for f in make_iterable(func))


def _get_x_max_func(func):
    """
    Get the maximum x value of the passed function(s)
    """
    return max(f.GetXmax() for f in make_iterabel(func))


def _get_x_min_func(func):
    """
    Get the minimum x value of the passed function(s)
    """
    return min(f.GetXmin() for f in make_iterabel(func))
