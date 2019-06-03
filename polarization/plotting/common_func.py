#!/usr/bin/env python
"""
Common functionality for scanning / plotting / combining PPDs
"""

from utils.hist_utils import rebin

def get_scaled_ppd(hfile, var, nbins=None):
    """
    Get the ppd scaled to unity
    """
    ppd = hfile.Get('ppd_1d_{}'.format(var))
    if nbins is not None:
        ppd = rebin(ppd, [(0, nbins)])
    ppd.Scale(1 / ppd.Integral())
    return ppd


def get_scaled_ppd_2d(hfile, var1, var2, nbins1=100, nbins2=100):
    """
    Get the 2d ppd scaled to unity
    """
    ppd = hfile.Get('ppd_2d_{}_{}'.format(var1, var2))
    ppd = rebin(ppd, [(0, nbins1), (1, nbins2)])
    ppd.Scale(1 / ppd.Integral())
    return ppd
