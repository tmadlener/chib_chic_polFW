#!/usr/bin/env python
"""
Common functionality for scanning / plotting / combining PPDs
"""

from utils.hist_utils import rebin
from utils.plot_helpers import mkplot, default_attributes
from utils.plot_decoration import YLABELS

YLABELS.update({'norm': 'N'})
ATTR = default_attributes(size=0)


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


def plot_lth(ppd):
    """
    Make the 1d plot of the lth ppd
    """
    can = mkplot(ppd, xLabel=YLABELS['lth'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[-1, 1])
    return can


def plot_dlth(ppd):
    """
    Make the 1d plot of the dlth ppd
    """
    ppdmax = ppd.GetBinCenter(ppd.GetMaximumBin())
    can = mkplot(ppd, xLabel=YLABELS['dlth'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[ppdmax - 2, ppdmax + 2])
    return can


def plot_norm(ppd):
    """
    Make the 1d plot of the norm ppd
    """
    can = mkplot(ppd, xLabel=YLABELS['norm'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[0.375, 0.625])
    return can


def plot_dlph(ppd):
    """
    Make the 1d plot of the dlph ppd
    """
    ppdmax = ppd.GetBinCenter(ppd.GetMaximumBin())
    can = mkplot(ppd, xLabel=YLABELS['dlph'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[ppdmax - 0.2, ppdmax + 0.2])
    return can


def plot_lph(ppd):
    """
    Make the 1d plot of the lph ppd
    """
    can = mkplot(ppd, xLabel=YLABELS['lph'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[-1, 1])
    return can
