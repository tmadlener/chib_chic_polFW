#!/usr/bin/env python
"""
Common functionality for scanning / plotting / combining PPDs
"""

import numpy as np

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


def plot_ltilde(ppd):
    """
    Make the 1d plot of the lambda_tilde ppd
    """
    can = mkplot(ppd, xLabel=YLABELS['ltilde'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[-1, 1])
    return can


def plot_dltilde(ppd):
    """
    Make the 1d plot of the delta_lambdatilde ppd
    """
    ppdmax = ppd.GetBinCenter(ppd.GetMaximumBin())
    can = mkplot(ppd, xLabel=YLABELS['dltilde'], yLabel='PPD [a.u.]',
                 drawOpt='hist', attr=ATTR, xRange=[ppdmax -2, ppdmax + 2])
    return can



def cond_chi1_lambdas(lth, lph, ltp):
    """
    Check whether the values of lth, lph and ltp satisfy the condition of
    Eq. 27 of PRD 83, 096001 (2011)
    """
    cond1 = (lth >= -1./3.) & (lth <= 1)
    cond2 = np.abs(lph) <= ((1 - lth) * 0.25)
    cond3 = (2.25 * (lth - 1./3.)**2 + 6 * ltp**2) <= 1
    cond4 = np.abs(ltp) <= np.sqrt(3) * 0.5 * (lph + 1./3.)
    cond5 = ((lph > 1./9.) & (( (6 * lph - 1)**2 + 6 * ltp**2 ) <= 1)) | (lph <= 1./9.)

    return cond1 & cond2 & cond3 & cond4 & cond5


def cond_chi2_lambdas(lth, lph, ltp):
    """
    Check whether the values of lth, lph and ltp satisfy the condition of
    Eq. 28 of PRD 83, 096001 (2011)
    """
    return (0.3125 * (lth - 0.2)**2 + lph**2 + ltp**2) <= 0.2


def cond_chi1_lth_lph(lth, lph):
    """
    Check whether the values of lth and lph satisfy the 2d condition of Eq. 27
    of PRD 83, 096001 (2011)
    """
    return (lth >= -1./3.) & (lth <= 1) & (np.abs(lph) <= ((1 - lth) * 0.25))


def cond_chi2_lth_lph(lth, lph):
    """
    Check whether the values of lth and lph satisfy the 2d condition of Eq. 28
    of PRD 83, 096001 (2011)
    """
    return (0.3125 * (lth - 0.2)**2 + lph**2) <= 0.2
