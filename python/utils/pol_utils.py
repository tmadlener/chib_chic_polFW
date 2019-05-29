#!/usr/bin/env python
"""
Module containing functions common to polarization studies
"""

import numpy as np
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.misc_helpers import create_random_str

def to_rad(ang):
    """
    Convert the passed angles in degrees to radian

    Args:
        ang (float): angle in degrees

    Returns:
        float: ang in radian
    """
    return ang * np.pi / 180.0


def ang_dist_2d(costh, phi, lam):
    """
    2D angular distribution pdf

    W (costh, phi | lambda) = 3 / (4 * pi * (3 + lth)) * \
        (1 + lth * costh**2 + \
         lph * sinth**2 * cos(2*phi) + \
         ltp * sin(2*th) * cosphi)

    Args:
        costh (float): cos(theta) in a given reference frame
        phi (float): phi in a given reference frame (NOTE: in degrees)
        lam (tuple): 3 element tuple containing lambda_theta, lambda_phi and
            lambda_theta,phi (in this order)

        returns:
            float: Value of W(costh, phi | lambda) at the passed parameters
    """
    # some useful info if you haven't seen this implementation before:
    ## cos^2 + sin^2 = 1
    ## 2 * cos x * sin x = sin 2x
    lth, lph, ltp = lam # unpack into three separate vars
    norm = 3.0 / (4 * np.pi * (3.0 + lth))
    phi_rad = to_rad(phi)
    costh2 = costh * costh
    sinth2 = 1.0 - costh2

    return norm * (1.0 + lth * costh2 + \
                   lph * sinth2 * np.cos(2.0 * phi_rad) + \
                   ltp * 2.0 * costh * np.sqrt(sinth2) * np.cos(phi_rad))


def costh_ratio_1d(costh, lam1, lam2):
    """
    Ratio of W(costh | lam1) / W(costh | lam2) without proper normalization

    R(costh | lam1, lam2) = (1 + lam1 * costh**2) / (1 + lam2 * costh**2)

    Args:
        costh (np.array): cos(theta) for all of the events for which the
            function should be evaluated
        lam1 (np.array): lambda1
        lam2 (np.array): lambda2

    Returns:
       np.array: Values of R(costh | lam1, lam2) at the passed values
    """
    costh2 = costh**2
    return (1 + lam1 * costh2) / (1 + lam2 * costh2)


def ang_dist_lth(costh, lth):
    """
    1D angular distribution pdf in costh

    W(costh | lth) = 3 / (2 * (3 + lth)) * (1 + lth * costh**2)

    Args:
        costh (float): cos(theta) in a given reference frame
        lth (float): lambda_theta parameter

    Returns:
        float: value of W(costh | lth) at the passed parameters
    """
    return 3.0 / (2.0 * (3 + lth)) * (1 + lth * costh * costh)


def w_costh_phi(set_vals=None, fix_norm=False, fix_lambdas=False):
    """
    Angular 2d distribution as TF2, with possibly fixed values.

    Args:
        set_vals (dict, optional): Dictionary containing the values to which
            lth (lamba_theta), lph (lambda_phi), ltp (lambda_theta phi) or norm
            (the normalization) should be set. All values not found will be set
            to 0, except for the normalization that will be set to 1
        fix_norm (boolean, optional): Fix the normalization to the set value (1
            by default)
        fix_lambdas (boolean, optional): Fix the lambdas to the set values (0 by
            default)
    """
    val_to_idx = {'norm': 0, 'lth': 1, 'lph': 2, 'ltp': 3}
    func_str = ('[0]  * 3 / (4 * pi * (3 + [1])) * ('
                '1 + [1] * x[0]*x[0] + '
                '[2] * (1 - x[0]*x[0]) * cos(2 * x[1] * pi / 180) + '
                '[3] * 2 * x[0] * sqrt(1 - x[0]*x[0]) * cos(x[1] * pi / 180)'
                ')')
    func = r.TF2(create_random_str(), func_str, -1, 1, -180, 180)
    for name, idx in val_to_idx.iteritems():
        func.SetParName(idx, name)

    default_vals = {'norm': 1, 'lth': 0, 'lph': 0, 'ltp': 0}
    if set_vals is not None:
        default_vals.update(set_vals)
    for name, val in default_vals.iteritems():
        func.SetParameter(val_to_idx[name], val)

    fix_params = []
    if fix_norm: fix_params.append(0)
    if fix_lambdas: fix_params.extend([1, 2, 3])
    for ipar in fix_params:
        func.FixParameter(ipar, func.GetParameter(ipar))

    return func
