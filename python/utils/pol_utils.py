#!/usr/bin/env python
"""
Module containing functions common to polarization studies
"""

import numpy as np

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
