#!/usr/bin/env python
"""
Module containing some common helper functionality
"""

import re

from utils.misc_helpers import stringify

def get_name(eta, base):
    """
    Get a name encoding the eta limits
    """
    return stringify('_'.join([base, 'eta_{:.1f}_{:.1f}'.format(*eta)]))


def get_eta_range(name):
    """
    Get the eta range back from the name created by the get_name function
    """
    flt_rgx = r'([0-9]+p?[0-9]*)'
    range_rgx = r'_'.join([flt_rgx] * 2)

    match = re.search(range_rgx, name)
    if match:
        return float(match.group(1).replace('p', '.')),\
            float(match.group(2).replace('p', '.'))
