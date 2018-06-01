#!/usr/bin/env python
"""
Module containing some common helper functionality
"""

from utils.misc_helpers import stringify

def get_name(eta, base):
    """
    Get a name encoding the eta limits
    """
    return stringify('_'.join([base, 'eta_{:.1f}_{:.1f}'.format(*eta)]))
