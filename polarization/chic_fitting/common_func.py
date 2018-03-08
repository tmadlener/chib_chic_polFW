#!/usr/bin/env python
"""
Module containing functionality shared by scripts in this folder
"""

import pickle

def get_bin_sel_info(pklfile, fitfile):
    """
    Get the binning and selection info
    """
    if not pklfile:
        pklfile = fitfile.replace('.root', '_bin_sel_info.pkl')

    with open(pklfile, 'r') as pklf:
        bin_sel_info = pickle.load(pklf)

    return bin_sel_info
