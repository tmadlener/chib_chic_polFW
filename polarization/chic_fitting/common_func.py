#!/usr/bin/env python
"""
Module containing functionality shared by scripts in this folder
"""

import json

def get_bin_sel_info(jsonfile, fitfile):
    """
    Get the binning and selection info
    """
    if not jsonfile:
        jsonfile = fitfile.replace('.root', '_bin_sel_info.json')

    with open(jsonfile, 'r') as jsonf:
        bin_sel_info = json.load(jsonf)

    return bin_sel_info
