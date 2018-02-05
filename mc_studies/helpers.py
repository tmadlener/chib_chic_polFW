#!/usr/bin/env python
"""
Module containing helper functions usable in mc studies
"""

import re

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.hist_utils import combine_cuts, draw_var_to_hist, set_hist_opts
from utils.misc_helpers import stringify

def get_hist_settings(var):
    """
    TODO
    """
    hist_settings = {
        'costh': {'n_bins': 16, 'min': 0, 'max': 1},
        'chicMass': {'n_bins': 200, 'min': 3.2, 'max': 3.65},
        'phi': {'n_bins': 20, 'min': 0, 'max': 90}
    }

    if 'costh' in var or 'cosalpha' in var:
        return hist_settings['costh']
    if 'phi' in var:
        return hist_settings['phi']

    return hist_settings[var]


def get_plot_varname(var):
    """
    Get the name from the variable that can be used in plot and histogram names.

    TODO
    """
    s_var = re.sub(r'TMath::Abs\((.*)\)', r'abs_\1', var)
    return s_var


def var_selection(var, v_min, v_max):
    """
    Get events in a min, max window for a given variable

    Args:
        var (str): branch name of the desired variable
        v_min, v_max (float): min and max value of the variable

    Returns:
        str: Selection string to be used in TTree::Draw()
    """
    low_cut = ' > '.join([var, str(v_min)])
    up_cut = ' < '.join([var, str(v_max)])

    return combine_cuts([low_cut, up_cut])


def select_trigger(brname='trigger', bit=2):
    """
    Select triggered events:

    Args:
        brname (str, optional): Name of the branch under which trigger info is
            stored
        bit (int, optional): bit that has to be set in order for being sleeted

    Returns
        str: Selection string to be used in TTree::Draw()
    """
    return '({0} & {1}) == {1}'.format(brname, bit)


def base_selection():
    """
    Get a basic selection that should be applied to all plots

    NOTE: Calls gets the histogram settings for 'chicMass' to determine min and
    max mass cuts

    Returns:
        str: Selection string to be used in TTree::Draw()
    """
    mass_hist_settings = get_hist_settings('chicMass')
    mass_sel = var_selection('chicMass', mass_hist_settings['min'],
                             mass_hist_settings['max'])
    vtx_prob_sel = 'vtxProb > 0.01'

    return combine_cuts([mass_sel, vtx_prob_sel])


def get_histogram(rfile, var, state, selection, treename='chic_mc_tuple'):
    """
    Get histogram of given state for variable and selection.
    TODO
    """
    tree = rfile.Get(treename)
    hist_set = get_hist_settings(var)

    hist_name = '_'.join([stringify(selection), state, get_plot_varname(var)])

    hist = r.TH1D(hist_name, '', hist_set['n_bins'],
                  hist_set['min'], hist_set['max'])

    set_hist_opts(hist)
    state_sel = {'chic1': 'wChic1', 'chic2': 'wChic2'}[state]

    draw_var_to_hist(tree, hist, var, selection, state_sel)

    return hist
