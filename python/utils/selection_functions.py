#!/usr/bin/env python
"""
Module holding some useful functionality to do selections on pandas.DataFrames

TODO: Make pylint happy (or at least happier) by cleaning up code
"""

import numpy as np

from utils.misc_helpers import get_bin_cut_df
from utils.data_handling import apply_selections

def get_cut_funcs(coords):
    """
    Get a cut function from the passed coordinates in the pt-eta plane
    """
    cut_funcs = []

    max_eta = max([c[1] for c in coords])
    # min_pt = min([c[0] for c in coords])

    # contains function, checking if x is somewhere in between l and h
    # without the need to order lower and upper bounds
    cont = lambda x, l, h: ((x > l) & (x < h)) | ((x > h) & (x < l))

    for c1, c2 in zip(coords[:-1], coords[1:]):
        pt1, eta1 = c1
        pt2, eta2 = c2
        # sorting to not have to worry about the order in which the coordinates have been defined
        # 3 cases
        if pt1 == pt2:
            # to capture the values of pt1, etc by "value" create a new scope in a function
            # See, e.g.: https://stackoverflow.com/a/13355291
            def pt_cut(pt, eta, pt1=pt1, pt2=pt2, eta1=eta1, eta2=eta2):
                return (pt > pt1) & cont(eta, eta1, eta2)
            cut_funcs.append(pt_cut)

        elif eta1 == eta2:
            def eta_cut(pt, eta, pt1=pt1, pt2=pt2, eta1=eta1, eta2=eta2):
                return (eta > eta1) & cont(pt, pt1, pt2) & (eta < max_eta)
            cut_funcs.append(eta_cut)

        elif pt1 != pt2 and eta1 != eta2:
            def lin_int(pt, eta, pt1=pt1, pt2=pt2, eta1=eta1, eta2=eta2):
                int_eta = eta1 + (pt - pt1) * (eta2 - eta1) / (pt2 - pt1)
                return (eta > int_eta) & cont(eta, eta1, eta2)
            cut_funcs.append(lin_int)

        else:
            print('Something went wrong while defining the selection')

    return cut_funcs


def pt_eta_sel(pt, eta, cuts):
    """Select in the 2D pt-eta plane applying the specified cuts"""
    # since get_cut_funcs returns multiple functions whose product is the or of all functions,
    # start with all events as not selected
    decision = np.zeros(pt.shape, dtype=bool)

    cut_funcs = get_cut_funcs(cuts)
    for cut in cut_funcs:
        decision |= cut(pt, eta)

    return decision


def get_gen_name(name, gen=False):
    """Helper function to get the gen level variable name from the passed variable"""
    if gen:
        return 'gen_' + name
    else:
        return name


def single_muon_sel(df, cuts, gen=False):
    """Apply a single muon selection (using the AND of both single muons)"""
    mu_name = get_gen_name('mu', gen)
    mup_pt, mun_pt = [mu_name + charge + 'Pt' for charge in ['P', 'N']]
    mup_eta, mun_eta = [mu_name + charge + 'Eta' for charge in ['P', 'N']]

    return pt_eta_sel(df[mup_pt], df[mup_eta].abs(), cuts) & \
        pt_eta_sel(df[mun_pt], df[mun_eta].abs(), cuts)


def photon_sel(df, cuts, gen=False):
    """Apply a photon selection"""
    phot_name = get_gen_name('photon', gen)
    return pt_eta_sel(df[phot_name + 'Pt'], df[phot_name + 'Eta'].abs(), cuts)


def trigger_sel(df):
    """Apply the trigger selection (reco and data only)"""
    trigger_bit = 1 # for Dimuon8_Jpsi; 2 for Dimuon10_Jpsi (2012 only)
    # NOTE: assuming that the trigger decision is stored in 'trigger' column
    return (df.trigger & trigger_bit) == trigger_bit


def vtx_prob_sel(df, prob=0.01):
    """Select only the events with vtx probability larger than 0.01"""
    return df.vtxProb > prob


def chic_mass_sel(df, min_mass=3.325, max_mass=3.725):
    """Select only those events with a chic mass between min_mass and max_mass"""
    return get_bin_cut_df(df, 'chicMass', min_mass, max_mass)


def jpsi_kin_sel(df, min_pt=8, max_pt=20, max_rap=1.2, gen=False):
    """Kinematic selection of jpsi"""
    jpsi_name = get_gen_name('Jpsi', gen)
    return (get_bin_cut_df(df, jpsi_name + 'Pt', min_pt, max_pt)) & \
        (df[jpsi_name + 'Rap'].abs() < max_rap)


def chic_state_sel(df, state):
    """(On MC) select either chic1 or chic2 using the pdgId"""
    pdgId = {'chic1': 20443, 'chic2': 445}
    return df.pdgId == pdgId[state]


def flat_pt(pt, eta):
    """Helper function to return a flat pt-eta cut in 2D coords"""
    return (pt, 0), (pt, eta)


def fiducial_cuts():
    """Get the coordinates representing the standard fiducial cuts"""
    coords = (
        (3.0, 1.6), (3.0, 1.4), (3.5, 1.4), (3.5, 1.2), (4.5, 1.2),
        (4.5, 0)
    )
    return coords


def loose_cuts():
    """Loosest cut that does not have any acceptance holes in pt-eta"""
    coords = (
        (2.0, 1.6), (3.5, 1.2), (3.5, 0)
    )

    return coords


def fiducial_muon_sel(df, gen=False):
    """
    The fiducial muon selection
    """
    return single_muon_sel(df, fiducial_cuts(), gen)


def loose_muon_sel(df, gen=False):
    """The loose muon selection"""
    return single_muon_sel(df, loose_cuts(), gen)


def get_n_events(data, selections=None, weight=None):
    """
    Get the number of events in the dataframe surviving the passed selection

    If weight is not None the corresponding column will be used as weights,
    unless weight is a function taking the dataframe as only input, than the sum
    of the array returned by that function call will be returned
    """
    sel_data = apply_selections(data, selections)
    if weight is None:
        return sel_data.shape[0]
    if hasattr(weight, '__call__'):
        return weight(sel_data).sum()
    return sel_data[weight].sum()
