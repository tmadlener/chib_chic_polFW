#!/usr/bin/env python
"""
Determine the smearing histograms from MC
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.data_handling import (
    get_dataframe, apply_selections, create_histogram
)
from utils.selection_functions import (
    trigger_sel, photon_sel, fiducial_muon_sel, jpsi_kin_sel, flat_pt,
    chic_mass_sel
)
from utils.misc_helpers import get_equi_pop_bins, flatten


def get_rel_residuals(dfr, var):
    """
    Get the relative residuals for a given variable where the residuals are
    calculated as reco vs gen and then are returned relative to gen
    """
    gen_var = 'gen_' + var
    return (dfr.loc[:, var] - dfr.loc[:, gen_var]) / dfr.loc[:, gen_var]



def add_rel_residuals(dfr, particles, variables):
    """
    Add the relative residuals for all passed combinations of particles and
    variables
    """
    for part in particles:
        for var in variables:
            full_var = ''.join([part, var])
            res_name = '_'.join([full_var, 'res', 'rel'])
            dfr.loc[:, res_name] = get_rel_residuals(dfr, full_var)


def create_res_v_gen_map(dfr, var, hist_sett):
    """
    Create the 2d (relative) residuals vs gen-level map that can be used in the
    SmearingProvider class
    """
    res_var = var + '_res_rel'
    gen_var = 'gen_' + var
    res_map = create_histogram(np.array([dfr.loc[:, gen_var].abs(),
                                         dfr.loc[:, res_var]]).T,
                               hist_sett)
    return res_map


def create_photon_res_map(dfr, n_bins_res, n_bins_p):
    """
    Create the photon residuals map for a given number of resolution bins (equal
    spacing) and a given number of P(x,y,z) bins (equal population)
    """
    # determine the P(x,y,x) binning
    px_bins = get_equi_pop_bins(dfr, lambda df: df.gen_photonPx.abs(), n_bins_p)
    # get_equi_pop_bins returns bins in an unsuitable format transform them to
    # proper format
    px_bins = np.array(sorted({b for b in flatten(px_bins)}))
    res_bins = np.linspace(-2, 2, n_bins_res)

    hist_sett = (len(px_bins) - 1, px_bins, len(res_bins) - 1, res_bins)

    px_map = create_res_v_gen_map(dfr, 'photonPx', hist_sett)
    py_map = create_res_v_gen_map(dfr, 'photonPy', hist_sett)
    pz_map = create_res_v_gen_map(dfr, 'photonPz', hist_sett)

    photon_map = px_map.Clone()
    photon_map.Add(py_map)
    photon_map.Add(pz_map)

    return photon_map


def create_muon_pxy_res_map(dfr, n_bins_res, n_bins_p):
    """
    Create the muon px and py combined residual maps (see photon for binning)
    """
    px_bins = get_equi_pop_bins(dfr, lambda df: df.gen_muPPx.abs(), n_bins_p)
    px_bins = np.array(sorted({b for b in flatten(px_bins)}))

    res_bins = np.linspace(-0.15, 0.15, n_bins_res)

    hist_sett = (len(px_bins) - 1, px_bins, len(res_bins) - 1, res_bins)

    px_p_map = create_res_v_gen_map(dfr, 'muPPx', hist_sett)
    px_n_map = create_res_v_gen_map(dfr, 'muNPx', hist_sett)
    py_p_map = create_res_v_gen_map(dfr, 'muPPy', hist_sett)
    py_n_map = create_res_v_gen_map(dfr, 'muNPy', hist_sett)

    muon_map = px_n_map.Clone()
    for mu_map in [px_n_map, py_p_map, py_n_map]:
        muon_map.Add(mu_map)

    return muon_map


def create_muon_pz_res_map(dfr, n_bins_res, n_bins_p):
    """
    See above, same thing for muon pz
    """
    pz_bins = get_equi_pop_bins(dfr, lambda df: df.gen_muPPz.abs(), n_bins_p)
    pz_bins = np.array(sorted({b for b in flatten(pz_bins)}))

    res_bins = np.linspace(-0.15, 0.15, n_bins_res)

    hist_sett = (len(pz_bins) - 1, pz_bins, len(res_bins) - 1, res_bins)

    pz_p_map = create_res_v_gen_map(dfr, 'muPPz', hist_sett)
    pz_n_map = create_res_v_gen_map(dfr, 'muNPz', hist_sett)

    muon_map = pz_p_map.Clone()
    muon_map.Add(pz_n_map)

    return muon_map


data = get_dataframe('/afs/hephy.at/data/tmadlener01/ChicPol/Chic2012/InputFiles/NewMCGen2012/flat_tuples/pGunRunI-chic-lowPt_new_tupling_add_more_vars.root')

add_rel_residuals(data, ['photon', 'muP', 'muN'], ['Px', 'Py', 'Pz'])

selections = (
    trigger_sel,
    # lambda df: photon_sel(df, flat_pt(0.4, 1.5), gen=True),
    jpsi_kin_sel,
    fiducial_muon_sel,
    lambda df: photon_sel(df, flat_pt(0.1, 1.5)),
    # lambda df: chic_mass_sel(df, 3.325, 3.725)
    lambda df: chic_mass_sel(df, 0, 3.65)
)

sel_data = apply_selections(data, selections)

outfile = r.TFile('res_maps.root', 'recreate')
photon_map = create_photon_res_map(sel_data, 100, 8)
photon_map.SetName('photon_rel_res_map')

muon_xy_map = create_muon_pxy_res_map(sel_data, 100, 4)
muon_xy_map.SetName('muon_xy_rel_res_map')

muon_z_map = create_muon_pz_res_map(sel_data, 100, 4)
muon_z_map.SetName('muon_z_rel_res_map')

photon_map.Write()
muon_xy_map.Write()
muon_z_map.Write()
outfile.Close()
