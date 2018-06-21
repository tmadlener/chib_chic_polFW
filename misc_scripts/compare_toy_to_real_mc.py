#!/usr/bin/env python
"""
Script to produce plots comparing toy mc to real mc chic2 / chic1 ratio
"""

from collections import OrderedDict

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
# r.gROOT.SetBatch()


from utils.data_handling import (
    get_dataframe, apply_selections, create_histogram
)
from utils.selection_functions import (
    loose_cuts, loose_muon_sel, fiducial_cuts, fiducial_muon_sel, flat_pt,
    photon_sel, pt_eta_sel, vtx_prob_sel, trigger_sel, chic_state_sel,
    jpsi_kin_sel
)
from utils.plot_helpers import (
    mkplot, default_colors, setup_legend, setup_latex, put_on_latex
)
from utils.hist_utils import divide
from utils.misc_helpers import get_bin_cut_df
from utils.setup_plot_style import set_TDR_style

def single_mu_sel_toy(dfr, cuts, smeared=True):
    """Apply the single muon selection to the toy mc"""
    mu_base = '{}_lep{}'
    mup_pt, mun_pt = mu_base.format('pT', 'P'), mu_base.format('pT', 'N')
    mup_eta, mun_eta = mu_base.format('eta', 'P'), mu_base.format('eta', 'N')

    if smeared:
        mup_pt += '_sm'
        mun_pt += '_sm'
        mup_eta += '_sm'
        mun_eta += '_sm'

    return pt_eta_sel(dfr.loc[:, mup_pt], dfr.loc[:, mup_eta].abs(), cuts) & \
        pt_eta_sel(dfr.loc[:, mun_pt], dfr.loc[:, mun_eta].abs(), cuts)


def photon_sel_toy(dfr, cuts, smeared=True):
    """Photon selection for the toy mc"""
    ppt, peta = 'pT_gamma', 'y_gamma'
    if smeared:
        ppt += '_sm'
        peta += '_sm'

    return pt_eta_sel(dfr.loc[:, ppt], dfr.loc[:, peta].abs(), cuts) & \
        (dfr.loc[:, ppt] < 7) # to be able to apply efficiencies


def jpsi_kin_sel_toy(dfr, min_pt=8, max_pt=20, max_rap=1.2, smeared=True):
    """Jpsi kinematic selection for toy mc"""
    jpt, jrap = 'pT_jpsi', 'y_jpsi'
    if smeared:
        jpt += '_sm'
        jrap += '_sm'

    return get_bin_cut_df(dfr, jpt, min_pt, max_pt) & \
        (dfr.loc[:, jrap].abs() < max_rap)


fiducial_toy = lambda df: single_mu_sel_toy(df, fiducial_cuts())
loose_toy = lambda df: single_mu_sel_toy(df, loose_cuts())
photon_toy = lambda df: photon_sel_toy(df, flat_pt(0.41, 1.5))
photon_mc = lambda df: photon_sel(df, flat_pt(0.41, 1.5))

# basic selection for the real MC (will always be used)
BASIC_SEL_MC = (
    trigger_sel,
    jpsi_kin_sel
)
# same for toy
BASIC_SEL_TOY = (
    jpsi_kin_sel_toy,
)

# selection sets for real MC
SELECTIONS_MC = OrderedDict()
SELECTIONS_MC['loose'] = BASIC_SEL_MC + (loose_muon_sel,)
SELECTIONS_MC['fiducial'] = BASIC_SEL_MC + (fiducial_muon_sel,)
SELECTIONS_MC['loose + photon'] = BASIC_SEL_MC + (loose_muon_sel, photon_mc)
SELECTIONS_MC['fiducial + photon'] = BASIC_SEL_MC + (fiducial_muon_sel, photon_mc)
# same for toy
SELECTIONS_TOY = OrderedDict()
SELECTIONS_TOY['loose'] = BASIC_SEL_TOY + (loose_toy,)
SELECTIONS_TOY['fiducial'] = BASIC_SEL_TOY + (fiducial_toy,)
SELECTIONS_TOY['loose + photon'] = BASIC_SEL_TOY + (loose_toy, photon_toy)
SELECTIONS_TOY['fiducial + photon'] = BASIC_SEL_TOY + (fiducial_toy, photon_toy)

# which variables to plot and how they can be obtained from the dataframe
VARIABLES = {
    'costh_HX': {
        'toy': lambda df: df.loc[:, 'costh_HX_sm'].abs(),
        'mc': lambda df: df.loc[:, 'costh_HX'].abs(),
        'sett': (8, 0, 1)
    }
}

# define which variables to read in
TOY_VARIABLES = [
    '{pT,eta}_lep{N,P}_sm', # smeared single muon kinematics
    'costh_HX_sm',
    '{y,pT}_{gamma,jpsi}_sm', # smeared jpsi and photon kinematics
    '*eff_sm', # smeared efficiencies
]
MC_VARIABLES = [
    'trigger',
    'Jpsi{Pt,Rap}',
    'photon{Pt,Eta}',
    'mu{P,N}{Pt,Eta}',
    'costh_HX',
    'pdgId'
]

PLOT_ATTRIBUTES = {
    'loose': {
        'mc': {'fillalpha': (default_colors()[4], 0.5), 'size': 0, 'marker': 22},
        'toy': {'color': default_colors()[4], 'size': 1.5, 'marker': 22}
    },
    'fiducial': {
        'mc': {'fillalpha': (default_colors()[3], 0.5), 'size': 0, 'marker': 23},
        'toy': {'color': default_colors()[3], 'size': 1.5, 'marker': 23}
    },
    'loose + photon': {
        'mc': {'fillalpha': (default_colors()[0], 0.5), 'size': 0, 'marker': 20},
        'toy': {'color': default_colors()[0], 'size': 1.5, 'marker': 20}
    },
    'fiducial + photon': {
        'mc': {'fillalpha': (default_colors()[1], 0.5), 'size': 0, 'marker': 21},
        'toy': {'color': default_colors()[1], 'size': 1.5, 'marker': 21}
    }
}

LEGPOS = {
    'costh_HX': (0.18, 0.16, 0.36, 0.3)
}
LABELS = {
    'costh_HX': '|cos#vartheta^{HX}|'
}

def get_weight_func(muon=True, photon=True, smeared=True):
    """Get the efficiency weight for an event"""
    peff = 'gamma_eff'
    mupeff, muneff = 'lepP_eff', 'lepN_eff'
    if smeared:
        peff += '_sm'
        mupeff += '_sm'
        muneff += '_sm'

    photon_eff = lambda df: df.loc[:, peff] * 0.01
    muon_eff = lambda df: df.loc[:, mupeff] * df.loc[:, muneff]

    if muon and photon:
        return lambda df: photon_eff(df) * muon_eff(df)
    if muon:
        return muon_eff
    if photon:
        return photon_eff
    return None


def create_hist(dfr, variable, hist_sett, get_weights=None):
    """
    Create a histogram from the passed dataframe
    """
    weights = None
    if get_weights is not None and hasattr(get_weights, '__call__'):
        weights = get_weights(dfr)

    if hasattr(variable, '__call__'):
        var = variable(dfr)
    else:
        var = dfr.loc[:, variable]
    hist = create_histogram(var, hist_sett, weights=weights)

    return hist


def get_ratio(dchi1, dchi2, variable, selections, hist_sett, get_weights=None):
    """
    Get the chic2 / chic1 ratio
    """
    hchi1 = create_hist(apply_selections(dchi1, selections), variable,
                        hist_sett, get_weights)
    hchi2 = create_hist(apply_selections(dchi2, selections), variable,
                        hist_sett, get_weights)

    # scale such that the integrated ratio is 1
    nchi1 = hchi1.Integral()
    nchi2 = hchi2.Integral()

    ratio = divide(hchi2, hchi1)

    ratio.Scale(nchi1 / nchi2)

    return ratio


def get_ratio_mc(dfr, variable, selections, hist_sett):
    """
    Get the chic2 / chic1 ratio for the real mc (necessary since chic1 and
    chic2 are in the same dataframe)
    """
    dchi1 = apply_selections(dfr, lambda df: chic_state_sel(df, 'chic1'))
    dchi2 = apply_selections(dfr, lambda df: chic_state_sel(df, 'chic2'))

    return get_ratio(dchi1, dchi2, variable, selections, hist_sett)


def main(args):
    """Main"""
    mcdata = get_dataframe(args.mcfile, columns=MC_VARIABLES)
    chi1data = get_dataframe(args.chic1file, columns=TOY_VARIABLES)
    chi2data = get_dataframe(args.chic2file, columns=TOY_VARIABLES)

    mc_ratios = OrderedDict()
    toy_ratios = OrderedDict()

    for name in SELECTIONS_TOY:
        toy_ratios[name] = OrderedDict()
        mc_ratios[name] = OrderedDict()

        toy_sel = SELECTIONS_TOY[name]
        mc_sel = SELECTIONS_MC[name]

        for var in VARIABLES:
            hist_sett = VARIABLES[var]['sett']
            toy_ratios[name][var] = get_ratio(chi1data, chi2data,
                                              VARIABLES[var]['toy'], toy_sel,
                                              hist_sett, get_weight_func())
            mc_ratios[name][var] = get_ratio_mc(mcdata, VARIABLES[var]['mc'],
                                                mc_sel, hist_sett)

    set_TDR_style()
    for name in SELECTIONS_TOY:
        for var in VARIABLES:
            leg = setup_legend(*LEGPOS[var])
            can = mkplot(mc_ratios[name][var],
                         attr=[PLOT_ATTRIBUTES[name]['mc']],
                         drawOpt='E2', legOpt='F',
                         xLabel=LABELS[var], yLabel='#chi_{c2} / #chi_{c1}',
                         leg=leg, legEntries=['real mc'])
            can = mkplot(toy_ratios[name][var],
                         attr=[PLOT_ATTRIBUTES[name]['toy']],
                         drawOpt='E1same', legOpt='PLE', can=can,
                         xLabel=LABELS[var], yLabel='#chi_{c2} / #chi_{c1}',
                         leg=leg, legEntries=['toy mc'])
            latex = setup_latex()
            can.add_tobject(latex)
            put_on_latex(latex, [(0.18, 0.96, name)])
            savename = '_'.join(['comp_real_toy', name, var]).replace(' + ', '_')
            can.SaveAs(savename + '.pdf')

    # can.SaveAs('test_mc_comp_toy.pdf')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to produce plot(s) '
                                     'comparing real mc to toy mc')
    parser.add_argument('mcfile', help='file containing the real mc file (raw) '
                        'data')
    parser.add_argument('chic1file', help='file containing the chic1 toy mc '
                        'data')
    parser.add_argument('chic2file', help='file containing the chic2 toy mc '
                        'data')


    clargs = parser.parse_args()
    main(clargs)
