#!/usr/bin/env python
"""
Script to produce some overview plots and store them into a root file
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.data_handling import get_dataframe, apply_selections
from utils.selection_functions import (
    fiducial_muon_sel, loose_muon_sel, photon_sel, jpsi_kin_sel, flat_pt,
    get_gen_name
)
from utils.hist_utils import create_histogram

def get_name(name, smeared):
    """
    Get the name of the variable respecting if they are smeared or not
    """
    if smeared:
        return name + '_sm'
    return name


# variables that will only be produced as they are here
SPECIAL_VARS = {
    'chicMass': ('chicMass', (100, 3.2, 4.0)),
    'JpsiMass': ('JpsiMass', (100, 2.6, 3.6)),
}
# variables that will be produced as they are here PLUS their generated version
SM_VARS = {
    'phi_HX': ('phi_HX', (20, -180, 180)),
    'JpsiPt': ('JpsiPt', (140, 0, 70)),
    'chicPt': ('chicPt', (140, 0, 70)),
    'photonPt': ('photonPt', (80, 0, 20)),
    'photonEta': ('photonEta', (40, -4, 4)),
    'muPPt': ('muPPt', (90, 0, 30)),
    'muNPt': ('muNPt', (90, 0, 30)),
    'muPEta': ('muPEta', (40, -4, 4)),
    'muNEta': ('muNEta', (40, -4, 4)),
    'costh_HX': (lambda df: df.loc[:, 'costh_HX'].abs(), (8, 0, 1)),
}
# generated version
GEN_VARS = {
    get_gen_name(v, True): (get_gen_name(v, True), s[1])
    for v, s in SM_VARS.iteritems()
}
GEN_VARS['gen_costh_HX'] = (lambda df: df.loc[:, 'gen_costh_HX'].abs(),
                            (8, 0, 1))

SELECTION_SETS = {
    'gen': None,
    # 'reco': ('reco_muon', 'reco_photon'),
    # 'jpsi_only': ('reco_muon', 'reco_photon', 'jpsi'),
    # 'jpsi_photon': ('reco_muon', 'reco_photon', 'jpsi', 'photon'),
    'jpsi_photon_loose': ('reco_muon', 'reco_photon', 'jpsi',
                          'photon', 'loose'),
    'jpsi_photon_fiducial': ('reco_muon', 'reco_photon', 'jpsi',
                             'photon', 'fiducial'),
}


def get_selection(cuts, gen):
    """
    Get the selection either for smeared or gen level variables
    """
    if cuts is None:
        return None

    jpsi = lambda df: jpsi_kin_sel(df, 8, 20, 1.2, gen)
    photon = lambda df: photon_sel(df, flat_pt(0.41, 1.5), gen)
    loose = lambda df: loose_muon_sel(df, gen)
    fiducial = lambda df: fiducial_muon_sel(df, gen)
    reco_photon = lambda df: df.loc[:, get_name('gamma_eff', not gen)] >= 0
    reco_muon = lambda df: (df.loc[:, get_name('lepP_eff', not gen)] >= 0) &\
                (df.loc[:, get_name('lepN_eff', not gen)] >= 0)

    poss_selections = {
        'jpsi': jpsi,
        'photon': photon,
        'fiducial': fiducial,
        'loose': loose,
        'reco_photon': reco_photon,
        'reco_muon': reco_muon
    }

    return [poss_selections[c] for c in cuts]


def create_store_dist(outfile, df, variable, hist_sett, get_weights=None,
                      **kwargs):
    """
    Create the distribution and store it to the output file
    """
    weights = None
    if get_weights is not None:
        weights = get_weights(df)

    if hasattr(variable, '__call__'):
        var = variable(df)
    else:
        var = df.loc[:, variable]
    hist = create_histogram(var, hist_sett, weights=weights,
                            **kwargs)

    outfile.cd()
    hist.Write()


def store_all_plots(outfile, df, selections, variables, store_base, gen,
                    get_effs=None):
    """
    Store all variable distributions for a given selection
    """
    sel_data = apply_selections(df, get_selection(selections, gen))

    for var_name, settings in variables.iteritems():
        store_name = '_'.join([store_base, var_name])
        ## debugging check
        if outfile.Get(store_name):
            print '{} is already present in file'.format(store_name)
        create_store_dist(outfile, sel_data, settings[0], settings[1],
                          get_weights=get_effs, name=store_name)


def main(args):
    """Main"""
    toy_data = get_dataframe(args.toydata)
    outfile = r.TFile.Open(args.outfile, 'recreate')

    # combined efficiencies calculated on smeared variables
    # NOTE: photon efficiencies are parametrized in percent
    comb_effs = lambda df: df.gamma_eff_sm * 0.01 * \
                df.lepP_eff_sm * df.lepN_eff_sm

    gen_effs = lambda df: df.gamma_eff * 0.01 * \
               df.lepP_eff * df.lepN_eff

    ## store gen level vars
    store_all_plots(outfile, toy_data, SELECTION_SETS['gen'], GEN_VARS,
                    'gen', True)
    ## store smeared vars at generation level
    store_all_plots(outfile, toy_data, SELECTION_SETS['gen'], SM_VARS,
                    'gen', False)
    store_all_plots(outfile, toy_data, SELECTION_SETS['gen'], SPECIAL_VARS,
                    'gen', False)

    ## store smeared level variables for different selection sets
    for sel in SELECTION_SETS:
        if sel == 'gen':
            continue # 'gen' already stored above
        store_base = sel
        for var_set in [SM_VARS, SPECIAL_VARS]:
            # store once without efficiencies
            store_all_plots(outfile, toy_data, SELECTION_SETS[sel], var_set,
                            store_base, False)
            ## ... and once with efficiencies
            store_all_plots(outfile, toy_data, SELECTION_SETS[sel], var_set,
                            store_base + '_eff', False, comb_effs)

        # store gen level selected and efficiency weighted histograms
        store_all_plots(outfile, toy_data, SELECTION_SETS[sel], SM_VARS,
                        'gen_' + store_base, True)
        store_all_plots(outfile, toy_data, SELECTION_SETS[sel], SM_VARS,
                        'gen_' + store_base + '_eff', True, gen_effs)

    outfile.Write()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to generate and store '
                                     'some overview distributions to a root '
                                     'file')
    parser.add_argument('toydata', help='toy data file')
    parser.add_argument('outfile', help='output file to which the created plots'
                        'are stored')

    clargs = parser.parse_args()
    main(clargs)
