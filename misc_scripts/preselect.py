#!/usr/bin/env python
"""
Script to preselect events (e.g. to feed them to the costh binned mass fits
afterwards)
"""

import pandas as pd
# Avoid the false positive warning of setting a value in a copy of a slice
# when adding 'mQ' to the final data frame
pd.options.mode.chained_assignment = None

import numpy as np
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from collections import OrderedDict

from utils.data_handling import (
    get_dataframe, apply_selections, store_dataframe
)
import utils.selection_functions as sf
from utils.plot_helpers import mkplot
from utils.constants import m_psiPDG


# variables that are in the file after the preselection
# variables necessary for the preslection will also be present!
VARIABLES = [
    '{costh,phi}_{HX,PX,CS}_fold',
    #'chic{Pt,Rap}',
    'chicMass',
    'mumugammaMass', 'JpsiMass',
    'JpsiPt', 'JpsiRap'
]

def get_jpsi_sel(jpsi_sel_str):
    """
    Get the jpsi selection function from the passed string
    """
    if jpsi_sel_str is None:
        return sf.all_sel()

    min_pt, max_pt, max_rap = [float(v) for v in jpsi_sel_str.split(':')]
    return sf.jpsi_kin_sel_(min_pt, max_pt, max_rap)


def get_lt_selection(mc=False, ct_sig=2.5):
    """
    Get the lifetime significance cut
    """
    if mc:
        return sf.all_sel()
    return sf.prompt_sel_(ct_sig)


def get_deta_sel(deta=False):
    """
    For MC additionally already
    """
    if deta:
        return sf.deta_sel
    return sf.all_sel()


def apply_selection_one_by_one(data, selections):
    """
    Apply the selections one by one and return the number of events after each
    selection
    """
    n_events = [sf.get_n_events(data)]
    sel_data = data
    for sel in selections:
        sel_data = apply_selections(sel_data, sel)
        n_events.append(sf.get_n_events(sel_data))

    return sel_data, n_events


def create_store_sel_hist(n_events, selections):
    """
    Create and store the selection histogram
    """
    hist = r.TH1D('sel_hist', '', len(n_events), 0, len(n_events))
    hist.SetStats(0)
    for i, n_ev in enumerate(n_events):
        hist.SetBinContent(i + 1, n_ev)

    axis = hist.GetXaxis()
    bin_labels = ['no selection'] + selections
    for i, label in enumerate(bin_labels):
        axis.SetBinLabel(i + 1, label)
    hist.SetYTitle('# events after cut')
    can = mkplot(hist)
    can.SaveAs('selection_hist.pdf')


def get_muon_sel(mu_sel):
    """
    Get the muon selection
    """
    if mu_sel == 'loose':
        return sf.loose_muon_sel()
    else:
        min_pt = float(mu_sel)
        return sf.single_muon_sel(sf.flat_pt(min_pt, 1.6))


def get_state_sel(is_mc, state):
    """
    Get a (sensible) state selection
    """
    if state == 'all':
        return sf.all_sel()

    if not is_mc:
        if state != 'all':
            logging.warning('Selecting chic1 or chic2 events is only possible'
                            ' on MC')
            return sf.all_sel()
    else:
        return sf.state_sel(state)


def main(args):
    """Main"""
    selections = OrderedDict()

    selections['trigger'] = sf.trigger_sel_(args.trigger)
    selections['muon'] = get_muon_sel(args.muon)
    selections['vtx prob'] = sf.vtx_prob_sel_(args.vtxprob)
    selections['photon sel'] = sf.photon_sel_(sf.flat_pt(0.4, 1.5))
    selections['jpsi kin sel'] = get_jpsi_sel(args.jpsi)
    selections['lifetime cut'] = get_lt_selection(args.mc, args.lifetime)
    selections['deta cut (MC only)'] = get_deta_sel(args.deta)
    # selections['chis mass cut'] = sf.chic_mass_sel # not strictly necessary from a PS point of view
    # To ensure rectangular region in costh-phi
    # selections['costh'] = lambda d: d.costh_HX_fold.abs() < 0.625
    selections['state_sel'] = get_state_sel(args.mc, args.state)

    global VARIABLES
    VARIABLES.extend(sf.collect_requirements(selections.values()))
    if args.deta:
        VARIABLES.append('gen_photonEta')
    if args.mc:
        VARIABLES.append('pdgId')
    VARIABLES = list(set(VARIABLES))

    data = get_dataframe(args.infile, columns=VARIABLES, where='trigger > 0')

    if not args.hist:
        sel_data = apply_selections(data, selections.values())
    else:
        sel_data, n_evts = apply_selection_one_by_one(data, selections.values())
        hist = create_store_sel_hist(n_evts, selections.keys())

    # compute the Q-value based mass
    sel_data.loc[:, 'mQ'] = sel_data.mumugammaMass - sel_data.JpsiMass + m_psiPDG

    store_dataframe(sel_data, args.outfile, 'chic_tuple')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to preselect ')
    parser.add_argument('infile', help='inputfile')
    parser.add_argument('outfile', help='preselected file')
    parser.add_argument('-t', '--trigger', help='trigger name',
                        default='Dimuon8_Jpsi')
    parser.add_argument('-j', '--jpsi', default=None,
                        help='jpsi selection. format: minpt:maxpt:maxrap')
    parser.add_argument('-m', '--muon', default='3.5', type=str,
                        help='Specify the muon selection. Either a number, which'
                        ' will then be used as a flat pt cut or \'loose\' for '
                        'the loose selection.')
    parser.add_argument('-l', '--lifetime', help='Specify the lifetime cut that '
                        'should be used (lifetime significance is used as cut '
                        'variable)', type=float, default=2.5)
    parser.add_argument('-v', '--vtxprob', help='Specify the desired vtx prob '
                        'cut value', default=0.01, type=float)
    parser.add_argument('--mc', help='do mc selection (no lifetime cut)',
                        action='store_true', default=False)
    parser.add_argument('--deta', help='Apply the deta selection (only possible'
                        ' on MC)', action='store_true', default=False)
    parser.add_argument('--hist', help='Produce a histogram that shows the '
                        'effects of the different cuts', default=False,
                        action='store_true')

    state_sel = parser.add_mutually_exclusive_group()
    state_sel.add_argument('--all', action='store_const', dest='state',
                           const='all', help='Use all events (default)')
    state_sel.add_argument('--chic1', action='store_const', dest='state',
                           const='chic1', help='Select only chic1 events '
                           '(MC only)')
    state_sel.add_argument('--chic2', action='store_const', dest='state',
                           const='chic2', help='Select only chic2 events '
                           '(MC only)')
    state_sel.set_defaults(state='all')


    clargs = parser.parse_args()
    main(clargs)
