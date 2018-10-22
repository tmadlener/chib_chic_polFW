#!/usr/bin/env python
"""
Script to preselect events (e.g. to feed them to the costh binned mass fits
afterwards)
"""

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

# variables that are in the file after the preselection
# variables necessary for the preslection will also be present!
VARIABLES = [
    '{costh,phi}_{HX,PX,CS}', 'chic{Pt,Rap}'
]

def get_jpsi_sel(jpsi_sel_str):
    """
    Get the jpsi selection function from the passed string
    """
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


def main(args):
    """Main"""
    selections = OrderedDict()

    selections['loose muon'] = sf.loose_muon_sel()
    selections['trigger'] = sf.trigger_sel_(args.trigger)
    selections['vtx prob'] = sf.vtx_prob_sel
    selections['jpsi kin sel'] = get_jpsi_sel(args.jpsi)
    selections['photon sel'] = sf.photon_sel_(sf.flat_pt(0.4, 1.5))
    selections['lifetime cut'] = get_lt_selection(args.mc, 2.5)
    selections['deta cut (MC only)'] = get_deta_sel(args.deta)
    selections['chis mass cut'] = sf.chic_mass_sel # not strictly necessary from a PS point of view

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

    store_dataframe(sel_data, args.outfile, 'chic_tuple')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to preselect ')
    parser.add_argument('infile', help='inputfile')
    parser.add_argument('outfile', help='preselected file')
    parser.add_argument('-t', '--trigger', help='trigger name',
                        default='Dimuon8_Jpsi')
    parser.add_argument('-j', '--jpsi', default='8:20:1.2',
                        help='jpsi selection. format: minpt:maxpt:maxrap')
    parser.add_argument('--mc', help='do mc selection (no lifetime cut)',
                        action='store_true', default=False)
    parser.add_argument('--deta', help='Apply the deta selection (only possible'
                        ' on MC)', action='store_true', default=False)
    parser.add_argument('--hist', help='Produce a histogram that shows the '
                        'effects of the different cuts', default=False,
                        action='store_true')

    clargs = parser.parse_args()
    main(clargs)
