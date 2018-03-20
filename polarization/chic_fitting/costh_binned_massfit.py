#!/usr/bin/env python
"""
Script for running costh binned mass fits
"""

import re
import pickle
import sys

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.data_handling import get_dataframe
from utils.misc_helpers import (
    get_costh_binning, get_bin_means, cond_mkdir_file, get_bin_cut_root,
    get_bin_cut_df
)
from utils.hist_utils import combine_cuts
from utils.roofit_utils import get_var, ws_import
from utils.chic_fitting import ChicMassModel
from utils.jpsi_fitting import JpsiMassModel
from utils.chib_fitting import ChibMassModel

def basic_sel_root(state):
    """
    Get the basic selection in a format that can be processed by root
    """
    logging.debug('Getting basic selection for {}'.format(state))
    cuts = {
        'chic': [
            'chicPt < 990', 'vtxProb > 0.01', 'TMath::Abs(JpsiRap) < 1.2',
            'chicMass > 3.325 && chicMass < 3.725'
        ],
        'jpsi': [
            'vtxProb > 0.01', 'TMath::Abs(JpsiRap) < 1.2',
            'JpsiMass > 2.8 && JpsiMass < 3.3'
        ],
        'chib': ['costh_HX > -1'] # Dummy condition
    }
    return combine_cuts(cuts[state])


def basic_sel_df(dfr, state):
    """
    Get the basic selection in a format suitable for data frames
    """
    logging.debug('Getting basic dataframe selection for {}'.format(state))
    if state == 'chic':
        cut = (dfr.chicPt < 990) & (dfr.vtxProb > 0.01) & \
              (np.abs(dfr.JpsiRap) < 1.2) & \
              (dfr.chicMass > 3.325) & (dfr.chicMass < 3.725)
    elif state == 'jpsi':
        cut = (np.abs(dfr.JpsiRap) < 1.2) & (dfr.vtxProb > 0.01) & \
              (dfr.JpsiMass > 2.8) & (dfr.JpsiMass < 3.3)
    elif state == 'chib':
        cut = np.array([True] * dfr.shape[0])

    return cut


def get_ws_vars(state, massmodel=None):
    """
    Get the workspace variables for the given state (jpsi or chic)
    """
    logging.debug('Getting workspace variables for {}'.format(state))
    wsvars = {
        'chic': [
            'chicMass[3.325, 3.725]', 'costh_HX[-1, 1]', 'Jpsict[-0.1, 1]',
            'vtxProb[0.01, 1]', 'JpsiRap[-1.2, 1.2]', 'JpsiPt[0, 100]',
            'chicPt[0, 990]', 'JpsictErr[0, 1]'
        ],
        'jpsi': [
            'JpsiMass[2.9, 3.25]', 'costh_HX[-1, 1]', 'Jpsict[-0.1, 1]',
            'JpsictErr[0, 1]', 'JpsiRap[-1.2, 1.2]', 'vtxProb[0.01, 1]',
            'JpsiPt[0, 100]'
        ],
        'chib' : [
            'chi_mass_rf1S[9.6,10.15]', 'costh_HX[-1, 1]', 'dimuon_pt[0,100]'
        ],
    }
    if state=='chib' and massmodel:
        return ['{}[{}, {}]'.format(massmodel.mname, massmodel.fitvarmin, massmodel.fitvarmax),
                'costh_HX[-1, 1]', 'dimuon_pt[0,1000]']

    return wsvars[state]


def import_data(wsp, datatree, state, massmodel=None):
    """
    Import all data that is necessary to the workspace
    """
    # definition of variables to be created in the workspace
    # NOTE: some of them act as preselection (e.g. vtxProb and JpsiRap)
    logging.info('Importing dataset')
    dvars = get_ws_vars(state, massmodel)

    # deactivate all branches and just activate the used ones
    datatree.SetBranchStatus('*', 0)
    varlist = []
    for dvar in dvars:
        wsp.factory(dvar)
        branch_var = re.sub(r'([A-Za-z_]*)\[.*\]', r'\1', dvar)
        datatree.SetBranchStatus(branch_var, 1)
        varlist.append(get_var(wsp, branch_var))

    # NOTE: probably breaks at some point with more variables
    var_argset = r.RooArgSet(*varlist)
    data = r.RooDataSet('full_data', 'full_data', datatree, var_argset)
    ws_import(wsp, data)


def get_bin_sel(basic_sel, costh_bin):
    """
    Get the selection string for the passed costh_bin
    """
    costh_sel = get_bin_cut_root('TMath::Abs(costh_HX)', *costh_bin)
    return combine_cuts([basic_sel, costh_sel])


def do_binned_fits(mass_model, wsp, basic_sel, costh_bins):
    """
    Do the fits in all costh bins
    """
    for ibin, ctbin in enumerate(costh_bins):
        bin_sel = get_bin_sel(basic_sel, ctbin)
        savename = 'costh_bin_{}'.format(ibin)
        mass_model.fit(wsp, savename, bin_sel)


def rw_bin_sel_pklfile(args):
    """
    Read and update the pickle file for J/psi or create it for chic
    """
    if args.state == 'jpsi' and not args.pklfile:
        logging.fatal('Need a pickle file containing the costh binning for '
                      'doing jpsi fits')
        sys.exit(1)

    ptvar = 'dimuon_pt' if args.state == 'chib' else 'JpsiPt'
    # get prompt events (state dependency in basic selection)
    basic_sel = combine_cuts([basic_sel_root(args.state),
                              get_bin_cut_root(ptvar, args.ptmin, args.ptmax)])
    if args.state != 'chib' :
        basic_sel = combine_cuts([basic_sel, 'TMath::Abs(Jpsict / JpsictErr) < 2.5'])

    # if chic: determine the costh binning and create the pkl file
    if args.state == 'chic' or args.state == 'chib':
        if not args.pklfile:
            pklfile = args.outfile.replace('.root', '_bin_sel_info.pkl')
        else:
            pklfile = args.pklfile

        treename = 'chic_tuple' if args.state == 'chic' else 'data'
        dfr = get_dataframe(args.datafile, treename)
        basic_sel_bin = (basic_sel_df(dfr, args.state)) & \
                        (get_bin_cut_df(dfr, ptvar, args.ptmin, args.ptmax))
        if args.state == 'chic':
            basic_sel_bin = basic_sel_bin & (np.abs(dfr.Jpsict / dfr.JpsictErr) < 2.5)

        costh_bins = get_costh_binning(dfr, args.nbins, selection=basic_sel_bin)
        costh_means = get_bin_means(dfr, lambda d: np.abs(d.costh_HX),
                                    costh_bins, basic_sel_bin)

        bin_sel_info = {
            'costh_bins': costh_bins,
            'costh_means': costh_means,
            'basic_sel': basic_sel
        }
        with open(pklfile, 'w') as pklf:
            pickle.dump(bin_sel_info, pklf)

    # if jpsi: read the pklfile, get the costh information and update the basic
    # selection and store a new pickle file
    else:
        with open(args.pklfile, 'r') as pklf:
            bin_sel_info = pickle.load(pklf)
            bin_sel_info['basic_sel'] = basic_sel # update info
            costh_bins = bin_sel_info['costh_bins']
        # write updated file
        with open(args.pklfile.replace('.pkl', '_jpsi.pkl'), 'w') as pklf:
            pickle.dump(bin_sel_info, pklf)

    return costh_bins, basic_sel


def main(args):
    """Main"""
    # Make output directory here, since next function wants to write to it
    cond_mkdir_file(args.outfile)
    costh_bins, basic_sel = rw_bin_sel_pklfile(args)

    dataf = r.TFile.Open(args.datafile)
    treename = {
        'chic' : 'chic_tuple',
        'chib' : 'data',
        'jpsi' : 'jpsi_tuple'
        }

    datat = dataf.Get(treename[args.state])

    # create the workspace
    ws = r.RooWorkspace('ws_mass_fit')
    if args.state == 'chic':
        mass_model = ChicMassModel('chicMass')
    elif args.state == 'chib':
        mass_model = ChibMassModel(args.configfile)
    else:
        mass_model = JpsiMassModel('JpsiMass')
    import_data(ws, datat, args.state, mass_model)
    mass_model.define_model(ws)
    ws.Print()

    do_binned_fits(mass_model, ws, basic_sel, costh_bins)

    ws.writeToFile(args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script for running costh '
                                     'binned mass fits')
    parser.add_argument('datafile', help='File containing the flat data tuple')
    parser.add_argument('outfile', help='Output file containing all the fit '
                        'results in a workspace')
    parser.add_argument('--ptmin', type=float, default=8, help='minimum pt')
    parser.add_argument('--ptmax', type=float, default=20, help='maximum pt')
    parser.add_argument('-n', '--nbins', type=int, default=4,
                        help='number of costh bins')
    parser.add_argument('-pf', '--pklfile', help='Pickle file containing the '
                        'costh binning. Required  as input for J/psi. Overrides'
                        ' default name for chic.', type=str, default='')
    parser.add_argument('--configfile', help='Config file in json format for chib model.',
                        type=str, default='config.json')

    state_sel = parser.add_mutually_exclusive_group()
    state_sel.add_argument('--chic', action='store_const', dest='state',
                           const='chic', help='Do mass fits for chic data')
    state_sel.add_argument('--jpsi', action='store_const', dest='state',
                           const='jpsi', help='Do mass fits for jpsi data')
    state_sel.add_argument('--chib', action='store_const', dest='state',
                           const='chib', help='Do mass fits for chib data')
    parser.set_defaults(state='chic')


    clargs = parser.parse_args()

    main(clargs)
