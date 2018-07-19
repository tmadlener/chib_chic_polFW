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
    chunks
)
from utils.roofit_utils import get_var, ws_import
from utils.chic_fitting import ChicMassModel
from utils.jpsi_fitting import JpsiMassModel
from utils.chib_fitting import ChibMassModel


def get_ws_vars(state, massmodel=None):
    """
    Get the workspace variables for the given state (jpsi or chic)
    """
    logging.debug('Getting workspace variables for {}'.format(state))
    wsvars = {
        'chic': [
            'chicMass[3.325, 3.725]', 'costh_HX[-1, 1]'
        ],
        'jpsi': [
            'JpsiMass[2.9, 3.25]', 'costh_HX[-1, 1]'
        ],
        'chib' : [
            'chi_mass_rf1S[9.6,10.15]', 'costh_HX[-1, 1]'
        ],
    }
    if state == 'chib' and massmodel:
        return [
            '{}[{}, {}]'.format(massmodel.mname, massmodel.fitvarmin, massmodel.fitvarmax),
            'costh_HX[-1, 1]'
        ]

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

    # The RooArgSet constructor is overloaded up to 9 arguments
    # Split the varlist in chunks of maximum 9 values
    varlists = chunks(varlist, 9)
    var_argset = r.RooArgSet()
    for varl in varlists:
        sub_set = r.RooArgSet(*varl)
        var_argset.add(sub_set)

    data = r.RooDataSet('full_data', 'full_data', datatree, var_argset)
    ws_import(wsp, data)


def do_binned_fits(mass_model, wsp, costh_bins):
    """
    Do the fits in all costh bins
    """
    for ibin, ctbin in enumerate(costh_bins):
        bin_sel = get_bin_cut_root('TMath::Abs(costh_HX)', *ctbin)
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

    # if chic: determine the costh binning and create the pkl file
    if args.state == 'chic' or args.state == 'chib':
        if not args.pklfile:
            pklfile = args.outfile.replace('.root', '_bin_sel_info.pkl')
        else:
            pklfile = args.pklfile

        treename = 'chic_tuple' if args.state == 'chic' else 'data'
        dfr = get_dataframe(args.datafile, treename)

        costh_bins = get_costh_binning(dfr, args.nbins)
        costh_means = get_bin_means(dfr, lambda d: np.abs(d.costh_HX),
                                    costh_bins)

        bin_sel_info = {
            'costh_bins': costh_bins,
            'costh_means': costh_means,
        }
        with open(pklfile, 'w') as pklf:
            pickle.dump(bin_sel_info, pklf)

    # if jpsi: read the pklfile, get the costh information and update the basic
    # selection and store a new pickle file
    else:
        with open(args.pklfile, 'r') as pklf:
            bin_sel_info = pickle.load(pklf)
            costh_bins = bin_sel_info['costh_bins']
        # write updated file
        with open(args.pklfile.replace('.pkl', '_jpsi.pkl'), 'w') as pklf:
            pickle.dump(bin_sel_info, pklf)

    return costh_bins


def main(args):
    """Main"""
    # Make output directory here, since next function wants to write to it
    cond_mkdir_file(args.outfile)
    costh_bins = rw_bin_sel_pklfile(args)

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

    do_binned_fits(mass_model, ws, costh_bins)

    ws.writeToFile(args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script for running costh '
                                     'binned mass fits')
    parser.add_argument('datafile', help='File containing the flat data tuple. '
                        'NOTE: All events that are in the tuple will be used so'
                        ' make sure to have the sample only contain what is '
                        'needed')
    parser.add_argument('outfile', help='Output file containing all the fit '
                        'results in a workspace')
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
