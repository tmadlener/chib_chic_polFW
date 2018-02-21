#!/usr/bin/env python
"""
Script for running costh binned mass fits
"""

import re
import pickle

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
from utils.chic_fitting import chi_mass_model, do_mass_fit


def basic_sel_root():
    """
    Get the basic selection in a format that can be processed by root
    """
    cuts = ['chicPt < 990', 'vtxProb > 0.01', 'TMath::Abs(JpsiRap) < 1.2',
            'chicMass > 3.325 && chicMass < 3.725']
    return combine_cuts(cuts)


def basic_sel_df(dfr):
    """
    Get the basic selection in a format suitable for data frames
    """
    cut = (dfr.chicPt < 990) & (dfr.vtxProb > 0.01) & \
          (np.abs(dfr.JpsiRap) < 1.2) & \
          (dfr.chicMass > 3.325) & (dfr.chicMass < 3.725)
    return cut


def import_data(wsp, datatree):
    """
    Import all data that is necessary to the workspace
    """
    # definition of variables to be created in the workspace
    # NOTE: some of them act as preselection (e.g. vtxProb and JpsiRap)
    logging.info('Importing dataset')
    dvars = ['chicMass[3.325, 3.725]', 'costh_HX[-1, 1]', 'Jpsict[-0.1, 1]',
             'vtxProb[0.01, 1]', 'JpsiRap[-1.2, 1.2]', 'JpsiPt[0, 100]',
             'chicPt[0, 990]', 'JpsictErr[0, 1]']

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


def create_workspace(name, datatree):
    """
    Create the workspace
    """
    wsp = r.RooWorkspace(name)
    import_data(wsp, datatree)
    chi_mass_model(wsp, 'chicMass')

    wsp.Print()

    return wsp


def do_fit_in_bin(wsp, basic_sel, costh_bin, bin_idx):
    """
    Do the fit in the passed costh bin
    """
    logging.info('Doing mass fit in costh bin {}: {:.2f}, {:.2f}'
                 .format(bin_idx, *costh_bin))
    costh_sel = get_bin_cut_root('TMath::Abs(costh_HX)', *costh_bin)
    savename = 'costh_bin_{}'.format(bin_idx)
    bin_sel = combine_cuts([basic_sel, costh_sel])

    do_mass_fit(wsp, savename, bin_sel)


def main(args):
    """Main"""
    dataf = r.TFile.Open(args.datafile)
    datat = dataf.Get('chic_tuple')
    cond_mkdir_file(args.outfile)

    df = get_dataframe(args.datafile, 'chic_tuple')
    basic_sel_bin = (basic_sel_df(df)) & \
                    (get_bin_cut_df(df, 'JpsiPt', args.ptmin, args.ptmax)) & \
                    (np.abs(df.Jpsict / df.JpsictErr) < 2.5)
    costh_bins = get_costh_binning(df, args.nbins, selection=basic_sel_bin)
    costh_means = get_bin_means(df, lambda d: np.abs(d.costh_HX),
                                costh_bins, basic_sel_bin)

    # get prompt events
    basic_sel = combine_cuts([basic_sel_root(), 'TMath::Abs(Jpsict / JpsictErr) < 2.5',
                              get_bin_cut_root('JpsiPt', args.ptmin, args.ptmax)])

    bin_sel_info = {
        'costh_bins': costh_bins,
        'costh_means': costh_means,
        'basic_sel': basic_sel
    }
    with open(args.outfile.replace('.root', '_bin_sel_info.pkl'), 'w') as pklf:
        pickle.dump(bin_sel_info, pklf)

    ws = create_workspace('ws_mass_fit', datat)
    for ibin, ct_bin in enumerate(costh_bins):
        do_fit_in_bin(ws, basic_sel, ct_bin, ibin)

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

    clargs = parser.parse_args()
    main(clargs)
