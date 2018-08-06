#!/usr/bin/env python
"""
Script for running costh binned mass fits
"""

import re
import json

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
from utils.roofit_utils import get_var, ws_import, get_args
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
            'chi_mass_rf1S[9.7,10.15]', 'costh_HX[-1, 1]'
        ],
    }
    if state == 'chib' and massmodel:
        return [
            '{}[{}, {}]'.format(massmodel.mname, massmodel.fitvarmin, massmodel.fitvarmax),
            'costh_HX[-1, 1]'
        ]

    return wsvars[state]


def get_shape_params(wsp, savename, mass_model):
    """
    Get the (free) shape parameters of the model and the fit result
    corresponding to the savename
    """
    fit_res = wsp.genobj('fit_res_{}'.format(savename))
    all_params = [p.GetName() for p in get_args(fit_res.floatParsFinal())]
    event_params = mass_model.nevent_vars

    for evpar in event_params:
        all_params.remove(evpar)

    return all_params


def import_data(wsp, datatree, datavars):
    """
    Import all data that is necessary to the workspace
    """
    logging.info('Importing dataset')

    # deactivate all branches and just activate the used ones
    datatree.SetBranchStatus('*', 0)
    varlist = []
    for dvar in datavars:
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


def do_binned_fits(mass_model, wsp, costh_bins, refit=False):
    """
    Do the fits in all costh bins
    """
    for ibin, ctbin in enumerate(costh_bins):
        bin_sel = get_bin_cut_root('TMath::Abs(costh_HX)', *ctbin)
        savename = 'costh_bin_{}'.format(ibin)
        mass_model.fit(wsp, savename, bin_sel)
        if refit:
            # fix all but the event number variables to their current values
            # in the workspace and refit them leaving only the event numbers
            # free (and don't forget to release the shape again afterwards)
            shape_par = get_shape_params(wsp, savename, mass_model)
            mass_model.fix_params(wsp, [(sp, None) for sp in shape_par])
            refit_sn = 'refit_costh_bin_{}'.format(ibin)
            mass_model.fit(wsp, refit_sn, bin_sel)
            mass_model.release_params(wsp, shape_par)


def create_bin_info_json(state, nbins, datafile, fitfile, bininfo_file=None,
                         fixedbinning=[]):
    """
    Determine the costh binning and store the information into a json
    """
    if bininfo_file is None:
        bininfo_file = fitfile.replace('.root', '_bin_sel_info.json')

    treename = 'chic_tuple' if state == 'chic' else 'data'
    dfr = get_dataframe(datafile, treename, columns=['costh_HX'])

    costh_bins = get_costh_binning(dfr, nbins) if not fixedbinning else fixedbinning
    costh_means = get_bin_means(dfr, lambda d: d.costh_HX.abs(), costh_bins)

    bin_sel_info = {
        'costh_bins': costh_bins,
        'costh_means': costh_means,
    }
    with open(bininfo_file, 'w') as info_file:
        json.dump(bin_sel_info, info_file, sort_keys=True, indent=2)

    return costh_bins


def rw_bin_sel_json(bininfo_file, datafile, updated_name=None):
    """
    Read the passed bin info file, get the costh_binning and calculate updated
    costh means for the bins and store the resulting content into a new file
    """
    if updated_name is None:
        updated_name = bininfo_file.replace('.json', '_jpsi.json')

    with open(bininfo_file, 'r') as info_file:
        bin_info = json.load(info_file)

    dfr = get_dataframe(datafile, 'jpsi_tuple', columns=['costh_HX'])

    costh_bins = bin_info['costh_bins']
    costh_means = get_bin_means(dfr, lambda d: d.costh_HX.abs(), costh_bins)
    bin_info['costh_means'] = costh_means

    with open(updated_name, 'w') as info_file:
        json.dump(bin_info, info_file, sort_keys=True, indent=2)

    return costh_bins


def run_fit(model, tree, costh_bins, datavars, outfile, refit=False):
    """Import data, run fits and store the results"""
    wsp = r.RooWorkspace('ws_mass_fit')
    import_data(wsp, tree, datavars)
    model.define_model(wsp)
    wsp.Print()

    do_binned_fits(model, wsp, costh_bins, refit)
    wsp.writeToFile(outfile)


def run_chic_fit(args):
    """Setup everything and run the chic fits"""
    logging.info('Running chic fits')
    cond_mkdir_file(args.outfile)

    model = ChicMassModel('chicMass')
    costh_binning = create_bin_info_json('chic', args.nbins, args.datafile,
                                         args.outfile, args.bin_info)

    dvars = get_ws_vars('chic')
    dataf = r.TFile.Open(args.datafile)
    tree = dataf.Get('chic_tuple')
    run_fit(model, tree, costh_binning, dvars, args.outfile, args.refit)


def run_chib_fit(args):
    """Setup everything and run the chic fits"""
    logging.info('Running chib fits')
    cond_mkdir_file(args.outfile)
    model = ChibMassModel(args.configfile)
    fixedbinning=[]
    if args.fixedbinning != '': 
        fixedbinning = [float(i) for i in args.fixedbinning.split(',')]    
    costh_binning = create_bin_info_json('chib', args.nbins, args.datafile,
                                         args.outfile, args.bin_info, fixedbinning)

    dvars = get_ws_vars('chib', model)
    dataf = r.TFile.Open(args.datafile)
    tree = dataf.Get('data')
    run_fit(model, tree, costh_binning, dvars, args.outfile, args.refit)


def run_jpsi_fit(args):
    """Setup everything and run the chic fits"""
    logging.info('Running jpsi fits')
    cond_mkdir_file(args.outfile)
    model = JpsiMassModel('JpsiMass')
    costh_binning = rw_bin_sel_json(args.costh_bin_file, args.datafile,
                                    args.bin_info)

    dvars = get_ws_vars('jpsi')
    dataf = r.TFile.Open(args.datafile)
    tree = dataf.Get('jpsi_tuple')
    run_fit(model, tree, costh_binning, dvars, args.outfile, args.refit)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script for running costh '
                                     'binned mass fits')

    subparsers = parser.add_subparsers(help='Mode to run', dest='mode')

    # parser for global flags
    global_parser = argparse.ArgumentParser(add_help=False)

    global_parser.add_argument('datafile', help='File containing the flat data '
                               'tuple. NOTE: All events that are in the tuple '
                               'will be used so make sure to have the sample '
                               'only contain what is needed')
    global_parser.add_argument('outfile', help='Output file containing all the '
                               'fit results in a workspace')
    global_parser.add_argument('-n', '--nbins', type=int, default=4,
                               help='number of costh bins')
    global_parser.add_argument('-b', '--bin-info', help='json file containing '
                               'the costh binning information (output). If not '
                               'provided a default name will be used',
                               default=None)
    global_parser.add_argument('--refit', help='Fit the shape parameters and '
                               'the number of events and then rerun the fit '
                               'after fixing the shape parameters',
                               default=False, action='store_true')

    # Add the chic parser
    chic_parser = subparsers.add_parser('chic', description='Run the fits using'
                                     ' chic mass model',
                                     parents=[global_parser])
    chic_parser.set_defaults(func=run_chic_fit)

    # Add the chib parser
    chib_parser = subparsers.add_parser('chib', description='Run the fits using'
                                     ' chib mass model',
                                     parents=[global_parser])
    chib_parser.add_argument('configfile', help='Config file in json format '
                             'for chib model.')
    chib_parser.add_argument('--fixedbinning', help='Use fixed binning defined here '
                             '(separated by commas, e.g. 0.1,0.2,0.4,1)',
                             type=str, default='')
    chib_parser.set_defaults(func=run_chib_fit)

    # Add the jpsi parser
    jpsi_parser = subparsers.add_parser('jpsi', description='Run the fits using'
                                     ' jpsi mass model',
                                     parents=[global_parser])
    jpsi_parser.add_argument('costh_bin_file', help='json file containing the '
                             'costh binning information')
    jpsi_parser.set_defaults(func=run_jpsi_fit)

    clargs = parser.parse_args()
    clargs.func(clargs)
