#!/usr/bin/env python
"""
Script for running costh binned mass fits
"""

import re
import json
import os
import sys

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import shutil

from utils.data_handling import get_dataframe
from utils.misc_helpers import (
    get_costh_binning, get_bin_means, cond_mkdir_file, get_bin_cut_root,
    chunks, parse_binning
)
from utils.roofit_utils import get_var, ws_import, get_args
from utils.chic_fitting import ChicMassModel
from utils.jpsi_fitting import JpsiMassModel
from utils.chib_fitting import ChibMassModel
from utils.config_fitting import ConfigFitModel


def get_ws_vars(state, massmodel=None, weights=None):
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

    if weights is not None:
        wsvars[state] += ['{}[0, 1e5]'.format(weights)]

    return wsvars[state]


def get_shape_params(wsp, savename, mass_model):
    """
    Get the (free) shape parameters of the model and the fit result
    corresponding to the savename
    """
    fit_res = wsp.genobj('fit_res_{}'.format(savename))
    all_params = [p.GetName() for p in get_args(fit_res.floatParsFinal())]
    event_params = mass_model.nevent_yields

    for evpar in event_params:
        # it is possible that the event yields are somehow fixed via another
        # variable. In that case they will not appear as floating parameters and
        # will also not be in the all_params list.
        # NOTE: In that case it is the responsibility of the user to make sure
        # that the parameters governing the yields are not accidentally fixing
        # it somehow.
        if evpar in all_params:
            all_params.remove(evpar)

    for float_par in mass_model.floating_costh:
        if not float_par in all_params:
            logging.warn('{} is not in the list of parameters for the model.'
                         'It is not possible to leave it floating'
                         .format(float_par))
        else:
            all_params.remove(float_par)

    return all_params


def import_data(wsp, datatree, datavars, weights=None):
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

    if weights is None:
        data = r.RooDataSet('full_data', 'full_data', datatree, var_argset)
    else:
        data = r.RooDataSet('full_data', 'full_data', datatree, var_argset, '',
                            weights)

    ws_import(wsp, data)


def do_binned_fits(mass_model, wsp, costh_bins, refit=False, weighted_fit=False):
    """
    Do the fits in all costh bins
    """
    for ibin, ctbin in enumerate(costh_bins):
        logging.info('Running fit for bin {}: {:.3f} < |costh| < {:.3f}'
                     .format(ibin, *ctbin))
        bin_sel = get_bin_cut_root('TMath::Abs(costh_HX)', *ctbin)
        savename = 'costh_bin_{}'.format(ibin)
        mass_model.fit(wsp, savename, bin_sel, weighted_fit)
        if refit:
            # fix all but the event number variables to their current values
            # in the workspace and refit them leaving only the event numbers
            # free (and don't forget to release the shape again afterwards)
            shape_par = get_shape_params(wsp, savename, mass_model)
            mass_model.fix_params(wsp, [(sp, None) for sp in shape_par])
            refit_sn = 'refit_costh_bin_{}'.format(ibin)
            mass_model.fit(wsp, refit_sn, bin_sel, weighted_fit)
            mass_model.release_params(wsp, shape_par)


def create_bin_info_json(state, bin_str, datafile, fitfile, bininfo_file=None):
    """
    Determine the costh binning and store the information into a json
    """
    if bininfo_file is None:
        bininfo_file = fitfile.replace('.root', '_bin_sel_info.json')

    treename = 'chic_tuple' if state == 'chic' else 'data'
    dfr = get_dataframe(datafile, treename, columns=['costh_HX'])

    auto_match = re.match(r'auto:(\d+)$', bin_str)
    if auto_match:
        costh_bins = get_costh_binning(dfr, int(auto_match.group(1)))
    else:
        binning = parse_binning(bin_str)
        if len(binning) == 0:
            sys.exit(1) # bail out, there is nothing we can do here
        costh_bins = zip(binning[:-1], binning[1:])

    costh_means = get_bin_means(dfr, lambda d: d.costh_HX.abs(), costh_bins)

    bin_sel_info = {
        'costh_bins': costh_bins,
        'costh_means': costh_means,
    }
    bin_sel_info['datafile'] = datafile
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


def run_fit(model, tree, costh_bins, datavars, outfile, refit=False,
            fix_shape=False, weights=None):
    """Import data, run fits and store the results"""
    wsp = r.RooWorkspace('ws_mass_fit')
    import_data(wsp, tree, datavars, weights)
    model.define_model(wsp)
    wsp.Print()

    if fix_shape:
        # Refitting doesn't make sense if we already fixed the shapes beforehand
        refit = False
        model.fit(wsp, 'costh_integrated', weighted_fit=weights is not None)
        shape_pars = get_shape_params(wsp, 'costh_integrated', model)
        model.fix_params(wsp, [(sp, None) for sp in shape_pars])

    do_binned_fits(model, wsp, costh_bins, refit, weights is not None)
    wsp.writeToFile(outfile)


def run_chic_fit(args):
    """Setup everything and run the chic fits"""
    logging.warn('The \'chic\' mode of costh_binned_massfits.py is deprecated.'
                 ' You should switch to using a configfile for the fitmodel')
    logging.info('Running chic fits')
    cond_mkdir_file(args.outfile)

    model = ChicMassModel('chicMass')
    costh_binning = create_bin_info_json('chic', args.binning, args.datafile,
                                         args.outfile, args.bin_info)

    dvars = get_ws_vars('chic')
    dataf = r.TFile.Open(args.datafile)
    tree = dataf.Get('chic_tuple')
    run_fit(model, tree, costh_binning, dvars, args.outfile, args.refit,
            args.fix_shape)


def run_chib_fit(args):
    """Setup everything and run the chic fits"""
    logging.warn('The \'chib\' mode of costh_binned_massfits.py is deprecated.'
                 ' You should switch to using a configfile for the fitmodel')
    logging.info('Running chib fits')
    cond_mkdir_file(args.outfile)
    model = ChibMassModel(args.configfile)
    costh_binning = create_bin_info_json('chib', args.binning, args.datafile,
                                         args.outfile, args.bin_info)

    dvars = get_ws_vars('chib', model)
    dataf = r.TFile.Open(args.datafile)
    tree = dataf.Get('data')
    run_fit(model, tree, costh_binning, dvars, args.outfile, args.refit,
            args.fix_shape)


def run_jpsi_fit(args):
    """Setup everything and run the chic fits"""
    logging.warn('The \'jpsi\' mode of costh_binned_massfits.py is deprecated.'
                 ' You should switch to using a configfile for the fitmodel')
    logging.info('Running jpsi fits')
    cond_mkdir_file(args.outfile)
    model = JpsiMassModel('JpsiMass')
    costh_binning = rw_bin_sel_json(args.costh_bin_file, args.datafile,
                                    args.bin_info)

    dvars = get_ws_vars('jpsi')
    dataf = r.TFile.Open(args.datafile)
    tree = dataf.Get('jpsi_tuple')
    run_fit(model, tree, costh_binning, dvars, args.outfile, args.refit,
            args.fix_shape)


def run_config_fit(args):
    """Setup everything and run the config file fits
    NOTE: currently assuming chic or chib input for this
    """
    logging.info('Running fits with model from config file')
    cond_mkdir_file(args.outfile)
    try:
        shutil.copy2(args.configfile,
                     os.path.join(os.path.split(args.outfile)[0],
                                  'fit_model.json'))
    except shutil.Error as err:
        logging.warn(str(err))

    model = ConfigFitModel(args.configfile)
    costh_binning = create_bin_info_json('chic', args.binning, args.datafile,
                                         args.outfile, args.bin_info)

    dvars = get_ws_vars('chic', weights=args.weight)
    dataf = r.TFile.Open(args.datafile)
    tree = dataf.Get('chic_tuple')
    run_fit(model, tree, costh_binning, dvars, args.outfile, args.refit,
            args.fix_shape, args.weight)


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
    global_parser.add_argument('--bin-info', help='json file containing '
                               'the costh binning information (output). If not '
                               'provided a default name will be used',
                               default=None)
    global_parser.add_argument('--refit', help='Fit the shape parameters and '
                               'the number of events and then rerun the fit '
                               'after fixing the shape parameters',
                               default=False, action='store_true')
    global_parser.add_argument('-s', '--fix-shape', help='fix the shape by doing'
                               ' a costh integrated fit first and fixing all '
                               'shape parameters to the ones obtained there',
                               default=False, action='store_true')
    global_parser.add_argument('-b', '--binning', help='Determine how the '
                               'variable should be binned. Either provide '
                               'a (comma separated) list of values that should '
                               ' be used as bin edges, a valid format for a '
                               'linear binning: (\'BEGIN:END:DELTA\' or '
                               '\'BEGIN:END,NSTEPS\'), or \'auto:N\' (default '
                               'N=4) which will bin the data into N equally '
                               'populated bins', default='auto:4')
    global_parser.add_argument('-w', '--weight', default=None, help='Use the '
                               'weight with the passed name present in the input'
                               ' data as per-event weight in the fit')

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
    chib_parser.set_defaults(func=run_chib_fit)

    # Add the jpsi parser
    jpsi_parser = subparsers.add_parser('jpsi', description='Run the fits using'
                                     ' jpsi mass model',
                                     parents=[global_parser])
    jpsi_parser.add_argument('costh_bin_file', help='json file containing the '
                             'costh binning information')
    jpsi_parser.set_defaults(func=run_jpsi_fit)

    # Add the config parser
    config_parser = subparsers.add_parser('config', description='Run the fits '
                                          'using a model specified in a json '
                                          'file', parents=[global_parser])
    config_parser.add_argument('configfile', help='Config file in json format')
    config_parser.set_defaults(func=run_config_fit)

    clargs = parser.parse_args()
    clargs.func(clargs)
