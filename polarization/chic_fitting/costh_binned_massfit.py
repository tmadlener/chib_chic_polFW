#!/usr/bin/env python
"""
Script for running costh binned mass fits
"""

import re
import json
import os
import sys

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import shutil

from root_numpy import array2tree

from utils.data_handling import get_dataframe, apply_selections
from utils.misc_helpers import (
    get_costh_binning, get_bin_means, cond_mkdir_file, get_bin_cut_root,
    parse_binning, parse_func_var, select_bin
)
from utils.roofit_utils import (
    get_var, ws_import, get_args, release_params, fix_params, try_factory
)
from utils.config_fitting import ConfigFitModel


def get_costh_bins(bin_str, binvar, data):
    """
    Get the bin edges from the binning_str and the data
    """
    auto_match = re.match(r'auto:(\d+)$', bin_str)
    if auto_match:
        costh_bins = get_costh_binning(data, int(auto_match.group(1)))
    else:
        binning = parse_binning(bin_str)
        if len(binning) == 0:
            sys.exit(1) # bail out, there is nothing we can do here
        costh_bins = zip(binning[:-1], binning[1:])

    costh_means = get_bin_means(data, binvar, costh_bins)

    return costh_bins, costh_means


def create_bin_info_json(filename, costh_bins, costh_means, bin_var, datafile):
    """
    Write the bin info json file to store the costh bins, mans and the binning
    variable
    """
    bin_sel_info = {
        'costh_bins': costh_bins,
        'costh_means': costh_means,
        'bin_variable': bin_var,
        'datafile': datafile
    }

    with open(filename, 'w') as info_file:
        json.dump(bin_sel_info, info_file, sort_keys=True, indent=2)


def create_workspace(model, datafile, binvar, binning, massrange, fitfile, weights=None):
    """
    Create the workspace with the data already imported and the model defined,
    also in charge of writing the bin info json file
    """
    wsp = r.RooWorkspace('ws_mass_fit')

    massrange = [float(v) for v in massrange.split(',')]

    # load the data and apply the mass selection of the fitting range immediately
    bin_var = parse_func_var(binvar) # necessary for loading
    variables = [model.mname, bin_var[0]]
    if weights is not None:
        variables.append(weights)
    data = apply_selections(get_dataframe(datafile, columns=variables),
                            select_bin(model.mname, *massrange))

    costh_bins, costh_means = get_costh_bins(binning, bin_var, data)
    create_bin_info_json(fitfile.replace('.root', '_bin_sel_info.json'),
                         costh_bins, costh_means, bin_var[0], datafile)


    # Create the variables in the workspace
    try_factory(wsp, '{}[{}, {}]'.format(model.mname, *massrange))
    if 'abs' in bin_var[1].__name__:
        try_factory(wsp, '{}[{}, {}]'.format(bin_var[0], -np.max(costh_bins), np.max(costh_bins)))
    else:
        try_factory(wsp, '{}[{}, {}]'.format(bin_var[0], np.min(costh_bins), np.max(costh_bins)))
    dset_vars = r.RooArgSet(get_var(wsp, model.mname), get_var(wsp, bin_var[0]))

    tree = array2tree(data.to_records(index=False))
    if weights is not None:
        try_factory(wsp, '{}[0, 1e5]'.format(weights))
        dataset = r.RooDataSet('full_data', 'full data sample', tree, dset_vars,
                               '', weights)
    else:
        dataset = r.RooDataSet('full_data', 'full data sample', tree, dset_vars)

    ws_import(wsp, dataset)

    return wsp, costh_bins


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



def do_binned_fits(mass_model, wsp, costh_bins, refit=False, weighted_fit=False,
                   bin_var=None):
    """
    Do the fits in all costh bins
    """
    if bin_var is None:
        logging.fatal('Trying to run a fit without specifying a bin variable. '
                      'It is also possible that this was not implemented for '
                      'the model that is trying to be used. Sorry :(')
        sys.exit(1)

    for ibin, ctbin in enumerate(costh_bins):
        logging.info('Running fit for bin {0}: {2:.3f} < |{1}| < {3:.3f}'
                     .format(ibin, bin_var, *ctbin))
        bin_sel = get_bin_cut_root('TMath::Abs({})'.format(bin_var), *ctbin)
        savename = 'costh_bin_{}'.format(ibin)
        mass_model.fit(wsp, savename, bin_sel, weighted_fit)
        if refit:
            # fix all but the event number variables to their current values
            # in the workspace and refit them leaving only the event numbers
            # free (and don't forget to release the shape again afterwards)
            shape_par = get_shape_params(wsp, savename, mass_model)
            fix_params(wsp, [(sp, None) for sp in shape_par])
            refit_sn = 'refit_costh_bin_{}'.format(ibin)
            mass_model.fit(wsp, refit_sn, bin_sel, weighted_fit)
            release_params(wsp, shape_par)


def run_fit(model, datafile, outfile, bin_var, binning, massrange,
            fix_shape=False, refit=False, weights=None):
    """Import data, run fits and store the results"""
    if bin_var is None:
        logging.fatal('Trying to run a fit without specifying a bin variable. '
                      'It is also possible that this was not implemented for '
                      'the model that is trying to be used. Sorry :(')
        sys.exit(1)

    wsp, costh_bins = create_workspace(model, datafile, bin_var, binning,
                                       massrange, outfile, weights)

    model.define_model(wsp)
    wsp.Print()

    if fix_shape:
        # Refitting doesn't make sense if we already fixed the shapes beforehand
        refit = False
        model.fit(wsp, 'costh_integrated', weighted_fit=weights is not None)
        shape_pars = get_shape_params(wsp, 'costh_integrated', model)
        fix_params(wsp, [(sp, None) for sp in shape_pars])

    bin_var = parse_func_var(bin_var)[0]

    do_binned_fits(model, wsp, costh_bins, refit, weights is not None, bin_var)
    wsp.writeToFile(outfile)


def run_chic_fit(args):
    """Setup everything and run the chic fits"""
    logging.error('The \'chic\' mode of costh_binned_massfits.py is deprecated.'
                  ' You should switch to using a configfile for the fitmodel')
    sys.exit(1)


def run_chib_fit(args):
    """Setup everything and run the chic fits"""
    logging.error('The \'chib\' mode of costh_binned_massfits.py is deprecated.'
                  ' You should switch to using a configfile for the fitmodel')
    sys.exit(1)


def run_jpsi_fit(args):
    """Setup everything and run the chic fits"""
    logging.error('The \'jpsi\' mode of costh_binned_massfits.py is deprecated.'
                  ' You should switch to using a configfile for the fitmodel')
    sys.exit(1)


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

    run_fit(model, args.datafile, args.outfile, args.binvariable, args.binning,
            args.massrange, args.fix_shape, args.refit, weights=args.weight)


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
    global_parser.add_argument('-bv', '--binvariable', help='The variable that '
                               'is used for defining the bins in which the fits '
                               'are done', default='abs(costh_HX_fold)')
    global_parser.add_argument('-m', '--massrange', help='Comma separated lower '
                               'and upper value of the mass window in which the '
                               'fits should be done', default='3.425,3.725')

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
