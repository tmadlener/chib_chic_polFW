#!/usr/bin/env python
"""
Script to run the simultaneous binned fit with a given input data sample
"""

import json
import os
import shutil
import sys

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from root_numpy import array2tree

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.two_dim_binned_fitting import BinnedFitModel
from utils.data_handling import get_dataframe
from utils.roofit_utils import ws_import, get_var, try_factory
from utils.misc_helpers import cond_mkdir_file


def create_workspace(workspace_name, datafile, model):
    """
    Create the workspace and already load the data into it
    """
    wsp = r.RooWorkspace(workspace_name)

    dset_vars = r.RooArgSet()
    variables = []

    for var, bounds in model.get_load_vars():
        try_factory(wsp, '{}[{}]'.format(var, bounds))
        dset_vars.add(get_var(wsp, var))
        variables.append(var)

    data = get_dataframe(datafile, columns=variables)
    tree = array2tree(data.to_records(index=False))

    dataset = r.RooDataSet('full_data', 'full data sample', tree, dset_vars)
    ws_import(wsp, dataset)

    wsp.Print()
    return wsp


def main(args):
    """Main"""
    with open(args.configfile, 'r') as configfile:
        config = json.load(configfile)

    cond_mkdir_file(args.outfile)
    try:
        shutil.copy2(args.configfile,
                     os.path.join(os.path.split(args.outfile)[0],
                                  'fit_model.json'))
    except shutil.Error as err:
        logging.warn(str(err))


    model = BinnedFitModel(config)
    wsp = create_workspace('ws_mass_fit', args.datafile, model)
    if not model.define_model(wsp):
        logging.error('Cannot define model in workspace. Check outputs to see '
                      'why the definition fails (Rerun with --verbosity 1 to '
                      'have more output)')
        sys.exit(1)

    model.fit(wsp, args.verbosity)

    wsp.writeToFile(args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script that runs the '
                                     'simultaneous binned fit with a given '
                                     'input data sample and stores the fit '
                                     'results into a root file that can then be '
                                     'used for plotting')
    parser.add_argument('datafile',  help='The file containing the data sample '
                        'to be used in the fit. All selections have to be '
                        'already applied to this sample')
    parser.add_argument('configfile', help='The json file containing the '
                        'configuration of the fit model and the data binning')
    parser.add_argument('-o', '--outfile', help='The name of the file that will '
                        'be created by this script and contains the fit results',
                        default='fit_results.root')
    parser.add_argument('-v', '--verbosity', help='The verbosity of the log '
                        'outputs. (-1 for minimal output, 0 for normal output, '
                        '1 for verbose output)', default=0, type=int,
                        choices={0, 1, -1})

    clargs = parser.parse_args()

    if clargs.verbosity == 1:
        logging.getLogger().setLevel(logging.DEBUG)
    elif clargs.verbosity == -1:
        logging.getLogger().setLevel(logging.WARNING)
    else:
        logging.getLogger().setLevel(logging.INFO)

    main(clargs)
