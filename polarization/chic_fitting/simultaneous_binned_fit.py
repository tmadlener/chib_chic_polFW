#!/usr/bin/env python
"""
Script to run the simultaneous binned fit with a given input data sample
"""

import json

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

def load_data(datafile, model):
    """
    Load the data into a DataFrame, but get only the variables that are
    necessary
    """
    variables = [model.fit_var] + model.bin_vars
    return get_dataframe(datafile, columns=variables)


def create_workspace(workspace_name, datafile, model):
    """
    Create the workspace and already load the data into it
    """
    wsp = r.RooWorkspace(workspace_name)

    dset_vars = r.RooArgSet()

    for ivar, var in enumerate(model.bin_vars):
        binning = model.binning[ivar]
        try_factory(wsp, '{}[{}, {}]'.format(var, min(binning), max(binning)))
        dset_vars.add(get_var(wsp, var))


    # TODO: Replace boundaries with the appropriate values for the chic mass
    # fit once switching to that data or even better put it into the config file
    try_factory(wsp, '{}[{}, {}]'.format(model.fit_var, 3.33, 3.72))
    dset_vars.add(get_var(wsp, model.fit_var))


    data = load_data(datafile, model)
    tree = array2tree(data.to_records(index=False))

    dataset = r.RooDataSet('full_data', 'full data sample', tree, dset_vars)
    ws_import(wsp, dataset)

    wsp.Print()
    return wsp


def main(args):
    """Main"""
    with open(args.configfile, 'r') as configfile:
        config = json.load(configfile)

    model = BinnedFitModel(config)
    wsp = create_workspace('ws_mass_fit', args.datafile, model)
    model.define_model(wsp)
    savename = 'twodim'
    model.fit(wsp, savename)

    cond_mkdir_file(args.outfile)
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

    clargs = parser.parse_args()
    main(clargs)
