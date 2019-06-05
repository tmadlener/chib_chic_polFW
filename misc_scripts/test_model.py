#!/usr/bin/env python
"""
Test if a json config is a valid model
"""

import json

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import numpy as np
import pandas as pd
from root_numpy import array2tree

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.two_dim_binned_fitting import BinnedFitModel
from utils.roofit_utils import try_factory, get_var, ws_import


def main(args):
    """Main"""
    with open(args.configfile, 'r') as conff:
        config = json.load(conff)

    model = BinnedFitModel(config)
    wsp = r.RooWorkspace('test_ws')

    test_vars = {}
    dset_vars = r.RooArgSet()

    for var, bounds in model.get_load_vars():
        try_factory(wsp, '{}[{}]'.format(var, bounds))
        low, high = [float(v) for v in bounds.split(',')]
        test_vars[var] = np.random.uniform(low, high, 10000)
        dset_vars.add(get_var(wsp, var))

    data = pd.DataFrame(test_vars)
    tree = array2tree(data.to_records(index=False))

    dataset = r.RooDataSet('full_data', 'test data', tree, dset_vars)
    ws_import (wsp, dataset)

    if model.define_model(wsp):
        print('configfile is a valid model definition')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to test if a model json'
                                     ' is valid')
    parser.add_argument('configfile', help='The configfile to be tested')
    parser.add_argument('-v', '--verbose', help='Make the define_model print out'
                        ' all expressions to see which one fails',
                        action='store_true', default=False)

    clargs = parser.parse_args()

    if clargs.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    main(clargs)
