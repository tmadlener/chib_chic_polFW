#!/usr/bin/env python
"""
Script to plot the simultaneous binned fit results
"""

import json

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.two_dim_binned_fitting import BinnedFitModel

def main(args):
    """Main"""
    with open(args.configfile, 'r') as configfile:
        config = json.load(configfile)

    model = BinnedFitModel(config)
    ffile = r.TFile.Open(args.fitfile)
    wsp = ffile.Get('ws_mass_fit')

    cans = model.plot(wsp)

    # TODO: Saving the plots


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the results of '
                                     'the simultaneous binned fit')
    parser.add_argument('fitfile', help='File containing the fit results')
    parser.add_argument('configfile', help='json file containing the '
                        'configuration of the fit model and the data binning')


    clargs = parser.parse_args()
    main(clargs)
