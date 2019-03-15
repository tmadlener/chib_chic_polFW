#!/usr/bin/env python
"""
Script to plot the simultaneous binned fit results
"""

import json

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from os.path import dirname

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
    canp = model.plot_fit_params(wsp)

    # TODO: conditionally making the outdir
    outdir = dirname(args.fitfile)

    # Saving the plots
    for bin_name, bin_borders in model.bins.iteritems():
        plotname = '/'.join([outdir, bin_name+'_massfit.pdf'])
        cans[bin_name].SaveAs(plotname)
        parname = '/'.join([outdir, bin_name+'_massfit_res.pdf'])
        canp[bin_name].SaveAs(parname)
    canp['full_model'].SaveAs('/'.join([outdir, 'full_massfit_res.pdf']))
    

    # TODO: save the plots within the utils file without returning the canvases
    # TODO: fix whatever is causing the canvas to appear on screen while running the code (not sure if above)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the results of '
                                     'the simultaneous binned fit')
    parser.add_argument('fitfile', help='File containing the fit results')
    parser.add_argument('configfile', help='json file containing the '
                        'configuration of the fit model and the data binning')


    clargs = parser.parse_args()
    main(clargs)
