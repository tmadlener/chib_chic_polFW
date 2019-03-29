#!/usr/bin/env python
"""
Script to plot the simultaneous binned fit results
"""

import json

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
import numpy as np

from os.path import dirname
from collections import OrderedDict

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.two_dim_binned_fitting import BinnedFitModel
from utils.roofit_utils import get_var


def store_proto_pars(wsp, model, outfile):
    #list of proto params with or without functional dependence
    simvars = []
    comvars = OrderedDict()
    for title, defn in model.proto_params.iteritems():
        if isinstance(defn, (tuple, list)):
            comvars[title] = defn
        else:
            simvars.append(title)

    outf = r.TFile.Open(outfile, 'recreate')

    sim_graphs = model.plot_simvar_graphs(wsp, simvars)
    for graph in sim_graphs:
        graph.Write('',r.TObject.kWriteDelete)

    com_funcs = model.plot_comvar_funcs(wsp, comvars)
    for func in com_funcs:
        func.Write('',r.TObject.kWriteDelete)

    outf.Write('', r.TObject.kWriteDelete)
    outf.Close()


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

    if args.graphs:
        outfile = '/'.join([outdir, 'proto_param_graphs.root'])
        store_proto_pars(wsp, model, outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the results of '
                                     'the simultaneous binned fit')
    parser.add_argument('fitfile', help='File containing the fit results')
    parser.add_argument('configfile', help='json file containing the '
                        'configuration of the fit model and the data binning')

    parser.add_argument('-g', '--graphs', help='Create graphs of the proto '
                        'parameters vs each bin var and store them into a root '
                        'file in the output directory', action='store_true',
                        default=False)


    clargs = parser.parse_args()
    main(clargs)
