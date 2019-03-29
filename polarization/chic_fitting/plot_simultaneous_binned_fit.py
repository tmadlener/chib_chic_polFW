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
from utils.graph_utils import assign_x


def plot_proto_pars(wsp, model, outfile):
    # auxiliary definitions
    bin_var_X = model.bin_vars[0]
    bin_var_Y = model.bin_vars[1]

    #auxiliary dict for getting the parameters in the right order of bins
    p_name = OrderedDict()
    p_name[bin_var_X] = ('j','i')
    p_name[bin_var_Y] = ('i','j')

    #list of proto params with or without functional dependence
    simvars = []
    comvars = OrderedDict()
    for title, defn in model.proto_params.iteritems():
        if isinstance(defn, (tuple, list)):
            comvars[title] = defn
        else:
            simvars.append(title)

    outf = r.TFile.Open(outfile, 'recreate')

    #cycle over all simple proto-parameters in each binning var that saves the TGraphs
    for bin_var in model.bin_vars:
        bintovar = model.bintovar[bin_var]
        #get parameter in the right binning order
        p_getname = '{param}_{x_var}_{'+p_name[bin_var][0]+'}_{y_var}_{'+p_name[bin_var][1]+'}'
        b_getname = p_getname[8:]

        for param in simvars:
            vals = []

            #get the values of the parameter for each bin of X,Y
            for i in xrange(len(model.binning[1-bintovar]) - 1):
                vals.append([get_var(wsp, p_getname.format(param=param, x_var=bin_var_X, y_var=bin_var_Y, i=i, j=j)) for j in xrange(len(model.binning[bintovar]) - 1)])

            #plot graphs as a function of binning var
            graph = model.plot_free_pars(wsp, bin_var, vals)

            #setting the correct mean (errors adjust accordingly), setting graph name
            for i in xrange(len(model.binning[1-bintovar]) - 1):
                g_mean = np.array([model.bin_mean(wsp, bin_var, b_getname.format(x_var=bin_var_X, y_var=bin_var_Y, i=i, j=j)) for j in xrange(len(model.binning[bintovar]) - 1)])
                graph[i] = assign_x(graph[i], g_mean)
                # graph name of the form [y_axis]_v_[x_axis]_bin_[bin of other bin_var]
                graph[i].SetName('{}_v_{}_bin_{}'.format(param, bin_var, i))
                graph[i].Write('',r.TObject.kWriteDelete)

    for el in comvars:
        # func stores the passed TF1, as well as the binning variable on which the function depends
        func = model.tf1_helper(wsp, comvars[el][0], comvars[el][1])
        func[0].SetName('{}_v_{}'.format(el, func[1]))
        func[0].Write('',r.TObject.kWriteDelete)

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
        plot_proto_pars(wsp, model, outfile)


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
