#!/usr/bin/env python
"""
Script to plot the simultaneous binned fit results
"""

import json

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from os.path import dirname
from collections import OrderedDict

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.two_dim_binned_fitting import BinnedFitModel
from utils.misc_helpers import cond_mkdir
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info

def store_proto_pars(wsp, model, outfile, sym_uncer):
    """
    list of proto params with or without functional dependence
    """
    simvars = []
    comvars = OrderedDict()
    for title, defn in model.proto_params:
        if isinstance(defn, (tuple, list)):
            v_dep = defn[1].split(', ')
            if not any([par_dep in v_dep for par_dep, _ in model.proto_params]):
                comvars[title] = defn
        else:
            simvars.append(title)

    for par in model.nevent_yields + ['r_chic2_chic1']:
        if par not in simvars:
            simvars.append(par)

    outf = r.TFile.Open(outfile, 'recreate')

    sim_graphs = model.plot_simvar_graphs(wsp, simvars)
    for graph in sim_graphs:
        graph.Write('',r.TObject.kWriteDelete)

    if sym_uncer:
        for graph in model.plot_simvar_graphs(wsp, simvars, sym_uncer):
            graph.Write('', r.TObject.kWriteDelete)

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


    if args.outdir is None:
        outdir = dirname(args.fitfile)
    else:
        outdir = args.outdir

    cond_mkdir(outdir)

    if args.publication:
        set_TDR_style()

    # Saving the plots
    if not args.no_plots:
        cans = model.plot(wsp, verbose=args.verbose, publication=args.publication,
                          preliminary=args.preliminary)
        canp = model.plot_fit_params(wsp)

        for bin_name in model.bins:
            plotname = '/'.join([outdir, bin_name+'_massfit.pdf'])
            if args.publication:
                add_auxiliary_info(cans[bin_name], 2012, prelim=True)
            cans[bin_name].SaveAs(plotname)

            parname = '/'.join([outdir, bin_name+'_massfit_res.pdf'])
            canp[bin_name].SaveAs(parname)

    if args.graphs:
        outfile = '/'.join([outdir, 'proto_param_graphs.root'])
        store_proto_pars(wsp, model, outfile, args.symmetric)


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
    parser.add_argument('-o', '--outdir', help='Output directory for the created'
                        ' plot files (defaults to the directory of the fitfile)',
                        default=None)
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Add some debug information to the produced plots')
    parser.add_argument('--symmetric', help='Also create graphs with symmetric '
                        'uncertainties in addition to the ones with asymmetric '
                        'uncertainties', action='store_true', default=False)
    parser.add_argument('--no-plots', help='Do not produce the fit result pdf '
                        'plots', action='store_true', default=False)
    parser.add_argument('--publication', action='store_true', default=False,
                        help='Make the mass fit plots such that they are fit for'
                        ' publication')
    parser.add_argument('--preliminary', action='store_true', default=False,
                        help='Add \'preliminary\' label to publication plots')


    clargs = parser.parse_args()
    main(clargs)
