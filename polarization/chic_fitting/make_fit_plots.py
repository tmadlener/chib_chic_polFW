#!/usr/bin/env python
"""
Make the plots of the fit results for each costh bin
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
import numpy as np

from os.path import dirname

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.chic_fitting import ChicMassModel
from utils.jpsi_fitting import JpsiMassModel
from utils.chib_fitting import ChibMassModel
from utils.config_fitting import ConfigFitModel
from utils.misc_helpers import cond_mkdir, get_bin_cut_root, get_bin_edges
from utils.roofit_utils import get_var_graph, get_args, get_var
from utils.graph_utils import has_sym_uncer

from common_func import get_bin_sel_info


def make_fit_res_plots(wsp, costh_bins, state, outdir, mass_model, bin_var=None, **kwargs):
    """
    Make the plots with the fit results

    kwargs forwarded to FitModel.plot
    """
    if bin_var is None:
        bin_var = 'costh_HX_fold' # backwards compatibility

    if kwargs.pop('fix_shape', False):
        snapname = 'snap_costh_integrated'
        pdfname = '/'.join([outdir, 'mass_fit_{}_costh_integrated.pdf'
                            .format(state)])
        mass_model.plot(wsp, pdfname, snapname, **kwargs)
        mass_model.plot_fit_params(wsp, pdfname.replace('.pdf', '_fit_res.pdf'),
                                   snapname)
        if kwargs.get('ppars', False):
            mass_model.print_fit_params(wsp, 'fit_res_costh_integrated')
        if kwargs.get('corr_matrix', False):
            mass_model.plot_corr_matrix(wsp, 'fit_res_costh_integrated',
                                        pdfname.replace('.pdf', '_corrmat.pdf'))

    for i, ctbin in enumerate(costh_bins):
        costh_cut = get_bin_cut_root('TMath::Abs({})'.format(bin_var), *ctbin)
        # costh_cut = get_bin_cut_root('TMath::Abs(costh_HX)', *ctbin)
        snapname = 'snap_costh_bin_{}'.format(i)

        pdfname = '/'.join([outdir, 'mass_fit_{}_costh_bin_{}.pdf'
                            .format(state, i)])

        mass_model.plot(wsp, pdfname, snapname, costh_cut, **kwargs)
        mass_model.plot_fit_params(wsp, pdfname.replace('.pdf', '_fit_res.pdf'),
                                   snapname)

        if kwargs.get('ppars', False):
            mass_model.print_fit_params(wsp, 'fit_res_costh_bin_{}'.format(i))

        if kwargs.get('corr_matrix', False):
            mass_model.plot_corr_matrix(wsp, 'fit_res_costh_bin_{}'.format(i),
                                        pdfname.replace('.pdf', '_corrmat.pdf'))

        if kwargs.get('refit', False):
            snapname = 'snap_refit_costh_bin_{}'.format(i)
            pdfname = '/'.join([outdir, 'mass_fit_{}_costh_bin_{}_refit.pdf'
                              .format(state, i)])

            mass_model.plot(wsp, pdfname, snapname, costh_cut, **kwargs)
            mass_model.plot_fit_params(wsp,
                                       pdfname.replace('.pdf', '_fit_res.pdf'),
                                       snapname)

            if kwargs.get('ppars', False):
                mass_model.print_fit_params(wsp, 'fit_res_refit_costh_bin_{}'
                                            .format(i))

            if kwargs.get('corr_matrix', False):
                mass_model.plot_corr_matrix(wsp, 'fit_res_refit_costh_bin_{}'
                                            .format(i),
                                            pdfname.replace('.pdf', '_corrmat.pdf'))


def get_free_params(wsp, mass_model):
    """
    Get the names of the free parameters
    """
    # NOTE: simply assuming here that all the bins have the same amount of free
    # parameters
    fit_res = wsp.genobj('fit_res_costh_bin_0')
    float_pars = [v.GetName() for v in get_args(fit_res.floatParsFinal())]
    # The yields and the ratio should always be produced
    for par in mass_model.nevent_yields + ['r_chic2_chic1']:
        if par not in float_pars:
            float_pars.append(par)

    return float_pars


def store_graphs(wsp, outfile, bin_info, mass_model):
    """
    Create and store the graphs for all the free parameters and store them
    into the (newly created) outfile
    """
    costh_binning = get_bin_edges(bin_info['costh_bins'])
    n_bins = len(costh_binning) - 1
    costh_means = np.array(bin_info['costh_means'])
    params = get_free_params(wsp, mass_model)

    outf = r.TFile.Open(outfile, 'recreate')
    for param in params:
        graph = get_var_graph(wsp, 'snap_costh_bin_{}', param, n_bins,
                              costh_binning, costh_means,
                              fit_res_base='fit_res_costh_bin_{}')
        if has_sym_uncer(graph):
            graph.SetTitle('symmetric')
        else:
            graph.SetTitle('asymmetric')

        graph.SetName('_'.join([param, 'v', 'costh']))
        graph.Write('', r.TObject.kWriteDelete)

    outf.Write('', r.TObject.kWriteDelete)
    outf.Close()


def main(args):
    """Main"""
    ffile = r.TFile.Open(args.fitfile)
    ws = ffile.Get('ws_mass_fit')

    bin_sel_info = get_bin_sel_info(args.bin_info, args.fitfile)

    outdir = args.outdir
    if not outdir:
        outdir = dirname(args.fitfile)
    cond_mkdir(outdir)

    if args.state == 'chic':
        mass_model = ChicMassModel('chicMass')
    elif args.state == 'chib':
        mass_model = ChibMassModel(args.configfile)
    elif args.state == 'config':
        mass_model = ConfigFitModel(args.configfile)
    else:
        mass_model = JpsiMassModel('JpsiMass')

    bin_var = None
    if 'bin_variable' in bin_sel_info:
        bin_var = bin_sel_info['bin_variable']

    if not args.no_fit_plots:
        make_fit_res_plots(ws, bin_sel_info['costh_bins'],
                           args.state, outdir, mass_model, bin_var,
                           logy=args.logy,
                           configfile=args.configfile, ppars=args.print_pars,
                           corr_matrix=args.corr_matrix, refit=args.refit,
                           fix_shape=args.fix_shape, weighted_fit=args.weight,
                           verbose=args.verbose)

    if args.graphs:
        outfile = '/'.join([outdir, 'free_fit_param_graphs.root'])

        store_graphs(ws, outfile, bin_sel_info, mass_model)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to produce the fit '
                                     'of the chic mass')
    parser.add_argument('fitfile', help='file containing the workspace with '
                        'the fit results and the dataset used for fitting')
    parser.add_argument('-b', '--bin-info', help='Json file containing the '
                        'costh binning. Use this to override the default which'
                        ' derives the name from the fitfile', default='', type=str)
    parser.add_argument('-o', '--outdir', help='Directory to which the plots '
                        'get stored (defaults to same directory as fitfile)',
                        default='', type=str)
    parser.add_argument('--logy', default=False, action='store_true',
                        help='Use log-scale on y-axis')
    parser.add_argument('-p', '--print-pars', help='Print the free parameters '
                        'to the terminal as well', action='store_true',
                        default=False)
    parser.add_argument('-c', '--corr-matrix', help='Make a plot of the corr '
                        'matrix of all the free parameters', default=False,
                        action='store_true')
    parser.add_argument('--refit', help='Plot the refitted results as well',
                        default=False, action='store_true')
    parser.add_argument('-s', '--fix-shape', help='Plot the results of the '
                        'costh integrated fit that was used to fix the shape '
                        'parameters', action='store_true', default=False)
    parser.add_argument('-g', '--graphs', help='Create graphs of the free fit '
                        'parameters vs costh and store them into a root file in'
                        'the output directory', action='store_true',
                        default=False)
    parser.add_argument('-w', '--weight', help='Assume that the fit has been '
                        'done using weights.', action='store_true',
                        default=False)
    parser.add_argument('-v', '--verbose', help='Put some more status info from '
                        'the fit onto the plots', action='store_true',
                        default=False)
    parser.add_argument('--no-fit-plots', help='Do not create the fit plots '
                        'overlaying the fitted function with the mass '
                        'distribution', action='store_true', default=False)

    state_sel = parser.add_mutually_exclusive_group()
    state_sel.add_argument('--chic', action='store_const', dest='state',
                           const='chic', help='Do mass fits for chic data')
    state_sel.add_argument('--jpsi', action='store_const', dest='state',
                           const='jpsi', help='Do mass fits for jpsi data')
    state_sel.add_argument('--chib', action='store_const', dest='state',
                           const='chib', help='Do mass fits for chib data')
    state_sel.add_argument('--config', help='Do mass fit plots for a model '
                           'defined via a config file',
                           const='config', action='store_const', dest='state')
    parser.set_defaults(state='chic')

    parser.add_argument('--configfile', help='Config file in json format for chib model.',
                        default="config.json", type=str)


    clargs = parser.parse_args()
    main(clargs)
