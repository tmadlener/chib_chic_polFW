#!/usr/bin/env python
"""
Make the plots of the fit results for each costh bin
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from os.path import dirname

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.chic_fitting import ChicMassModel
from utils.jpsi_fitting import JpsiMassModel
from utils.chib_fitting import ChibMassModel
from utils.hist_utils import combine_cuts
from utils.misc_helpers import cond_mkdir, get_bin_cut_root

from common_func import get_bin_sel_info

def make_fit_res_plots(wsp, costh_bins, state, outdir, **kwargs):
    """
    Make the plots with the fit results

    kwargs forwarded to FitModel.plot
    """
    if state == 'chic':
        mass_model = ChicMassModel('chicMass')
    elif state == 'chib':
        mass_model = ChibMassModel(kwargs.pop("configfile"))
    else:
        mass_model = JpsiMassModel('JpsiMass')

    for i, ctbin in enumerate(costh_bins):
        costh_cut = get_bin_cut_root('TMath::Abs(costh_HX)', *ctbin)
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


def main(args):
    """Main"""
    ffile = r.TFile.Open(args.fitfile)
    ws = ffile.Get('ws_mass_fit')

    bin_sel_info = get_bin_sel_info(args.bin_info, args.fitfile)

    outdir = args.outdir
    if not outdir:
        outdir = dirname(args.fitfile)
    cond_mkdir(outdir)

    make_fit_res_plots(ws, bin_sel_info['costh_bins'],
                       args.state, outdir, logy=args.logy,
                       configfile=args.configfile, ppars=args.print_pars,
                       corr_matrix=args.corr_matrix, refit=args.refit)


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

    state_sel = parser.add_mutually_exclusive_group()
    state_sel.add_argument('--chic', action='store_const', dest='state',
                           const='chic', help='Do mass fits for chic data')
    state_sel.add_argument('--jpsi', action='store_const', dest='state',
                           const='jpsi', help='Do mass fits for jpsi data')
    state_sel.add_argument('--chib', action='store_const', dest='state',
                           const='chib', help='Do mass fits for chib data')
    parser.set_defaults(state='chic')

    parser.add_argument('--configfile', help='Config file in json format for chib model.',
                        default="config.json", type=str)


    clargs = parser.parse_args()
    main(clargs)
