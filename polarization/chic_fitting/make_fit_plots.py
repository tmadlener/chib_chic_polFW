#!/usr/bin/env python
"""
Make the plots of the fit results for each costh bin
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
from os.path import dirname

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.chic_fitting import ChicMassModel
from utils.jpsi_fitting import JpsiMassModel
from utils.chib_fitting import ChibMassModel
from utils.hist_utils import combine_cuts
from utils.misc_helpers import cond_mkdir, get_bin_cut_root

from common_func import get_bin_sel_info


def make_fit_res_plots(wsp, costh_bins, base_sel, state, outdir, **kwargs):
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
        full_selection = combine_cuts([costh_cut, base_sel])

        pdfname = '/'.join([outdir, 'mass_fit_{}_costh_bin_{}.pdf'
                            .format(state, i)])

        mass_model.plot(wsp, pdfname, snapname, full_selection, **kwargs)
        mass_model.plot_fit_params(wsp, pdfname.replace('.pdf', '_fit_res.pdf'),
                                   snapname)


def main(args):
    """Main"""
    ffile = r.TFile.Open(args.fitfile)
    ws = ffile.Get('ws_mass_fit')

    bin_sel_info = get_bin_sel_info(args.pklfile, args.fitfile)

    outdir = args.outdir
    if not outdir:
        outdir = dirname(args.fitfile)
    cond_mkdir(outdir)

    make_fit_res_plots(ws, bin_sel_info['costh_bins'],
                       bin_sel_info['basic_sel'], args.state, outdir,
                       logy=args.logy, configfile=args.configfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to produce the fit '
                                     'of the chic mass')
    parser.add_argument('fitfile', help='file containing the workspace with '
                        'the fit results and the dataset used for fitting')
    parser.add_argument('-pf', '--pklfile', help='Pickle file containing the '
                        'costh binning and selection informations. Use this to '
                        'override the default which derives the name from the '
                        'fitfile', default='', type=str)
    parser.add_argument('-o', '--outdir', help='Directory to which the plots '
                        'get stored (defaults to same directory as fitfile)',
                        default='', type=str)
    parser.add_argument('--logy', default=False, action='store_true',
                        help='Use log-scale on y-axis')

    state_sel = parser.add_mutually_exclusive_group()
    state_sel.add_argument('--chic', action='store_const', dest='state',
                           const='chic', help='Do mass fits for chic data')
    state_sel.add_argument('--jpsi', action='store_const', dest='state',
                           const='jpsi', help='Do mass fits for jpsi data')
    state_sel.add_argument('--chib', action='store_const', dest='state',
                           const='chib', help='Do mass fits for chib data')
    parser.set_defaults(state='chic')

    parser.add_argument('--configfile', help='Config file in json format for chib model.', type=str)


    clargs = parser.parse_args()
    main(clargs)
