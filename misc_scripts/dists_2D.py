#!/usr/bin/env python
"""
Script to generate 2D distribution plots and store them into a root file
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from collections import OrderedDict

from utils.data_handling import apply_selections, get_dataframe
from utils.hist_utils import hist2d
from utils.misc_helpers import _get_var, get_bin_cut_df

import utils.selection_functions as sf

# variables that will always be loaded
LOAD_VARIABLES = [
    'JpsiPt', 'JpsiRap',
    'costh_HX', 'phi_HX',
    'photonPt', 'photonEta',
    'muNPt', 'muNEta', 'muPPt', 'muPEta',
    'lepP_eff_sm', 'lepN_eff_sm', 'gamma_eff_sm',
]

# some helper functions
GAMMA_SEL = lambda d: sf.photon_sel(d, sf.flat_pt(0.41, 1.5))
GAMMA_EFF = lambda d: d.gamma_eff_sm * 0.01 * (d.gamma_eff_sm > 0)
MUON_EFF = lambda d: d.lepP_eff_sm * d.lepN_eff_sm
FULL_EFF = lambda d: MUON_EFF(d) * GAMMA_EFF(d)

# some predefined selections
SELECTIONS = OrderedDict()
SELECTIONS['no_sel'] = (None, None)
SELECTIONS['jpsi'] = (sf.jpsi_kin_sel, None)
SELECTIONS['jpsi_muon_sel'] = ((sf.jpsi_kin_sel, sf.loose_muon_sel), None)
SELECTIONS['jpsi_gamma_sel'] = ((sf.jpsi_kin_sel, GAMMA_SEL), None)
SELECTIONS['jpsi_gamma_muon_sel'] = ((sf.jpsi_kin_sel, sf.loose_muon_sel,
                                      GAMMA_SEL), None)
SELECTIONS['jpsi_muon_eff'] = ((sf.jpsi_kin_sel, sf.loose_muon_sel), MUON_EFF)
SELECTIONS['jpsi_gamma_eff'] = ((sf.jpsi_kin_sel, GAMMA_SEL), GAMMA_EFF)
SELECTIONS['jpsi_gamma_muon_eff'] = ((sf.jpsi_kin_sel, sf.loose_muon_sel,
                                      GAMMA_SEL), FULL_EFF)

# variables that need some functions to be obtained
# NOTE: variables used here need to be present in the LOAD_VARIABLES
FUNCVARS = {
    'lowMuPt': lambda d: d.muPPt * (d.muPPt < d.muNPt) + d.muNPt * (d.muNPt < d.muPPt),
    'highMuPt': lambda d: d.muPPt * (d.muPPt > d.muNPt) + d.muNPt * (d.muNPt > d.muPPt)
}


def get_df_and_weights(dfr, selection=None, weight=None):
    """
    Get the dataframe after selections and the according weights
    """
    sdfr = apply_selections(dfr, selection)
    return sdfr, _get_var(sdfr, weight)


def make_2D_hist(dfr, varx, vary, selection=None, weight=None, **kwargs):
    """
    Make the 2D histogram and return it
    """
    logging.debug('Making histograms for {}:{}'.format(varx, vary))
    dfr, weights = get_df_and_weights(dfr, selection, weight)
    return hist2d(_get_var(dfr, varx), _get_var(dfr, vary), weights=weights,
                  **kwargs)


def store_hists(rfile, hists, basename):
    """
    Store the histograms into the passed rootfile
    """
    rfile.cd()
    for hname, hist in hists.iteritems():
        hist.SetName('_'.join([basename, hname]))
        hist.Write()


def create_hists(dfr, hist_func, selection_weights):
    """
    Create the histograms
    """
    hists = OrderedDict()
    for name in selection_weights:
        sel_f, sel_w = selection_weights[name]
        hists[name] = hist_func(dfr, sel_f, sel_w)

    return hists


def create_costh_binned_hists(dfr, hist_func, selection_weights, binning):
    """
    Create the histograms in costh bins
    """
    hists = OrderedDict()
    for name in selection_weights:
        sel_f, sel_w = selection_weights[name]
        hists[name] = OrderedDict()

        # Apply the selections as soon as possible and don't repeat them
        sel_dfr = apply_selections(dfr, sel_f)

        for bin_low, bin_high in zip(binning[:-1], binning[1:]):
            binname = '{:.2f}_{:.2f}'.format(bin_low, bin_high)
            bin_f = lambda d: get_bin_cut_df(d, lambda d: d.costh_HX.abs(),
                                             bin_low, bin_high),
            hists[name][binname] = hist_func(sel_dfr, bin_f, sel_w)

    # Rearrange histograms so that the dictionary is flat
    flat_hists = OrderedDict()
    for sel in hists:
        for cbin in hists[sel]:
            flat_hists['_'.join([sel, cbin])] = hists[sel][cbin]
    return flat_hists


def main(args):
    """Main"""
    varx, vary = args.variablex, args.variabley
    if varx not in LOAD_VARIABLES and varx not in FUNCVARS:
        LOAD_VARIABLES.append(args.variablex)
    if vary not in LOAD_VARIABLES and vary not in FUNCVARS:
        LOAD_VARIABLES.append(args.variabley)


    selections = OrderedDict()
    for sel in args.selections.split(','):
        if sel not in SELECTIONS:
            logging.warn('Cannot produce plots for selection \'{}\', '
                         'because it is not defined'.format(sel))
        else:
            selections[sel] = SELECTIONS[sel]


    def histf(data, sel, weight):
        """Wrap the function to conform to interface"""
        vrx = varx
        vry = vary
        if varx in FUNCVARS:
            vrx = FUNCVARS[varx]
        if vary in FUNCVARS:
            vry = FUNCVARS[vary]

        return make_2D_hist(data, vrx, vry, selection=sel, weight=weight,
                            nbinsx=args.nbinsx, nbinsy=args.nbinsy,
                            minx=args.minx, maxx=args.maxx,
                            miny=args.miny, maxy=args.maxy)


    data = get_dataframe(args.inputfile, columns=LOAD_VARIABLES)
    outfile = r.TFile(args.outfile, 'recreate' if args.recreate else 'update')

    if args.costh_int:
        hists = create_hists(data, histf, selections)
    else:
        hists = create_costh_binned_hists(data, histf, selections,
                                          np.linspace(0, 1, 11))

    store_hists(outfile, hists, args.basename)
    outfile.Write('', r.TObject.kWriteDelete)
    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to create 2D '
                                     'distribution plots and store them into a '
                                     'root file.')
    parser.add_argument('inputfile', help='Input file for which the plots should'
                        ' be created')
    parser.add_argument('outfile', help='Output file containing the created '
                        'plots.')
    parser.add_argument('-r', '--recreate', help='Recreate the output file '
                        'instead of updating it', action='store_true',
                        default=False)
    parser.add_argument('-n', '--basename', help='Base name used to store the '
                        '2D distributions', default='dist2D')
    parser.add_argument('-s', '--selections', help='comma separated list of '
                        'selections for which histograms should be created. '
                        'Available: {}'.format(SELECTIONS.keys()),
                        default='jpsi')
    parser.add_argument('-vx', '--variablex', help='x-variable to plot. Has to '
                        'be available in the input file as a branch or one of '
                        'pre-defined variables: {}'.format(FUNCVARS.keys()),
                        default='lowMuPt')
    parser.add_argument('-vy', '--variabley', help='y-variable to plot. See '
                        'x-variable for more information',
                        default='highMuPt')
    parser.add_argument('-nx', '--nbinsx', help='Number of bins in x direction',
                        default=80, type=int)
    parser.add_argument('-ny', '--nbinsy', help='Number of bins in y direction',
                        default=80, type=int)
    parser.add_argument('--maxy', help='Max value in y direction', default=20,
                        type=float)
    parser.add_argument('--miny', help='Max value in y direction', default=0,
                        type=float)
    parser.add_argument('--maxx', help='Max value in x direction', default=10,
                        type=float)
    parser.add_argument('--minx', help='Max value in x direction', default=0,
                        type=float)
    parser.add_argument('--costh-int', help='Do the plots integrated in costh',
                        action='store_true', default=False)

    clargs = parser.parse_args()
    main(clargs)
