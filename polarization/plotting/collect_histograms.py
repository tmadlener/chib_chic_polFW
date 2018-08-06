#!/usr/bin/env python
"""
Create histograms and put them into a root file for later retrieval by plotting
script.

Attributes:
    default_hist_set (dict): Dictionary containing default settings for histos
    proto_draw_expr (dict): Dictionary containing the prototypes for the draw
        expressions for the different variables
    weight_branches (dict): Names of the chic1 and chic2 signal weight branches
"""

import re
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')


from utils.hist_utils import (
    draw_var_to_hist, set_hist_opts, set_bins_to_zero
)
from utils.data_handling import check_branch_available
from utils.PlotServer import PlotServer


default_hist_set = {
    'costh': {'n_bins': 16, 'min': 0, 'max': 1},
    'phi': {'n_bins': 20, 'min': 0, 'max': 90},
    'cosalpha': {'n_bins': 16, 'min': 0, 'max': 1}
}
proto_draw_expr = {
    'costh': 'TMath::Abs(costh_{})',
    'phi': 'phi_{}_fold',
    'cosalpha': 'TMath::Abs(cosalpha_{})'
}
weight_branches = {'chic1': 'wChic1', 'chic2': 'wChic2',
                   'chib1': 'N1_sw', 'chib2': 'N2_sw'}



def get_key_from_var(var):
    """
    Get the key that is used in the global variables from the variable

    TODO: doc
    """
    var_rgx = r'_?(cos(alpha|th)|phi)_[A-Z]{2}'
    var_match = re.search(var_rgx, var)
    if var_match:
        hist_var = var_match.group(1)
        return hist_var
    return ''


def get_proto_hist(var, name, nbins=None, binning=None):
    """
    Get prototype histogram for a given variable

    TODO: doc
    """
    ## Basic default settings, which will be used
    ## NOTE: In the future it is planned to make it possible to override these
    ## via a JSON file
    logging.debug('Getting prototype histogram for var: {}'.format(var))

    # there were root versions where this lead to automatic binning -> default
    hist = r.TH1D('', '', 100, 1, -1)

    hist_var = get_key_from_var(var)
    if hist_var:
        histset = default_hist_set[hist_var]
    # This way, because otherwise logging is not correct
        if nbins: histset['n_bins'] = nbins
        logging.debug('Using histogram settings {}'.format(histset))

        if binning is not None: hist =  r.TH1D(name, '', len(binning)-1, binning)
        else: hist =  r.TH1D(name, '', histset['n_bins'],
                             histset['min'], histset['max'])
        set_hist_opts(hist)
    else:
        logging.warning('Could not get histogram settings for var: {}'
                        .format(var))

    return hist


def get_draw_expr(var, frame):
    """
    Get the draw expression from the variable
    """
    return proto_draw_expr[var].format(frame)


def divide_hists(hnum, hdenom, cutsigma=None):
    """
    Divide histograms return new histogram

    Args:
        hnum (ROOT.TH1): numerator histogram
        hdenom (ROOT.TH1): denominator histogram
        cutsigma (float, optional): significance above 0, for which bins should
           be set to 0. (I.e. sqrt(N), where N is entries in bin, so for setting
           bins below 10 to 0, set cutsigma to sqrt(10))

    Returns:
        ROOT.TH1: ratio of hnum / hdenom
    """
    if cutsigma is not None:
        set_bins_to_zero(hnum, cutsigma**2)
        set_bins_to_zero(hdenom, cutsigma**2)

    # NOTE: If there is no 'chic2' in the name this will probably overwrite the
    # original
    # TODO: make slightly more versatile (and move to hist_utils probably)
    rationame = re.sub(r'chi[b|c]\d', r'ratio', hnum.GetName())
    ratio = hnum.Clone(rationame)
    ratio.Divide(hdenom)
    return ratio


def get_hists_for_frame(tree, frame, pvars, cutsigma=None, 
                        states=['chic1', 'chic2'], nbins=None,
                        binning=None):
    """
    Get all histograms for a frame

    TODO: doc
    """
    hists = {}
    for var in pvars:
        pvar = '_'.join([var, frame])
        if not check_branch_available(tree, pvar):
            continue
        hists[var] = {}
        for state in states:
            hists[var][state] = get_proto_hist(pvar, '_'.join([state, pvar]),
                                               nbins, binning)
            draw_expr = get_draw_expr(var, frame)
            weight = weight_branches[state]
            draw_var_to_hist(tree, hists[var][state], draw_expr, '', weight)

        # assuming that second element in states is numerator
        hists[var]['ratio'] = divide_hists(hists[var][states[1]],
                                           hists[var][states[0]],
                                           cutsigma)
    return hists


def get_hists_from_file(rfile, treename, frames, cutsigma=None,
                        states=['chic1', 'chic2'], nbins=None, 
                        binning=None):
    """
    Get all histograms from file

    TODO: doc
    """
    tree = rfile.Get(treename)
    hists = {}
    for frame in frames:
        hists[frame] = get_hists_for_frame(tree, frame,
                                           ['phi', 'costh', 'cosalpha'],
                                           cutsigma, states, nbins, binning)

    return hists


def save_hists_to_file(hists, filen, year, trigger, top_dir, pt):
    """
    Save the histograms to file (so that they can later be retrieved via the
    PlotServer).
    """
    logging.info('Saving histograms to \'{}\''.format(filen))
    pserver = PlotServer(filen, 'update')
    for frame in hists:
        for var in hists[frame]:
            for state in hists[frame][var]:
                pserver.store_hist(hists[frame][var][state], top_dir, year,
                                   trigger, pt, var, frame, state)


def main(args):
    """Main"""
    logging.info('Processing file: {}'.format(args.inputfile))
    infile = r.TFile.Open(args.inputfile)
    frames = args.frames.split(',')

    states = ['chic1', 'chic2']
    if args.state == 'chib':
        states = ['chib1', 'chib2']

    binning=None
    if args.binning:
        binning = np.array(args.binning.split(',')).astype(np.float)
        if binning[-1]==-1:
            from utils.data_handling import get_dataframe
            dfr = get_dataframe(args.inputfile, args.treename, columns=['costh_HX'])
            binning[-1]=dfr['costh_HX'].abs().max()
    if args.jsonbinning:
        import json
        with open(args.jsonbinning, 'r') as jsonf:
            tmpbins = json.load(jsonf)['costh_bins']
            binning = [b[0] for b in tmpbins]
            binning.append(tmpbins[-1][-1])
            binning = np.array(binning).astype(np.float)
    print(binning)
        
    file_hists = get_hists_from_file(infile, args.treename, frames,
                                     args.cutsigma, states, args.nbins,
                                     binning)

    top_dir = args.topdir
    if not top_dir:
        top_dir = 'mc' if args.mc else 'data'

    save_hists_to_file(file_hists, args.outfile,
                       args.year, args.triggerpath, top_dir, args.ptbin)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for creating the '
                                     'histograms and putting them into a root '
                                     'file for latter retrieval by plotting '
                                     'script.')
    parser.add_argument('inputfile', help='Data file from which the histograms '
                        'should be produced')
    parser.add_argument('outfile', help='Output file to which the created '
                        'histograms should be added')
    parser.add_argument('-t', '--treename', default='chic_tuple',
                        help='Name of the Tyree to use for creating the histos')
    parser.add_argument('-f', '--frames', default='HX,CS',
                        help='Comma separated list of desired frames')
    parser.add_argument('-y', '--year', default='2012',
                        help='Year of data taking (default 2012)')
    parser.add_argument('-mc', action='store_true', default=False,
                        help='MC or data sample (default data)')
    parser.add_argument('-trg', '--triggerpath', default='Dimuon8_Jpsi',
                        help='Essential parts of the trigger path')
    parser.add_argument('-pt', '--ptbin', default=0, type=int,
                        help='pT bin of this plot')
    parser.add_argument('-c', '--cutsigma', default=None, type=float,
                        help='Apply the passed cut on the distribution hists '
                        'to get rid of bins that are compatible to zero within '
                        'n sigmas. This is done before doing the ratio')
    parser.add_argument('-n', '--nbins', default=None, type=int,
                        help='Set number of bins.')
    parser.add_argument('-d', '--topdir', help='Top level directory name to be '
                        'used in the directory structure created in the output '
                        'file. If not empty, overrides \'mc\' or \'data\' which'
                        ' are the default', default='', type=str)
    parser.add_argument('--jsonbinning', help='json file containing bin information'
                       ' as produced by create_bin_info_json in costh_binned_massfit.py',
                       type=str, default='')    
    parser.add_argument('--binning', help='Use fixed binning defined here '
                        '(separated by commas, e.g. 0.1,0.2,0.4,-1), '
                        'if last bin equals -1, then the maximum costh value is taken',
                        type=str, default='')

    state_sel = parser.add_mutually_exclusive_group()
    state_sel.add_argument('--chic', action='store_const', dest='state',
                           const='chic', help='collect histograms and store '
                           'them assuming that data is chic data')
    state_sel.add_argument('--chib', action='store_const', dest='state',
                           const='chib', help='collect histograms and store '
                           'them assuming that data is chib data')
    parser.set_defaults(state='chic')

    clargs = parser.parse_args()
    main(clargs)
