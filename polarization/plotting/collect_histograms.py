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
from utils.misc_helpers import get_full_trigger

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
weight_branches = {'chic1': 'wChic1', 'chic2': 'wChic2'}



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


def get_proto_hist(var, name):
    """
    Get prototype histogram for a given variable

    TODO: doc
    """
    ## Basic default settings, which will be used
    ## NOTE: In the future it is planned to make it possible to override these
    ## via a Jason file
    logging.debug('Getting prototype histogram for var: {}'.format(var))

    # there were root versions where this lead to automatic binning -> default
    hist = r.TH1D('', '', 100, 1, -1)

    hist_var = get_key_from_var(var)
    if hist_var:
        histset = default_hist_set[hist_var]
        logging.debug('Using histogram settings {}'.format(histset))
        hist =  r.TH1D(name, '',
                       histset['n_bins'], histset['min'], histset['max'])
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
    ratio = hnum.Clone(hnum.GetName().replace('chic2', 'ratio'))
    ratio.Divide(hdenom)
    return ratio


def get_hists_for_frame(tree, frame, pvars, cutsigma=None):
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
        for state in ['chic1', 'chic2']:
            hists[var][state] = get_proto_hist(pvar, '_'.join([state, pvar]))
            draw_expr = get_draw_expr(var, frame)
            weight = weight_branches[state]
            draw_var_to_hist(tree, hists[var][state], draw_expr, '', weight)

        hists[var]['ratio'] = divide_hists(hists[var]['chic2'],
                                           hists[var]['chic1'],
                                           cutsigma)
    return hists


def get_hists_from_file(rfile, treename, frames, cutsigma=None):
    """
    Get all histograms from file

    TODO: doc
    """
    tree = rfile.Get(treename)
    hists = {}
    for frame in frames:
        hists[frame] = get_hists_for_frame(tree, frame,
                                           ['phi', 'costh', 'cosalpha'],
                                           cutsigma)

    return hists


def save_hists_to_file(hists, filen, subdir):
    """
    Save histograms to file (using unique name)
    """
    logging.info('Saving histograms to \'{}\''.format(filen))
    outfile = r.TFile(filen, 'update')
    outfile.mkdir(subdir)
    outfile.cd(subdir)

    for frame in hists:
        for var in hists[frame]:
            for state in hists[frame][var]:
                hists[frame][var][state].Write()

    outfile.Write()
    outfile.Close()


def get_unique_subdir(year, trigger, mc=False, pt=''):
    """
    Create a subdir that uniquely identifies the conditions that were used to
    get the histograms

    TODO: doc
    """
    data_mc = 'mc' if mc else 'data'
    subdirs = [p for p in [data_mc, year, trigger, pt] if p]
    return '/'.join(subdirs)


def main(args):
    """Main"""
    logging.info('Processing file: {}'.format(args.inputfile))
    infile = r.TFile.Open(args.inputfile)
    frames = args.frames.split(',')

    file_hists = get_hists_from_file(infile, args.treename, frames,
                                     args.cutsigma)

    trigger = get_full_trigger(args.triggerpath)
    subdir = get_unique_subdir(args.year, trigger, args.mc, args.ptbin)

    save_hists_to_file(file_hists, args.outfile, subdir)


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
    parser.add_argument('-pt)', '--ptbin', default='0',
                        help='pT bin of this plot')
    parser.add_argument('-c', '--cutsigma', default=None, type=float,
                        help='Apply the passed cut on the distribution hists '
                        'to get rid of bins that are compatible to zero within '
                        'n sigmas. This is done before doing the ratio')

    clargs = parser.parse_args()
    main(clargs)
