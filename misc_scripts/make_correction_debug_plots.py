#!/usr/bin/env python
"""
Script to make some plots from the debug information that is stored when
correcting the ratios
"""


import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, get_y_max
from utils.misc_helpers import cond_mkdir
from utils.setup_plot_style import set_basic_style
from utils.data_handling import list_obj


def determine_bin_variable(rfile):
    """
    Determine whether the ratio has been obtained as a function of costh or phi
    from the names of the histograms in the file

    NOTE: very brittle, but works for now
    """
    all_hists = list_obj(rfile, 'TH1D')
    if all('costh_HX_fold' in h for h in all_hists):
        return 'costh_HX_fold'

    if all('phi_HX_fold' in h for h in all_hists):
        return 'phi_HX_fold'


def get_xlab(state):
    """
    Get the state probability label
    """
    return 'P({})'.format({'chi1': '#chi_{c1}', 'chi2': '#chi_{c2}'}[state])


def make_rej_evs_w_state_prob_plot(rfile, state, var):
    """
    Make a plot showing the state probability of the rejected events in all the
    bins
    """
    hist_base = r'{}_JpsiPt_0_{}_([0-9]+)_state_prob_rej_ev'.format(state, var)
    hists = [rfile.Get(n) for n in list_obj(rfile, 'TH1D', hist_base)]

    can = mkplot(hists, drawOpt='PEX0', xRange=[0, 1.25], xLabel=get_xlab(state),
                 legPos=[0.83, 0.2, 0.94, 0.2 + 0.04 * len(hists)],
                 legEntries=['bin {}'.format(i) for i in xrange(len(hists))],
                 yLabel='rejected events', logy=True, yRange=[0.1, None])
    mkplot(r.TLine(1, 0, 1, get_y_max(hists) * 1.1), can=can, drawOpt='same',
           attr=[{'color': 12, 'line': 7, 'width': 2}])

    return can


def make_state_prob_plot(rfile, state, var):
    """
    Make a plot comparing (the normalized) state probability distributions in
    all the bins
    """
    hist_base = r'^{}_JpsiPt_0_{}_([0-9]+)_state_prob$'.format(state, var)
    hists = [rfile.Get(n) for n in list_obj(rfile, 'TH1D', hist_base)]
    [h.Scale(1.0 / h.Integral()) for h in hists]

    can = mkplot(hists, drawOpt='hist', xRange=[0, 1.25], xLabel=get_xlab(state),
                 legPos=[0.83, 0.2, 0.94, 0.2 + 0.04 * len(hists)], legOpt='l',
                 legEntries=['bin {}'.format(i) for i in xrange(len(hists))],
                 yLabel='normalized distributions', logy=True, yRange=[0.0003, None])
    mkplot(r.TLine(1, 0, 1, get_y_max(hists) * 1.1), can=can, drawOpt='same',
           attr=[{'color': 12, 'line': 7, 'width': 2}])

    return can


def main(args):
    """Main"""
    set_basic_style()
    dbfile = r.TFile.Open(args.debugfile)
    var = determine_bin_variable(dbfile)

    cond_mkdir(args.outdir)

    for state in ['chi1', 'chi2']:
        can = make_rej_evs_w_state_prob_plot(dbfile, state, var)
        can.SaveAs('{}/{}_state_probs_rej_evs_{}_bins.pdf'
                   .format(args.outdir, state, var))

        can = make_state_prob_plot(dbfile, state, var)
        can.SaveAs('{}/{}_state_prob_{}_bins.pdf'.format(args.outdir, state, var))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make some plots from'
                                     ' the debug information that can be stored '
                                     'when obtaining corrected ratios')
    parser.add_argument('debugfile', help='File created by correct_ratio.py when'
                        ' when used with the \'--debug\' flag')
    parser.add_argument('-o', '--outdir', help='Output directory into which the '
                        'plots will be stored', default='.')

    clargs = parser.parse_args()
    main(clargs)
