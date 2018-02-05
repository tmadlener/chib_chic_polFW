#!/usr/bin/env python
"""
Script to check the effect of fiducial single muon cuts on chic2 / chic1
polarization ratio

"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.hist_utils import combine_cuts, set_bins_to_zero
from utils.plot_helpers import mkplot
from helpers import (
    get_plot_varname, var_selection, get_histogram, select_trigger
)
from utils.misc_helpers import cond_mkdir

def fid_cuts(ptname, etaname):
    """
    Fiducial cuts for one singlem muon

    Args:
        ptname (str): branch name of muon pt
        etaname (str): branch name of muon eta

    Returns:
        str: cut expression usable in TTree::Draw()
    """
    cuts = []
    cuts.append(combine_cuts([ptname + ' > 4.5',
                              'TMath::Abs(' + etaname + ') < 1.2']))
    cuts.append(combine_cuts([ptname + ' > 4.0',
                              var_selection('TMath::Abs('+etaname+')', 1.2, 1.4)
                              ]))
    cuts.append(combine_cuts([ptname + ' > 3.5',
                              var_selection('TMath::Abs('+etaname+')', 1.4, 1.6)
                              ]))
    return combine_cuts(cuts, ' || ')


def get_fidcuts():
    """
    Complete single muon fiducial cuts for both muons

    Returns:
        str: cut expression usable in TTree::Draw()
    """
    return combine_cuts([fid_cuts('muN_pt', 'muN_eta'),
                         fid_cuts('muP_pt', 'muP_eta')])


def get_chic_ratio(rfile, var, selection):
    """
    Get the ratio histogram of chic2 / chic1 for the passed var
    """
    h_chic1 = get_histogram(rfile, var, 'chic1', selection)
    h_chic2 = get_histogram(rfile, var, 'chic2', selection)

    ## apply same minimum bin content cut as on data
    set_bins_to_zero(h_chic1, thresh=r.TMath.Sqrt(10))
    set_bins_to_zero(h_chic2, thresh=r.TMath.Sqrt(10))

    ratio = h_chic2.Clone(h_chic2.GetName().replace('chic2', 'ratio'))
    ratio.Divide(h_chic1)

    return ratio


def make_plot(rfile, var, x_label, savename):
    """
    Make plot of a specific variable
    """
    r_nocut = get_chic_ratio(rfile, var, select_trigger())
    r_fidcuts = get_chic_ratio(rfile, var,
                               combine_cuts([select_trigger(), get_fidcuts()]))

    leg = r.TLegend(0.5, 0.91, 0.9, 0.94)
    leg.SetNColumns(2)

    can = mkplot([r_nocut, r_fidcuts], leg=leg, yRange=[0, 1],
                 yLabel='#chi_{c2} / #chi_{c1}', xLabel=x_label,
                 legEntries=['no cuts', 'std fid cuts'])
    can.SaveAs(savename)


def get_savename(outdir, plot_name, extension):
    """
    Get the full filename for saving
    """
    save_name = '.'.join([plot_name, extension])
    save_name = '/'.join([outdir, save_name])
    return save_name


def main(args):
    """Main"""
    mcfile = r.TFile.Open(args.mcfile)

    frames = ['CS', 'HX', 'PX']
    variables = ['TMath::Abs(costh_{})', 'phi_{}_fold']
    x_labels = ['|cos#theta^{{{0}}}|', '#phi^{{{0}}}_{{folded}}']


    cond_mkdir(args.outdir)

    for frame in frames:
        for i, var in enumerate(variables):
            plot_var = var.format(frame)
            x_label = x_labels[i].format(frame)
            plot_name = '_'.join([get_plot_varname(plot_var),
                                  'fidcut', 'ratio'])
            save_name = get_savename(args.outdir, plot_name, args.extension)

            make_plot(mcfile, plot_var, x_label, save_name)

    make_plot(mcfile, 'TMath::Abs(cosalpha_HX)', '|cos#alpha^{HX}|',
              get_savename(args.outdir, 'abs_cosalpha_fidcuts_ratio',
                           args.extension))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for producing plots '
                                     'for studying effects of single muon '
                                     'fiducial cuts')
    parser.add_argument('mcfile', help='File containing reconstructed MC '
                        'events')
    parser.add_argument('-o', '--outdir', type=str, default='.',
                        help='output directory of plots')
    parser.add_argument('-e', '--extension', type=str, default='pdf',
                        help='extension of the created plots')

    clargs = parser.parse_args()
    main(clargs)
