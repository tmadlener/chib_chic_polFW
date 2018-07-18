#!/usr/bin/env python
"""
Script for checking the effects of single muon acceptance cuts on the costh
distribution in genlevel mc
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gStyle.SetPadRightMargin(0.15)
r.gStyle.SetPadLeftMargin(0.12)
# r.gROOT.SetBatch()

from gen_level_ratios import costh_phi
from data_single_muon_acc import (
    basic_sel, jpsi_kin_sel, fiducial_cuts, loose_cuts, flat_pt
)

from utils.data_handling import get_dataframe, apply_selections
form utils.hist_utils import create_histogram
from utils.misc_helpers import get_bin_cut_df, cond_mkdir
from utils.plot_helpers import mkplot, default_colors
from utils.selection_functions import get_cut_funcs

# which selections to plot?
# NOTE: basic_sel does not select anything for gen mc
plot_selections = [
    'no_sel', 'no_sel_fiducial', 'no_sel_loose', 'no_sel_flat',
    'jpsi_kin_sel', 'jpsi_kin_sel_loose', 'jpsi_kin_sel_fiducial', 'jpsi_kin_sel_flat',
    # 'basic_sel', 'basic_sel_fiducial', 'basic_sel_loose',
    # 'basic_jpsi_kin_sel', 'basic_jpsi_kin_sel_fiducial', 'basic_jpsi_kin_sel_loose'
]

# costh binning as on 2012 data (without any selection)
costh_binning_data = np.array(
    [0,
     0.07517317766951126,
     0.15925103304632882,
     0.26912794046647831,
     0.70230451686879192,
     1.0] # overflow bin
)


def single_muon_sel(df, cuts):
    """
    Single muon cuts applied to both Mons
    """
    return apply_singlemuon_cuts(df.muP_pt, df.muP_eta.abs(), cuts) \
        & apply_singlemuon_cuts(df.muN_pt, df.muN_eta.abs(), cuts)


def get_selections(jpsiPt):
    # some shorthands for less typing below
    base_sel = lambda df: basic_sel(df, True, True)
    jpsi_sel = lambda df: jpsi_kin_sel(df, jpsiPt)
    loose = lambda df: single_muon_sel(df, loose_cuts())
    fiducial = lambda df: single_muon_sel(df, fiducial_cuts())
    flat = lambda df: single_muon_sel(df, flat_pt(3.5))

    selection_sets = {
        'no_sel': None,
        'no_sel_fiducial': (fiducial,),
        'no_sel_loose': (loose,),
        'no_sel_flat': (flat,),

        'basic_sel': (base_sel,),
        'basic_sel_fiducial': (base_sel, fiducial),
        'basic_sel_loose': (base_sel, loose),

        'basic_jpsi_kin_sel': (base_sel, jpsi_sel),
        'basic_jpsi_kin_sel_fiducial': (base_sel, jpsi_sel, fiducial),
        'basic_jpsi_kin_sel_loose': (base_sel, jpsi_sel, loose),

        'jpsi_kin_sel': (jpsi_sel,),
        'jpsi_kin_sel_loose': (jpsi_sel, loose),
        'jpsi_kin_sel_fiducial': (jpsi_sel, fiducial),
        'jpsi_kin_sel_flat': (jpsi_sel, flat)
    }

    return selection_sets


# def get_attributes():
#     cols = default_colors()
#     attribute_sets = {
#         'no_sel': {'color': cols[0], 'marker': 24, 'size': 1.5},
#         'no_sel_fiducial': {'color': cols[0], 'marker': 25, 'size': 1.5},
#         'no_sel_loose': {'color': cols[0], 'marker': 26, 'size': 1.5},
#         'no_sel_flat': {'color': cols[0], 'marker': 32, 'size': 1.5},

#         'jpsi_kin_sel': {'color': cols[1], 'marker': 24, 'size': 1.5},
#         'jpsi_kin_sel_fiducial': {'color': cols[1], 'marker': 25, 'size': 1.5},
#         'jpsi_kin_sel_loose': {'color': cols[1], 'marker': 26, 'size': 1.5},
#         'jpsi_kin_sel_flat': {'color': cols[1], 'marker': 32, 'size': 1.5},

#         'basic_jpsi_kin_sel': {'color': cols[2], 'marker': 24, 'size': 1.5},
#         'basic_jpsi_kin_sel_fiducial': {'color': cols[2], 'marker': 25, 'size': 1.5},
#         'basic_jpsi_kin_sel_loose': {'color': cols[2], 'marker': 26, 'size': 1.5},

#         'basic_sel': {'color': cols[3], 'marker': 24, 'size': 1.5},
#         'basic_sel_fiducial': {'color': cols[3], 'marker': 25, 'size': 1.5},
#         'basic_sel_loose': {'color': cols[3], 'marker': 26, 'size': 1.5},
#     }

#     return attribute_sets


def get_costh_phi_hist(df, frame, selections, gen=False):
    """
    Make a 2D costh phi hist from the passed data
    """
    return create_histogram(costh_phi(apply_selections(df, selections), frame, gen),
                            (8, 0, 1, 10, 0, 90))


def get_costh_hist(df, frame, selections, gen=False):
    """
    Make a 1D costh hist from the passed data
    """
    return create_histogram(costh_phi(apply_selections(df, selections), frame, gen)[:,0],
                            # (8, 0, 1),
                            (len(costh_binning_data) - 1, costh_binning_data),
                            x_axis='|cos#theta^{{{}}}|'.format(frame))


def get_phi_hist(df, frame, selections, gen=False):
    """
    Make a 1D phi hist from the passed data
    """
    return create_histogram(costh_phi(apply_selections(df, selections), frame, gen)[:,1],
                            (10, 0, 90),
                            x_axis='|#phi^{{{}}}|'.format(frame))


def get_all_hists(df, frame, selections):
    """
    Get all 1D and 2D histograms for a given selection
    """
    return {
        'costh_phi': get_costh_phi_hist(df, frame, selections),
        'costh': get_costh_hist(df, frame, selections),
        'phi': get_phi_hist(df, frame, selections),
        'jpsipt': get_jpsipt_hist(df, selections)
    }



def get_jpsipt_hist(df, selections):
    return create_histogram(apply_selections(df, selections).jpsiPt,
                            (30, 5, 35))


def get_ratios(chi1_hists, chi2_hists):
    """
    Get the chi2 / chi1 ratios
    """
    ratios = {}
    for hist in chi1_hists:
        ratios[hist] = chi2_hists[hist].Clone()
        ratios[hist].Divide(chi1_hists[hist])

    return ratios


def save_histogram(hists, savename, **kwargs):
    """Save histogram(s) via mkplot"""
    mkplot(hists, **kwargs).SaveAs(savename)
    # return mkplot(hists, **kwargs)

def save_2D_histograms(hists, outdir, prefix, **kwargs):
    """
    Save all histograms to the output directory
    """
    basename = '{}/{}'.format(outdir, prefix)

    for selection in hists:
        savename = '{}_{}_costh_phi.pdf'.format(basename, selection)
        save_histogram(hists[selection]['costh_phi'], savename, **kwargs)


def save_1D_histograms(hists, outdir, prefix, **kwargs):
    """
    Make a comparison plot with all 1D histograms
    """
    xlabels = {'costh': '|cos#theta|', 'phi': '#phi', 'jpsipt': 'p_{T}^{J/#psi}'}

    def select_hists(hists, h):
        # sel_hists = [
        #     hists['no_sel'][h], hists['jpsi_kin_sel'][h],
        #     hists['no_sel_loose'][h], hists['jpsi_kin_sel_loose'][h],
        #     hists['no_sel_flat'][h], hists['jpsi_kin_sel_flat'][h],
        #     hists['no_sel_fiducial'][h], hists['jpsi_kin_sel_fiducial'][h]
        # ]
        sel_hists = [
            hists['jpsi_kin_sel'][h],
            hists['jpsi_kin_sel_loose'][h],
            hists['jpsi_kin_sel_flat'][h],
            hists['jpsi_kin_sel_fiducial'][h]
        ]
        return sel_hists


    def get_attributes():
        """attributes for plotting nicely in accordance with order of histograms as
        returned by select_histsk"""
        defc = default_colors()
        # return [
        #     {'color': defc[0], 'marker': 20, 'size': 1.5}, {'color': defc[0], 'marker': 24, 'size': 1.5},
        #     {'color': defc[1], 'marker': 21, 'size': 1.5}, {'color': defc[1], 'marker': 25, 'size': 1.5},
        #     {'color': defc[2], 'marker': 22, 'size': 1.5}, {'color': defc[2], 'marker': 26, 'size': 1.5},
        #     {'color': defc[3], 'marker': 23, 'size': 1.5}, {'color': defc[3], 'marker': 32, 'size': 1.5}
        # ]
        return [
            {'color': defc[0], 'marker': 20, 'size': 1.5},
            {'color': defc[1], 'marker': 21, 'size': 1.5},
            {'color': defc[2], 'marker': 22, 'size': 1.5},
            {'color': defc[3], 'marker': 23, 'size': 1.5}
        ]


    def get_legentries():
        """Get leg entries with same order as select_hists"""
        return ['no cuts', 'loose', 'flat', 'fiducial']
        # return ['no sel', 'jpsi kin sel',
        #         '', 'loose', '', 'flat', '', 'fiducial']

    leg_entries =  get_legentries()

    leg_pos = (0.12, 0.1, 0.42, 0.3) if 'ratio' in prefix else (0.5, 0.4, 0.8, 0.6)
    leg = r.TLegend(*leg_pos)
    leg.SetNColumns(2)

    plot_attr = get_attributes()
    for hist in ['costh', 'phi']:
    # for hist in ['jpsipt']:
        leg.Clear()
        plot_hists = select_hists(hists, hist) # very good naming
        # return save_histogram(plot_hists, '{}/{}_{}_gen_mc.pdf'.format(outdir, prefix, hist),
        save_histogram(plot_hists, '{}/{}_{}_gen_mc.pdf'.format(outdir, prefix, hist),
                       leg=leg, legEntries=leg_entries, legOpt='PLE', attr=plot_attr,
                       xLabel=xlabels[hist],
                       **kwargs)


def make_comp_plots(df, frame, outdir, jpsiPt):
    """
    Make a set of histograms and put them onto some plots for comparison
    """
    cond_mkdir(outdir)
    selection_sets = get_selections(jpsiPt)

    chi1_data = df[df.wChic1 == 1]
    chi2_data = df[df.wChic2 == 1]

    chi1_hists, chi2_hists, ratio_hists = {}, {}, {}
    for selection in plot_selections:
        print selection
        chi1_hists[selection] = get_all_hists(chi1_data, frame,
                                              selection_sets[selection])
        chi2_hists[selection] = get_all_hists(chi2_data, frame,
                                              selection_sets[selection])
        ratio_hists[selection] = get_ratios(chi1_hists[selection],
                                            chi2_hists[selection])

    # save_2D_histograms(chi1_hists, outdir, 'chic1', drawOpt='colz',
    #                    xLabel='|cos#theta^{{{}}}|'.format(frame),
    #                    yLabel='#phi^{{{}}}'.format(frame))
    # save_2D_histograms(chi2_hists, outdir, 'chic2', drawOpt='colz',
    #                    xLabel='|cos#theta^{{{}}}|'.format(frame),
    #                    yLabel='#phi^{{{}}}'.format(frame))
    # save_2D_histograms(ratio_hists, outdir, 'ratio', drawOpt='colz',
    #                    xLabel='|cos#theta^{{{}}}|'.format(frame),
    #                    yLabel='#phi^{{{}}}'.format(frame))

    save_1D_histograms(ratio_hists, outdir, 'ratio', drawOpt='E1',
                       yRange=[0, None], yLabel='#chi_{c2} / #chi_{c1}')
    # save_1D_histograms(chi1_hists, outdir, 'chic1', drawOpt='E1',
    #                    yRange=[0, None], yLabel='N^{#chi_{c1}}')
    # save_1D_histograms(chi2_hists, outdir, 'chic2', drawOpt='E1',
    #                    yRange=[0, None], yLabel='N^{#chi_{c2}}')

    # can_chi1 = save_1D_histograms(chi1_hists, outdir, 'jpsipt', drawOpt='E1',
    #                               yRange=[0, None])
    # can_chi2 = save_1D_histograms(chi2_hists, outdir, 'jpsipt', drawOpt='E1',
    #                               yRange=[0, None])





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for checking effects '
                                     'of single muon acceptance cuts on costh '
                                     'distribution at gen level')
    parser.add_argument('genfile', help='File containing the flat genlevel tuple')
    parser.add_argument('--ptmin', type=float, default=8, help='minimum jpsi pt')
    parser.add_argument('--ptmax', type=float, default=20, help='maximum jpsi pt')
    parser.add_argument('-t', '--treename', default='chic_mc_tuple', type=str,
                        help='name of the tree in which the original variables are '
                        '(used for storing output file).')
    parser.add_argument('-f', '--frame', help='reference frame', default='HX')
    parser.add_argument('-o', '--outdir', help='output directory', default='.')

    args = parser.parse_args()

    data = get_dataframe(args.genfile, args.treename)
    make_comp_plots(data, args.frame, args.outdir, (args.ptmin, args.ptmax))
