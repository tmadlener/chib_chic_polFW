#!/usr/bin/env python
"""
Script for checking the chic2 / chic1 ratio at gen level after applying
different cuts to mimick cuts applied either explicitly or implicitly at
analysis level
"""

import numpy as np

# import matplotlib.pyplot as plt
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True


from root_numpy import fill_hist

from utils.data_handling import get_dataframe
from utils.misc_helpers import create_random_str, get_bin_cut_df
from utils.plot_helpers import mkplot
from utils.hist_utils import set_hist_opts


def create_histogram(plotvar, hist_sett, weights=None):
    """
    Create the n dimensional histogram using
    """
    name = create_random_str()
    # use the number of dimensions from the plotvar to determine which sort of
    # histogram to use
    nD = plotvar.shape
    if len(nD) == 1:
        nD = 1
    else:
        nD = nD[1]

    hist_type = 'TH{}D'.format(nD)
    hist = getattr(r, hist_type)(name, '', *hist_sett)
    set_hist_opts(hist)

    fill_hist(hist, plotvar, weights=weights)

    return hist


def make_plot(plotvar, histset, savename='', **kwargs):
    """
    Create (and return) the histogram and save it immediately
    """
    hist = create_histogram(plotvar, histset, kwargs.pop('weights', None))
    can = mkplot(hist, **kwargs)
    if (savename):
        can.SaveAs(savename)

    return hist


# get N x 2 array for folded costh & phi
costh_phi = lambda df, f: np.array([df['costh_{}'.format(f)].abs(), df['phi_{}_fold'.format(f)]]).T


def basic_sel(df, jpsiPt, photonPt=0):
    """
    Do a basic selection that mimics analysis and acceptance cuts
    """
    return (df.muP_pt > 3) & (df.muN_pt > 3) \
        & (df.muP_eta.abs() < 1.6) & (df.muN_eta.abs() < 1.6) \
        & (df.jpsiRap.abs() < 1.2) \
        & get_bin_cut_df(df, 'jpsiPt', *jpsiPt) \
        & (df.photonPt > photonPt)


if __name__ == '__main__':
    tuple_file = '/afs/hephy.at/work/t/tmadlener/ChiPol/chic2_chic1_ratios/tuple_batch7.root'
    gen_data = get_dataframe(tuple_file, 'chic_mc_tuple')
    chi1_data = gen_data[gen_data.wChic1 == 1]
    chi2_data = gen_data[gen_data.wChic2 == 1]

    c1_h = make_plot(costh_phi(chi1_data, 'HX'), (8, 0, 1, 10, 0, 90), 'costh_phi_chic1_nocuts_gen.svg', drawOpt='colz')
    c2_h = make_plot(costh_phi(chi2_data, 'HX'), (8, 0, 1, 10, 0, 90), 'costh_phi_chic2_nocuts_gen.svg', drawOpt='colz')
    ratio = c2_h.Clone('costh_phi_ratio')
    ratio.Divide(c1_h)
    can = mkplot(ratio, drawOpt='colz')
    can.SaveAs('costh_phi_ratio_nocuts_gen.svg')
