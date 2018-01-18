"""
Script for producing plots for studying acceptance effects for the chic
"""

import re

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.hist_utils import (
    draw_var_to_hist, set_hist_opts, combine_cuts, set_labels, get_y_max,
    set_range_hists, get_quantiles
)
from utils.plot_helpers import default_colors, set_attributes, mkplot
from utils.misc_helpers import stringify, cond_mkdir

def get_hist_settings(var):
    """
    TODO
    """
    hist_settings = {
        'costh': {'n_bins': 16, 'min': 0, 'max': 1},
        'chicMass': {'n_bins': 200, 'min': 3.2, 'max': 3.65},
        'phi': {'n_bins': 20, 'min': 0, 'max': 90}
    }

    if 'costh' in var or 'cosalpha' in var:
        return hist_settings['costh']
    if 'phi' in var:
        return hist_settings['phi']

    return hist_settings[var]


def get_plot_varname(var):
    """
    Get the name from the variable that can be used in plot and histogram names.

    TODO
    """
    s_var = re.sub(r'TMath::Abs\((.*)\)', r'abs_\1', var)
    return s_var


def var_selection(var, v_min, v_max):
    """
    Get events in a min, max window for a given variable

    Args:
        var (str): branch name of the desired variable
        v_min, v_max (float): min and max value of the variable

    Returns:
        str: Selection string to be used in TTree::Draw()
    """
    low_cut = ' > '.join([var, str(v_min)])
    up_cut = ' < '.join([var, str(v_max)])

    return combine_cuts([low_cut, up_cut])


def select_trigger(brname='trigger', bit=2):
    """
    Select triggered events:

    Args:
        brname (str, optional): Name of the branch under which trigger info is
            stored
        bit (int, optional): bit that has to be set in order for being sleeted

    Returns
        str: Selection string to be used in TTree::Draw()
    """
    return '({0} & {1}) == {1}'.format(brname, bit)


def base_selection():
    """
    Get a basic selection that should be applied to all plots

    NOTE: Calls gets the histogram settings for 'chicMass' to determine min and
    max mass cuts

    Returns:
        str: Selection string to be used in TTree::Draw()
    """
    mass_hist_settings = get_hist_settings('chicMass')

    return var_selection('chicMass', mass_hist_settings['min'],
                         mass_hist_settings['max'])


def get_histogram(rfile, var, state, selection, treename='chic_mc_tuple'):
    """
    Get histogram of given state for variable and selection.
    TODO
    """
    tree = rfile.Get(treename)
    hist_set = get_hist_settings(var)

    hist_name = '_'.join([stringify(selection), state, get_plot_varname(var)])

    hist = r.TH1D(hist_name, '', hist_set['n_bins'],
                  hist_set['min'], hist_set['max'])

    set_hist_opts(hist)
    state_sel = {'chic1': 'wChic1', 'chic2': 'wChic2'}[state]

    draw_var_to_hist(tree, hist, var, selection, state_sel)

    return hist


def get_gen_splitted_hists(rfile, var, base_sel=[], treename='chic_mc_tuple'):
    """
    Get histogram of desired variable split into chic2 and chic1 using generator
    information.

    Args:
        rfile (ROOT.TFile): root file containing the data to be plotted
        var (str): String expression of the variable to be plotted
        hist_settings (dict): Necessary keys: 'n_bins', 'min', 'max' for setting
            the binning of the histogram
        base_sel (list of str, optional): List of cut expressions to be applied
            while drawing the histograms
        treename (str, optional): Name of the TTree in the passed root file

    Returns:
        (ROOT.TH1D, ROOT.TH1D): histogram of chic1 and chic2 separated
            distribution of the desired variable
    """
    selection = combine_cuts(base_sel)
    chic1_h = get_histogram(rfile, var, 'chic1', selection, treename)
    chic2_h = get_histogram(rfile, var, 'chic2', selection, treename)

    return (chic1_h, chic2_h)


def make_var_split_plot(rfile, var, savename, sum_dist=False, leg=False,
                        **kwargs):
    """
    Generate and save a gen_splitted plot

    TODO
    """
    chic1_h, chic2_h = get_gen_splitted_hists(rfile, var,
                                              [select_trigger(),
                                               base_selection()])

    x_label, y_label = kwargs.pop('xlabel', ''), kwargs.pop('ylabel', '')
    set_labels(chic1_h, x_label, y_label)
    set_labels(chic2_h, x_label, y_label)

    plot_hists = [chic1_h, chic2_h]
    leg_entries = ['#chi_{c1}', '#chi_{c2}']
    if sum_dist:
        sum_h = chic1_h.Clone(chic1_h.GetName().replace('chic1', 'sum'))
        sum_h.Add(chic2_h)
        plot_hists = [sum_h] + plot_hists
        leg_entries = ['sum'] + leg_entries

    y_max = get_y_max(plot_hists)
    set_range_hists(plot_hists, y_range=[0, y_max * 1.1])

    if leg:
        legend = r.TLegend(0.6, 0.91, 0.9, 0.95)
        legend.SetNColumns(3)
        legend.SetBorderSize(0)
        can = mkplot(plot_hists, leg=legend, legEntries=leg_entries, **kwargs)
    else:
        can = mkplot(plot_hists, **kwargs)

    can.SaveAs(savename)


def make_var_ratio_plot(rfile, var, savename, **kwargs):
    """
    TODO
    """
    chic1_h, chic2_h = get_gen_splitted_hists(rfile, var,
                                              [select_trigger(),
                                               base_selection()])

    ratio_h = chic2_h.Clone(chic2_h.GetName().replace('chic2', 'ratio'))
    ratio_h.Divide(chic1_h)

    ratio_h.GetYaxis().SetRangeUser(0, .75)
    x_label, y_label = kwargs.pop('xlabel', ''), kwargs.pop('ylabel', '')
    set_labels(ratio_h, x_label, y_label)

    can = mkplot([ratio_h])

    can.SaveAs(savename)


def get_quantile_lines(x_vals, y_vals, color, style):
    """
    Return TLine objects to be plotted
    """
    line1 = r.TLine(x_vals[0], y_vals[0], x_vals[0], y_vals[1])
    line2 = r.TLine(x_vals[1], y_vals[0], x_vals[1], y_vals[1])

    set_attributes(line1, color=color, line=style, width=2)
    set_attributes(line2, color=color, line=style, width=2)

    return line1, line2


def make_mass_quantile_plot(hist, quantile_pairs, color, savename):
    """
    Quantile plots for one state.

    TODO:
    """
    y_range = [0, hist.GetMaximum() * 1.05]
    set_range_hists([hist], y_range=y_range)

    can = mkplot([hist], colors=[color], drawOpt='H')

    leg = r.TLegend(0.1, 0.91, 0.9, 0.95)
    leg.SetNColumns(len(quantile_pairs))
    leg.SetBorderSize(0)

    line_styles = [2, 3, 5, 6]
    # store all lines in this list so that they live long enough to be plotted
    lines = []
    for i, pair in enumerate(quantile_pairs):
        quantiles = get_quantiles(hist, pair)
        qline_1, qline_2 = get_quantile_lines(quantiles, y_range,
                                              default_colors()[0], line_styles[i])
        qline_2.Draw()
        qline_1.Draw()

        leg.AddEntry(qline_1, '[{}, {}]'.format(*pair), 'l')

        lines.append(qline_1)
        lines.append(qline_2)

    leg.Draw()
    can.Draw()
    can.SaveAs(savename)


def make_mass_quantile_plots(mcfile, quantile_pairs, savename_base):
    """
    TODO
    """
    chic1_h, chic2_h = get_gen_splitted_hists(mcfile, 'chicMass',
                                              [select_trigger(),
                                               base_selection()])
    set_labels(chic1_h, xlabel='M^{#mu#mu#gamma}')
    set_labels(chic2_h, xlabel='M^{#mu#mu#gamma}')

    chic1_col = default_colors()[1]
    chic2_col = default_colors()[2]

    make_mass_quantile_plot(chic1_h, quantile_pairs, chic1_col,
                            savename_base + '_chic1.pdf')
    make_mass_quantile_plot(chic2_h, quantile_pairs, chic2_col,
                            savename_base + '_chic2.pdf')


def make_dist_ratio_plot_set(rfile, var, savename_base, sum_dist=False,
                             **kwargs):
    """
    TODO
    """
    save_split = savename_base + '_split.pdf'
    save_ratio = savename_base + '_ratio.pdf'

    make_var_split_plot(rfile, var, save_split, leg=True, **kwargs)
    make_var_ratio_plot(rfile, var, save_ratio, **kwargs)


def get_mass_cut_hists(rfile, var, state, mass_cuts, selection,
                            treename='chic_mc_tuple'):
    """
    Get list of histograms with mass cut applied according to the mass_cuts list
    (one histogram per mass_cut)
    """
    hists = []
    for cut in mass_cuts:
        full_sel = combine_cuts([cut, selection])
        hists.append(get_histogram(rfile, var, state, full_sel, treename))

    return hists


def get_ratio_hists(hists, baseline):
    """
    Get the ratio histograms w.r.t. to baseline
    """
    ratio_h = [h.Clone(h.GetName() + '_ratio') for h in hists]
    for hist in ratio_h:
        hist.Divide(baseline)

    return ratio_h


def make_quantile_ratio_plots(rfile, var, state, quantile_pairs, savename_base):
    """
    Make ratio plots w.r.t. no selection and quantile selection

    TODO
    """
    mass_hist = get_histogram(rfile, 'chicMass', state,
                              combine_cuts([select_trigger(),
                                            base_selection()]))

    # get the mass values for the quantiles
    q_masses = [get_quantiles(mass_hist, q) for q in quantile_pairs]
    mass_cuts = [var_selection('chicMass', q[0], q[1]) for q in q_masses]
    leg_entries = ['[{}, {}]'.format(*q) for q in quantile_pairs]
    mass_cut_hists = get_mass_cut_hists(rfile, var, state, mass_cuts,
                                        combine_cuts([base_selection(), select_trigger()]))

    # get the baseline histogram as well as the full histogram (w/o cuts
    # apart from trigger)
    base_hist = get_histogram(rfile, var, state,
                              combine_cuts([base_selection(),
                                            select_trigger()]))
    full_hist = get_histogram(rfile, var, state, select_trigger())

    leg = r.TLegend(0.1, 0.91, 0.9, 0.95)
    leg.SetNColumns(len(quantile_pairs) + 2)
    leg.SetBorderSize(0)

    can = mkplot(mass_cut_hists + [full_hist, base_hist],
                 leg=leg, legEntries=leg_entries + ['no cut', 'baseline'],
                 autoRange=True)
    can.SaveAs(savename_base + '_dists.pdf')

    full_ratios = get_ratio_hists(mass_cut_hists, full_hist)
    leg.Clear()
    leg.SetNColumns(len(quantile_pairs))
    can = mkplot(full_ratios, leg=leg, legEntries=leg_entries, autoRange=True)
    can.SaveAs(savename_base + '_ratio_full.pdf')


    base_ratios = get_ratio_hists(mass_cut_hists, base_hist)
    leg.Clear()
    can = mkplot(base_ratios, leg=leg, legEntries=leg_entries, autoRange=True)
    can.SaveAs(savename_base + '_ratio_base.pdf')




def main(args):
    """Main"""
    mcfile = r.TFile.Open(args.mcfile)
    frames = ['CS', 'HX', 'PX']
    frame_vars = ['TMath::Abs(costh_{})', 'phi_{}']
    frame_indep_vars = ['TMath::Abs(cosalpha_HX)']
    quantile_pairs = [[0.05, 0.99], [0.1, 0.95], [0.2, 0.9]]

    outdir = args.outdir
    cond_mkdir(outdir)

    # mass overview plots
    make_var_split_plot(mcfile, 'chicMass',
                        '/'.join([outdir, 'chicMass_MC_overview.pdf']),
                        leg=True, sum_dist=True, drawOpt='H')

    make_mass_quantile_plots(mcfile, quantile_pairs,
                             '/'.join([outdir, 'mass_quantile_overview']))

    # frame independent and frame dependent plots split by state
    for state in ['chic1', 'chic2']:
        for var in frame_indep_vars:
            save_base = '_'.join([get_plot_varname(var), state, 'quantile'])
            make_quantile_ratio_plots(mcfile, var, state, quantile_pairs,
                                      '/'.join([outdir, save_base]))

        for frame in frames:
            for var in frame_vars:
                fvar = var.format(frame)
                save_base = '_'.join([get_plot_varname(fvar), state, 'quantile'])
                make_quantile_ratio_plots(mcfile, fvar, state, quantile_pairs,
                                          '/'.join([outdir, save_base]))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for producing plots '
                                     'for studying acceptance effects for the '
                                     'chic')
    parser.add_argument('mcfile', help='File containing the MC generated events'
                        ' that should be used for plotting')
    parser.add_argument('-o', '--outdir', type=str, default='.',
                        help='output directory for plots')

    args = parser.parse_args()
    main(args)
