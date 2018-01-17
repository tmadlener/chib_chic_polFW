"""
Script for producing plots for studying acceptance effects for the chic
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.hist_utils import (
    draw_var_to_hist, set_hist_opts, combine_cuts, set_labels
)
from utils.plot_helpers import plot_on_canvas


def get_hist_settings(var):
    """
    TODO
    """
    hist_settings = {
        'costh': {'n_bins': 16, 'min': 0, 'max': 1},
        'chicMass': {'n_bins': 100, 'min': 3.275, 'max': 3.725},
        'phi': {'n_bins': 20, 'min': 0, 'max': 90}
    }

    if 'costh' in var:
        return hist_settings['costh']
    if 'phi' in var:
        return hist_settings['phi']

    return hist_settings[var]


def select_state(state, varname='gen_chic_p4.M()'):
    """
    Select the state via its generated mass.

    Args:
        state (str): Can be either 'chic1', or 'chic2'
        varname (str, optional): The name of the branch in which the generator
            4-vector is stored

    Returns:
        str: Selection string to be used in TTree::Draw()
    """
    # All chic are generated with the same mass (distinct values for chi1 and
    # chic2). For deciding if something is a chic1 or chic2 at generator level
    # define a mass inbetween the two and check if it's above or below
    mass_split = 3.53 # GE
    operator = {'chic1': ' < ', 'chic2': ' > '}

    return operator[state].join([varname, str(mass_split)])


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


def get_gen_splitted_hists(rfile, var, hist_settings, base_sel=[],
                           treename='chic_mc_tuple'):
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
    def get_name(state, var):
        """Get unique name for histogram from state and var.

        If the var contains function calls, etc. they are replaced with a
        syntax that doesn't confuse the TTree::Draw() command

        TODO
        """
        nvar = var.replace('()', '')
        nvar = nvar.replace('.', '_')
        nvar = nvar.replace('TMath::', '')
        nvar = nvar.replace('(', '_').replace('', '_')

        return '__'.join([nvar, state])

    n_bins, h_min, h_max = [hist_settings[v] for v in ['n_bins', 'min', 'max']]

    tree = rfile.Get(treename)
    chic1_h = r.TH1D(get_name('chic1', var), '', n_bins, h_min, h_max)
    chic2_h = r.TH1D(get_name('chic2', var), '', n_bins, h_min, h_max)

    set_hist_opts(chic1_h)
    set_hist_opts(chic2_h)

    chic1_sel = combine_cuts(base_sel)
    chic2_sel = combine_cuts(base_sel)

    draw_var_to_hist(tree, chic1_h, var, chic1_sel, 'wChic1')
    draw_var_to_hist(tree, chic2_h, var, chic2_sel, 'wChic2')

    return (chic1_h, chic2_h)


def make_var_split_plot(rfile, var, savename, sum_dist=False, leg=False,
                        **kwargs):
    """
    Generate and save a gen_splitted plot

    TODO
    """
    chic1_h, chic2_h = get_gen_splitted_hists(rfile, var,
                                              get_hist_settings(var),
                                              [select_trigger()])


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

    can = r.TCanvas('_'.join([var, 'can']), '', 50, 50, 600, 600)

    if leg:
        legend = r.TLegend(0.6, 0.91, 0.9, 0.95)
        legend.SetNColumns(3)
        legend.SetBorderSize(0)
        can = plot_on_canvas(can, plot_hists, leg=legend,
                             legEntries=leg_entries, **kwargs)
    else:
        can = plot_on_canvas(can, plot_hists, **kwargs)

    can.SaveAs(savename)


def make_var_ratio_plot(rfile, var, savename, **kwargs):
    """
    TODO
    """
    chic1_h, chic2_h = get_gen_splitted_hists(rfile, var,
                                              get_hist_settings(var),
                                              [select_trigger()])


    ratio_h = chic2_h.Clone(chic2_h.GetName().replace('chic2', 'ratio'))
    ratio_h.Divide(chic1_h)

    x_label, y_label = kwargs.pop('xlabel', ''), kwargs.pop('ylabel', '')
    set_labels(ratio_h, x_label, y_label)

    can = r.TCanvas('_'.join([var, 'ratio', 'can']), '', 50, 50, 600, 600)
    can = plot_on_canvas(can, [ratio_h])

    can.SaveAs(savename)


def main(args):
    """Main"""
    mcfile = r.TFile.Open(args.mcfile)
    # make_mass_plot(mcfile, 'refit_mass_reco_gen_split.pdf')
    make_var_split_plot(mcfile, 'TMath::Abs(costh_HX)', 'costh_HX_split_gen.pdf',
                        leg=True, xlabel='|cos#theta^{HX}|')
    make_var_split_plot(mcfile, 'chicMass', 'refit_mass_reco_gen_split.pdf',
                        True, True, drawOpt='H',
                        xlabel = 'M^{#mu#mu#gamma}')
    make_var_ratio_plot(mcfile, 'TMath::Abs(costh_HX)', 'costh_HX_ratio_gen.pdf',
                        xlabel='|cos#theta^{HX}|', ylabel='#chi_{c2} / #chi_{c1}')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for producing plots '
                                     'for studying acceptance effects for the '
                                     'chic')
    parser.add_argument('mcfile', help='File containing the MC generated events'
                        ' that should be used for plotting')

    args = parser.parse_args()
    main(args)
