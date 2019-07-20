#!/usr/bin/env python
"""
Script to get a corrected ratio from the fitted ratio
"""
import json

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np

from utils.two_dim_binned_fitting import BinnedFitModel
from utils.roofit_utils import eval_pdf, get_var, set_var
from utils.data_handling import get_dataframe, apply_selections
from utils.misc_helpers import select_bin
from utils.EfficiencyProvider import AcceptanceCorrectionProvider, eval_corrmap
from utils.hist_utils import project, divide, get_array, rebin, hist1d, hist2d
from utils.selection_functions import collect_requirements
from utils.graph_utils import get_errors, divide_graphs


MAP_NAME = 'fold_costh_phi_JpsiPt_JpsiRap_{}_HX'
CHI1_MODEL = 'M_chic1'
CHI2_MODEL = 'M_chic2'


def get_acc_mask(cmfile, use_pt, min_acc):
    """
    Define an acceptance (only) map mask to exclude events with very low
    acceptance values
    """
    logging.info('Masking all bins of correction map with acceptance < {}'
                 .format(min_acc))
    if use_pt:
        gen_dist = project(cmfile.Get(MAP_NAME.format('gen')), [0, 1, 2])
        acc_dist = project(cmfile.Get(MAP_NAME.format('acc')), [0, 1, 2])
    else:
        gen_dist = project(cmfile.Get(MAP_NAME.format('gen')), [1, 0])
        acc_dist = project(cmfile.Get(MAP_NAME.format('acc')), [1, 0])

    accmap = divide(acc_dist, gen_dist)

    mask = get_array(accmap) < min_acc
    return mask


def get_correction_map(cmfile, use_pt=False, acc_only=False, min_acc=None):
    """
    Get the correction map from the correction map file
    """
    # TODO: Potentially rebin this thing
    reco = 'acc' if acc_only else 'reco'
    if use_pt:
        gen_dist = project(cmfile.Get(MAP_NAME.format('gen')), [0, 1, 2])
        reco_dist = project(cmfile.Get(MAP_NAME.format(reco)), [0, 1, 2])
    else:
        gen_dist = project(cmfile.Get(MAP_NAME.format('gen')), [1, 0])
        reco_dist = project(cmfile.Get(MAP_NAME.format(reco)), [1, 0])

    accmap = divide(reco_dist, gen_dist)

    if min_acc is not None:
        acc_mask = get_acc_mask(cmfile, use_pt, min_acc)
    else:
        acc_mask = None

    return AcceptanceCorrectionProvider(accmap, mask=acc_mask)


def eval_abs_costh_phi_fold_HX(data):
    """Evaluate the correction map using abs(costh) and phi (folded)"""
    return data.costh_HX_fold.abs(), data.phi_HX_fold


def eval_abs_costh_phi_fold_pt_HX(data):
    """
    Evaluate the correction map using abs(costh) and phi (folded) as well as
    JpsiPt
    """
    return data.costh_HX_fold.abs(), data.phi_HX_fold, data.JpsiPt


def eval_costh_phi_fold_HX(data):
    """
    Evaluate the correction map using costh and phi (folded)
    """
    return data.costh_HX_fold, data.phi_HX_fold


def eval_costh_phi_fold_pt_HX(data):
    """Evaluate the correction map using costh, phi (folded) and JpsiPt"""
    return data.costh_HX_fold, data.phi_HX_fold, data.JpsiPt


def get_corr_weights(dfr, filen, use_pt=False, acc_only=False, min_acc=None):
    """
    Get the correction weights for all events using the correction map
    constructed from the histograms in the passed file
    """
    cmapfile = r.TFile.Open(filen)
    cmap = get_correction_map(cmapfile, use_pt, acc_only, min_acc)

    # Check how far down the costh goes in the correction map
    is_abs_costh = np.min(cmap.var_binnings[0]) == 0

    # Determine the appropriate eval function for the correction map
    if use_pt:
        if is_abs_costh:
            logging.info('Using corrections in abs(costh), phi and Jpsi pT')
            eval_f = eval_abs_costh_phi_fold_pt_HX
        else:
            logging.info('Using corrections in costh, phi and Jpsi pT')
            eval_f = eval_costh_phi_fold_pt_HX
    else:
        if is_abs_costh:
            logging.info('Using corrections in abs(costh) and phi')
            eval_f = eval_abs_costh_phi_fold_HX
        else:
            logging.info('Using corrections in costh and phi')
            eval_f = eval_costh_phi_fold_HX

    return eval_corrmap(cmap, eval_f)(dfr)


def load_data(filen, model):
    """
    Load the raw data making sure that only the data specified in the fitting
    range as well as in the range defined by the binning to avoid loading
    events that are not used in the fit.

    NOTE: This still does not beat a proper pre-selection!
    """
    mass_sel = select_bin(model.fit_var, *model.fit_range)
    selections = [mass_sel]
    for var, bounds in model.get_load_vars():
        selections.append(
            select_bin(var, *[float(v) for v in bounds.split(',')]))

    load_vars = ['{costh,phi}_HX_fold'] + collect_requirements(selections)

    return apply_selections(get_dataframe(filen, columns=load_vars),
                            selections)


def get_yields(wsp, binname, yield_names):
    """
    Get the sum of all the yields for a given bin
    """
    yields = []
    for name in yield_names:
        var_name = '_'.join([name, binname])
        yields.append(get_var(wsp, var_name).getVal())

    return np.sum(yields)


def get_state_fractions(bin_data, wsp, model, binname):
    """
    Get the probabilities for each event to be a chi1 or a chi2 event for a
    given bin
    """
    full_pdf = wsp.pdf('_'.join([model.full_model, binname]))
    chi1_pdf = wsp.pdf('_'.join([CHI1_MODEL, binname]))
    chi2_pdf = wsp.pdf('_'.join([CHI2_MODEL, binname]))

    full_yields = get_yields(wsp, binname, model.nevent_yields)
    chi1_yields = get_yields(wsp, binname, ['Nchic1'])
    chi2_yields = get_yields(wsp, binname, ['Nchic2'])

    mname = model.fit_var
    mvar = get_var(wsp, mname)
    full_pdf_vs = eval_pdf(full_pdf, mvar, bin_data.loc[:, mname]) * full_yields
    chi1_pdf_vs = eval_pdf(chi1_pdf, mvar, bin_data.loc[:, mname]) * chi1_yields
    chi2_pdf_vs = eval_pdf(chi2_pdf, mvar, bin_data.loc[:, mname]) * chi2_yields

    return chi1_pdf_vs / full_pdf_vs, chi2_pdf_vs / full_pdf_vs


def get_graph(wsp, model, var, sym_uncer=False):
    """
    Get the uncorrected graph to get the x-axis coordinates, since that is less
    work then redoing the whole graphing

    WARNING: This is pretty brittle and only works for costh and phi graphs
    """
    graphs = model.plot_simvar_graphs(wsp, [var], sym_uncer)
    names = [g.GetName() for g in graphs]
    idx = next((i for i, n in enumerate(names) if ('costh' in n) or ('phi' in n)), None)
    if idx is None:
        logging.error('Cannot find the prototype graph for costh or phi using '
                      'the passed model.')

    return graphs[idx]


def print_info(state, bin_dfr, state_prob, print_f=logging.info):
    """
    Print some information about rejected events, etc.
    """
    corr_w = bin_dfr.loc[:, 'corr_{}'.format(state)]
    disc_ev = corr_w == 0
    print_f('{}: {} / {} ({:.2f} %) events with 0 correction weight'.
            format(state, np.sum(disc_ev), len(disc_ev),
                   100.0 * np.sum(disc_ev) / len(disc_ev)))
    print_f('Min observed costh-phi value of all events: ({:.3f}, {:.3f})'
            .format(bin_dfr.costh_HX_fold.abs().max(),
                    bin_dfr.phi_HX_fold.max()))
    print_f('Max observed costh-phi value of used events: ({:.3f}, {:.3f})'
            .format(bin_dfr.loc[~disc_ev, 'costh_HX_fold'].abs().max(),
                    bin_dfr.loc[~disc_ev, 'phi_HX_fold'].max()))


def debug_plots(basename, rfile, state_prob, corr_w, mass):
    """
    Make some debug plots and store them into a root file so that they can be
    viewed later
    """
    rfile.cd()
    prob_h = hist1d(state_prob, min=0, max=1, name=basename + '_state_prob')
    corr_h = hist1d(corr_w, log=True, name=basename + '_corr_w')
    w_h = hist1d(state_prob * corr_w, log=True, name=basename + '_final_weight')

    # Keep the 2d histograms consistent between different bins by setting ranges
    corr_v_prob = hist2d(state_prob, corr_w, minx=0, maxx=1,
                         miny=1, maxy=1e5, logy=True, y_axis='corr_weight',
                         name=basename + '_corr_w_v_state_prob')
    w_v_prob = hist2d(state_prob, corr_w * state_prob, minx=0, maxx=1,
                      miny=0.01, maxy=1e4, logy=True, y_axis='final_weight',
                      name=basename + '_final_weight_v_state_prob')

    mass_h = hist1d(mass, min=3.2, max=3.75, weights=state_prob,
                    name=basename + '_mass_w_state_prob')

    mass_h_w = hist1d(mass, min=3.2, max=3.75, weights=state_prob * corr_w,
                      name=basename + '_mass_w_final_weight')

    rej_sp_h = hist1d(state_prob[corr_w == 0], min=0, max=1,
                      name=basename + '_state_prob_rej_ev')

    for h in [prob_h, corr_h, w_h, corr_v_prob, w_v_prob, mass_h, mass_h_w,
              rej_sp_h]:
        h.Write()


def get_corrected_ratio(data, wsp, model, sym_uncer=False, dbg_file=None):
    """
    Get the corrected ratio in all bins
    """
    logging.info('Getting corrected ratio with {} errors'.
                 format('HESSE' if sym_uncer else 'MINOS'))
    corr_ratio = []
    # NOTE: Assuming here that the bins are ordered correctly AND that the
    for label, bounds in model.bins.iteritems():
        selections = []
        for ivar, var in enumerate(model.bin_cut_vars):
            selections.append(select_bin(var, *bounds[ivar]))
        bin_data = apply_selections(data, selections)

        chi1_prob, chi2_prob = get_state_fractions(bin_data, wsp, model, label)

        print_info('chi1', bin_data, chi1_prob)
        print_info('chi2', bin_data, chi2_prob)

        if dbg_file is not None:
            debug_plots('chi1_{}'.format(label), dbg_file, chi1_prob,
                        bin_data.corr_chi1, bin_data.chicMass)
            debug_plots('chi2_{}'.format(label), dbg_file, chi2_prob,
                        bin_data.corr_chi2, bin_data.chicMass)

        chi1_w = bin_data.loc[:, 'corr_chi1'] * chi1_prob
        chi2_w = bin_data.loc[:, 'corr_chi2'] * chi2_prob
        chi1_corr = np.sum(chi1_w)
        chi2_corr = np.sum(chi2_w)
        corr_ratio.append(chi2_corr / chi1_corr)

    # Assume that the relative uncertainties are unchanged for the corrected and
    # the uncorrected graph and use them to determine the uncertainties of the
    # corrected graph
    uncorr_graph = get_graph(wsp, model, 'r_chic2_chic1', sym_uncer)
    xlo, xhi, err_lo, err_hi = get_errors(uncorr_graph)
    xvals, yvals = np.array(uncorr_graph.GetX()), np.array(uncorr_graph.GetY())
    corr_ratio = np.array(corr_ratio)

    return r.TGraphAsymmErrors(len(corr_ratio), xvals, corr_ratio, xlo, xhi,
                               err_lo / yvals * corr_ratio,
                               err_hi / yvals * corr_ratio)


def main(args):
    """Main"""
    fitfile = r.TFile.Open(args.fitresfile)
    wsp = fitfile.Get('ws_mass_fit')
    wsp.loadSnapshot('snap_two_dim')

    with open(args.fitmodelfile) as conff:
        config = json.load(conff)
    model = BinnedFitModel(config)

    data = load_data(args.datafile, model)

    n_ev_data = data.shape[0]
    n_ev_wsp = wsp.data('full_data').numEntries()
    if n_ev_data != n_ev_wsp:
        logging.warning('Number of events read in from input data file and '
                        'number of events stored in workspace are not the same')
    logging.debug('Loaded {} events, workspace contains {} events'
                  .format(n_ev_data, n_ev_wsp))

    chi1_corr_w = get_corr_weights(data, args.chi1corrfile, args.pt,
                                   args.acceptance, args.min_acc)
    chi2_corr_w = get_corr_weights(data, args.chi2corrfile, args.pt,
                                   args.acceptance, args.min_acc)
    data['corr_chi1'] = chi1_corr_w
    data['corr_chi2'] = chi2_corr_w
    # data['eff_corr'] = get_eff_corr_weights(data)

    if args.debug:
        dbgfile = r.TFile(args.outfile.replace('.root', '.dbg.root'), 'recreate')
    else:
        dbgfile = None

    # NOTE: This only works only as long as the first variable is not also used
    # for binning!
    variable = model.bin_vars[-1]
    corr_ratio = get_corrected_ratio(data, wsp, model, args.symmetric, dbgfile)
    corr_ratio.SetName('r_chic2_chic1_v_{}_bin_0'.format(variable))

    outfile = r.TFile.Open(args.outfile, 'recreate')
    outfile.cd()
    corr_ratio.Write()
    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to obtain a corrected '
                                     'ratio from a fitted ratio using correction'
                                     ' maps')
    parser.add_argument('fitresfile', help='File containing the fit results of a'
                        ' simultaneous mass fit')
    parser.add_argument('datafile', help='The file containing the (preselected) '
                        'data that have been used for the fit')
    parser.add_argument('fitmodelfile', help='Config file describing the fit '
                        'model')
    parser.add_argument('chi1corrfile', help='File containing the histograms to '
                        'construct a chic1 correction map')
    parser.add_argument('chi2corrfile', help='File containing the histograms to '
                        'construct a chic2 correction map')
    parser.add_argument('-o', '--outfile', help='File into which the corrected '
                        'ratio should be stored', default='corr_ratio.root')
    parser.add_argument('--symmetric', action='store_true', default=False,
                        help='Use HESSE uncertainties instead of MINOS '
                        'uncertainties on the ratio')
    parser.add_argument('--pt', help='Also use Jpsi pT to derive correction '
                        'weights', action='store_true', default=False)
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Create a file containing some debug plots')
    parser.add_argument('-a', '--acceptance', help='Use the acceptance only map '
                        'to derive correction weights (i.e. neglect '
                        'efficiencies)', action='store_true', default=False)
    parser.add_argument('--min-acc', type=float, help='Use a minimum ACCEPTANCE '
                        'value cut to reject events that.', default=None)


    clargs = parser.parse_args()
    main(clargs)
