#!/usr/bin/env python
"""
Quick script to test whether the ratios (after corrections) look like they should
"""
import sys

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')


from collections import OrderedDict

from utils.data_handling import get_dataframe, apply_selections
from utils.hist_utils import hist1d, divide, get_array
from utils.misc_helpers import (
    _get_var, create_random_str, make_iterable, parse_func_var, cond_mkdir,
    parse_sel_expr
)
from utils.plot_helpers import mkplot, default_attributes, setup_legend
from utils.EfficiencyProvider import AcceptanceCorrectionProvider
from utils.symbolic import lth_1, lth_2

RATIO_SHAPES = './ratio_shapes_from_gen_samples.root'
RATIO_SHAPES_F = r.TFile.Open(RATIO_SHAPES)


pfloat = lambda x: '{:.2f}'.format(float(x)).replace('.', 'p') if x is not None else 'None'

def reco_w(data):
    reco = (data.gamma_eff_sm > 0) & (data.lepN_eff_sm > 0) & (data.lepP_eff_sm)
    return data.gamma_eff_sm * 0.01 * data.lepN_eff_sm * data.lepP_eff_sm * reco


def corr_2d(corrmap, frame):
    def corr_w(data):
        corr = corrmap.eval(_get_var(data, 'costh_' + frame), _get_var(data, 'phi_' +frame))
        print('Percentage in acceptance: ', np.sum(corr > 0) / float(len(corr)))
        return (corr > 0) * corr
    return corr_w


def corr_3d(corrmap, frame, var):
    def corr_w(data):
        corr = corrmap.eval(_get_var(data, 'costh_' + frame),
                            _get_var(data, 'phi_' + frame),
                            _get_var(data, *parse_func_var(var)))
        print('Percentage in acceptance: ', np.sum(corr > 0) / float(len(corr)))
        return (corr > 0) * corr
    return corr_w


def rcorr(corrmap, frame):
    return lambda d: reco_w(d) * corr_2d(corrmap, frame)(d)


def rcorr_3d(corrmap, frame, var):
    return lambda d: reco_w(d) * corr_3d(corrmap, frame, var)(d)


def costh_h(data, corrmap=None, frame='HX', corr_var=None):
    if corrmap is None:
        return hist1d(data.costh_HX.abs(), min=0, max=1, nbins=8,
                      weights=reco_w(data))

    if corr_var is None:
        return hist1d(data.costh_HX.abs(), min=0, max=1, nbins=8,
                      weights=rcorr(corrmap, frame)(data))

    return hist1d(data.costh_HX.abs(), min=0, max=1, nbins=8,
                  weights=rcorr_3d(corrmap, frame, corr_var)(data))


def analytical_ratio(chi1_scen, chi2_scen):
    lth1 = lth_1(R=chi1_scen.split('_')[-1].replace('o', '/'))

    lth2 = lth_2(R1=chi2_scen.split('_')[2].replace('o', '/'),
                 R2=chi2_scen.split('_')[4].replace('o', '/'))

    func_ratio = r.TF1(create_random_str(4),
                       '[0] * (3 + [1]) / (3 + [2]) * (1 + [2] * x[0]*x[0]) / (1 + [1] * x[0]*x[0])',
                       0, 1)
    func_ratio.FixParameter(1, float(lth1))
    func_ratio.FixParameter(2, float(lth2))
    func_ratio.SetParameter(0, 1)

    return func_ratio, lth2 - lth1


def ratio_shape(chi1_scen, chi2_scen):
    lth1 = lth_1(R=chi1_scen.split('_')[-1].replace('o', '/'))

    lth2 = lth_2(R1=chi2_scen.split('_')[2].replace('o', '/'),
                 R2=chi2_scen.split('_')[4].replace('o', '/'))

    shape = RATIO_SHAPES_F.Get('__'.join(['CS', chi1_scen, chi2_scen]))
    return shape, lth2 - lth1


def create_selection(sel_str):
    """
    Create a list (tuple) of selection functions
    """
    if sel_str is None:
        return None

    selections = [parse_sel_expr(s) for s in sel_str.split(',')]
    if not all(selections):
        logging.fatal('Cannot parse the passed selection: \'{}\''
                      .format(sel_str))
        sys.exit(1)

    return selections


def samples(toydatadir, scenario, sel=None):
    chi1, chi2 = scenario

    chi1_file = '/'.join([toydatadir, chi1, 'toy_data.root'])
    chi2_file = '/'.join([toydatadir, chi2, 'toy_data.root'])

    return (apply_selections(get_dataframe(chi1_file), sel),
            apply_selections(get_dataframe(chi2_file), sel))


def scale_0(hists):
    for h in make_iterable(hists):
        h.Scale(1.0 / h.GetBinContent(1))
        h.SetStats(0)
    return hists


def scale_func(func):
    norm = func.Eval(0.0625) # middle of first bin in ratios
    func.SetParameter(0, 1.0 / norm)
    return func


def get_corrmap(cfile, corr_frame, corr_var=None, mask_acc=False, acc_only=False, **kwargs):
    # acc map only
    if not mask_acc or acc_only:
        if corr_var is None:
            # ACCMAP only
            if acc_only:
                acc_map = cfile.Get('acc_map_costh_phi_{}'.format(corr_frame))
            else:
                acc_map = cfile.Get('acc_map_costh_phi_{}_eff'.format(corr_frame))
        else:
            # ACCMAP only
            if acc_only:
                acc_map = cfile.Get('acc_map_costh_phi_{}_{}'.format(corr_var, corr_frame))
            else:
                acc_map = cfile.Get('acc_map_costh_phi_{}_{}_eff'.format(corr_var, corr_frame))

        return AcceptanceCorrectionProvider(acc_map, **kwargs)

    if corr_var is None:
        acc_map = get_array(cfile.Get('acc_map_costh_phi_{}'.format(corr_frame)))
        acc_eff_name = 'acc_map_costh_phi_{}_eff'.format(corr_frame)
    else:
        acc_map = get_array(cfile.Get('acc_map_costh_phi_{}_{}'.format(corr_var, corr_frame)))
        acc_eff_name = 'acc_map_costh_phi_{}_{}_eff'.format(corr_var, corr_frame)

    acc_mask = acc_map < kwargs.pop('min_acc')
    return AcceptanceCorrectionProvider(cfile.Get(acc_eff_name), mask=acc_mask, **kwargs)


def pdfname(outbase, gen_frame, delta_lam, min_prec, min_acc):
    return '/'.join([outbase,
                     '_'.join(['comp_costh_ratio_dlam', pfloat(delta_lam),
                               'min_prec', pfloat(min_prec),
                               'min_acc', pfloat(min_acc),
                               'gen', gen_frame])]) + '.pdf'


def main(args):
    c1f = r.TFile.Open(args.corrmaps1)
    c2f = r.TFile.Open(args.corrmaps2)

    corrmap1 = get_corrmap(c1f, 'HX', args.corr_var, args.mask_acc, args.use_acc_only,
                           min_acc=args.min_acc, mask_prec=args.min_prec)
    corrmap2 = get_corrmap(c2f, 'HX', args.corr_var, args.mask_acc, args.use_acc_only,
                           min_acc=args.min_acc, mask_prec=args.min_prec)

    corrmap1_PX = get_corrmap(c1f, 'PX', args.corr_var, args.mask_acc, args.use_acc_only,
                           min_acc=args.min_acc, mask_prec=args.min_prec)
    corrmap2_PX = get_corrmap(c2f, 'PX', args.corr_var, args.mask_acc, args.use_acc_only,
                           min_acc=args.min_acc, mask_prec=args.min_prec)

    selections = create_selection(args.selection)


    scenarios = [
        ('chic1_R_2o3', 'chic2_R1_2o5_R2_2o5'),
        ('chic1_R_1', 'chic2_R1_0_R2_1'),
        ('chic1_R_0', 'chic2_R1_0_R2_0'),
    ]


    attr = default_attributes(size=1.0, open_markers=True)[:3]
    attr += [{'color': 13, 'line': 7, 'width': 2, 'marker': 1}]

    cond_mkdir(args.outbase)
    outfile = r.TFile('/'.join([args.outbase,
                                'costh_ratios_min_acc_{}_min_prec_{}_gen_{}.root'.format(pfloat(args.min_acc), pfloat(args.min_prec), args.gen_frame)]),
                      'recreate')

    for scen in scenarios:
        chi1_d, chi2_d = samples(args.toydatadir, scen, selections)

        print(scen)
        hists = OrderedDict()
        print('uncorr')
        hists['uncorr'] = (costh_h(chi2_d), costh_h(chi1_d))
        print('corr HX')
        hists['corr HX'] = (costh_h(chi2_d, corrmap2, 'HX', args.corr_var),
                            costh_h(chi1_d, corrmap1, 'HX', args.corr_var))
        print('corr PX')
        hists['corr PX'] = (costh_h(chi2_d, corrmap2_PX, 'PX', args.corr_var),
                            costh_h(chi1_d, corrmap1_PX, 'PX', args.corr_var))


        if args.gen_frame == 'HX':
            ratio_a, delta_lth = analytical_ratio(*scen)
        if args.gen_frame == 'CS':
            ratio_a, delta_lth = ratio_shape(*scen)

        ratios = OrderedDict()
        for case in hists:
            ratios[case] = divide(*hists[case])

            leg = setup_legend(0.15, 0.15, 0.55, 0.25)
            leg.SetNColumns(2)


        plots = scale_0(ratios.values())
        leg_entries = ratios.keys()
        if args.gen_frame == 'HX':
            plots.append(scale_func(ratio_a))
            leg_entries.append('analytical')

        can = mkplot(plots,
                     attr=attr,
                     xLabel='|cos#vartheta^{HX}|', yLabel='#chi_{c2} / #chi_{c1}', xRange=[0, 1],
                     leg=leg, legEntries=leg_entries,
                     drawOpt='E1',
                     legOpt='PLE')
        if args.gen_frame == 'CS':
            can = mkplot(scale_0(ratio_a), attr=[attr[-1]], leg=leg, legEntries=['analytical'],
                         drawOpt='HISTsame', legOpt='PLE',
                         can=can)

        can.SaveAs(pdfname(args.outbase, args.gen_frame, delta_lth, args.min_prec, args.min_acc))

        outfile.cd()
        for corr, ratio in ratios.iteritems():
            ratio.SetName('_'.join([
                scen[0], scen[1], 'ratio_costh_HX', corr.replace(' ', '_'), 'gen', args.gen_frame
            ]))
            ratio.Write()

        for corr in hists:
            hist1, hist2 = hists[corr]
            hist1.SetName('_'.join([
                scen[0], scen[1], 'chi1_costh_HX', corr.replace(' ', '_'), 'gen', args.gen_frame
            ]))
            hist1.Write()

            hist2.SetName('_'.join([
                scen[0], scen[1], 'chi2_costh_HX', corr.replace(' ', '_'), 'gen', args.gen_frame
            ]))
            hist2.Write()


    outfile.Close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Script to produce corrected '
                                     'costh ratios')
    parser.add_argument('corrmaps1', help='correction maps for chic1')
    parser.add_argument('corrmaps2', help='correction maps for chic2')
    parser.add_argument('toydatadir', help='toy chic1/2 data directory that is '
                        'to be used for testing (must contain the three '
                        'scenarios where each folder must contain a '
                        'toy_data.root file)')
    # parser.add_argument('toydata2dir', help='toy chic2 data directory that is '
    #                     'to be used for testing (must contain the three '
    #                     'scenarios where each folder must contain a '
    #                     'toy_data.root file)')
    parser.add_argument('-f', '--gen-frame', help='gen frame', default='HX')
    parser.add_argument('-p', '--min-prec', help=' min precision', default=None,
                        type=float)
    parser.add_argument('-a', '--min-acc', help='minimal acceptance', default=0.0,
                        type=float)
    parser.add_argument('-v', '--corr-var', help='correction variable (3d only)',
                        default=None, type=str)
    parser.add_argument('-o', '--outbase', help='output base directory',
                        default='.')
    parser.add_argument('-ma', '--mask-acc', help='Use the acceptance only map'
                        ' for masking but the reco map for corrections',
                        default=False, action='store_true')
    parser.add_argument('--use-acc-only', action='store_true', default=False,
                        help='Use the acceptance only maps (no eff weighting)')
    parser.add_argument('-s', '--selection', help='Comma separated list of '
                        'selections that should be applied to the data before '
                        'filling the histograms. NOTE: you can only use the '
                        '\'<\' sign to formulate your selection (but it can be '
                        'two sided)', default=None)

    clargs = parser.parse_args()
    main(clargs)
