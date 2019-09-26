#!/usr/bin/env python
"""
Script to make a comparison of the ratio with some analytical shapes
"""
from collections import OrderedDict

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, default_attributes, setup_legend
from utils.setup_plot_style import add_auxiliary_info, set_TDR_style
from utils.plot_decoration import YLABELS, VAR_PLOT
from utils.misc_helpers import cond_mkdir_file, create_random_str, fmt_float
from utils.graph_utils import pull_graph

import utils.misc_helpers
utils.misc_helpers.MAX_N_DIGITS = 1

from common_utils import get_graph

ANALYTIC_ATTRS = default_attributes(size=0, line=7)
LINE_STYLES = [7, 9, 2, 3, 5, 1]
for iatt, att in enumerate(ANALYTIC_ATTRS):
    att['line'] = LINE_STYLES[iatt % len(LINE_STYLES)]



def r_costh(set_vals, fix_lambdas=True):
    """
    R(costh) with possible fixed values
    """
    val_to_idx = {'norm': 0, 'lth1': 1, 'lth2': 2}
    func = r.TF1(create_random_str(),
                 '[0] * (1 + [2] * x[0]*x[0]) / (1 + [1] * x[0]*x[0])',
                 0, 1)
    def_vals = {'norm': 1.0, 'lth1': 0, 'lth2': 0}
    def_vals.update(set_vals)

    for name, val in def_vals.iteritems():
        func.SetParameter(val_to_idx[name], val)

    if fix_lambdas:
        for ipar in [1, 2]:
            func.FixParameter(ipar, func.GetParameter(ipar))

    return func


def get_costh_label(lam1, lam2):
    """
    Get a label for costh plots
    """
    fstr = '#Delta#lambda_{{#vartheta}} = {}, #lambda_{{#vartheta}}(#chi_{{c1}}) = {}'
    return fstr.format(fmt_float(lam2 - lam1), fmt_float(lam1))


def fit_to_dist(func, dist):
    """Fit the function to the distribution and return the chi2 and ndf values"""
    fitres = dist.Fit(func, 'SEX0q0')
    return fitres.Chi2(), fitres.Ndf()


def get_analytic_scenarios_costh(dlth=None, lam1=None, lam2=None,
                                 no_extremes=False):
    """
    Get analytic shapes for different polarization scenarios

    If dlth is not None, the value of dlth will be used with different values
    for lambda 1. If it is None, then the extreme polarization scenarios will be
    used instead

    If lam1 AND lam2 are passed a scenario specifically for these values will
    also be added
    """
    if dlth is not None:
        lth1 = [0, 0.5, 1.0]
        lth2 = [dlth + v for v in lth1]
    else:
        if no_extremes:
            lth1, lth2 = [], []
        else:
            lth1 = [0, 1, -1./3.]
            lth2 = [0, -3./5., 1]

    if lam1 is not None and lam2 is not None:
        lth1 = [0, lam1]
        lth2 = [0, lam2]

    funcs = OrderedDict()
    for l1, l2 in zip(lth1, lth2):
        label = get_costh_label(l1, l2)
        funcs[label] = r_costh({'lth1': l1, 'lth2': l2})

    return funcs


def get_chi2_ndf(ratio, analytic_funcs):
    """
    Fit all analytic functions to the ratio and get the chi2 / ndf for all of
    them
    """
    results = OrderedDict()
    for label, func in analytic_funcs.iteritems():
        results[label] = fit_to_dist(func, ratio)
        print label, results[label]
    return results


def make_ratio_plot(ratio, analytic_scens, variable):
    """
    Make the comparison plot of the ratio and the analytical scenarios
    """
    can = mkplot(ratio, drawOpt='PE',
                 attr=[{'color': 1, 'size': 1, 'marker': 20}],
                 yLabel=YLABELS.get('r_chic2_chic1'), yRange=[0, 0.75],
                 **VAR_PLOT[variable])

    chi2_ndfs = get_chi2_ndf(ratio, analytic_scens)

    leg = setup_legend(0.2, 0.15, 0.65, 0.3)
    mkplot(analytic_scens.values(), drawOpt='sameL', can=can, attr=ANALYTIC_ATTRS,
           leg=leg, legEntries=analytic_scens.keys(), legOpt='L')

    add_auxiliary_info(can, 2012, prelim=True)
    return can


def make_ratio_pull_plot(ratio, analytic_scens, variable):
    """
    Make a comparison plot with an additional pull panel at the bottom
    """
    # Scale all analytic scenarios first
    chi2_ndfs = get_chi2_ndf(ratio, analytic_scens)

    can = mkplot([]) # get a canvas that is already a TCanvasWrapper
    can.cd()
    rpad = r.TPad('ratio_pad', 'ratio pad', 0, 0.3, 1, 0.98)
    r.SetOwnership(rpad, False)
    rpad.SetBottomMargin(0)
    rpad.Draw()

    rpad = mkplot(ratio, drawOpt='PE', can=rpad,
                  attr=[{'color': 1, 'size': 1, 'marker': 20}],
                  yLabel=YLABELS.get('r_chic2_chic1'), yRange=[0.001, 0.75],
                  **VAR_PLOT[variable])

    leg = setup_legend(0.2, 0.15, 0.65, 0.32)
    mkplot(analytic_scens.values(), drawOpt='sameL', can=rpad,
           attr=ANALYTIC_ATTRS,
           leg=leg, legEntries=analytic_scens.keys(), legOpt='L')


    can.cd()
    ppad = r.TPad('pull_pad', 'pull pad', 0, 0, 1, 0.3)
    r.SetOwnership(ppad, False)
    ppad.SetTopMargin(0)
    ppad.SetBottomMargin(0.3)
    ppad.Draw()

    # Draw grid first so that it is not overlaying the plotted points
    xran = VAR_PLOT[variable]['xRange']
    ppad = mkplot(r.TLine(xran[0], 0, xran[1], 0), can=ppad,
                  attr=[{'color': 1, 'width': 2, 'line': 1}],
                  yLabel='pulls', yRange=[-4.99, 4.99],# grid=True,
                  **VAR_PLOT[variable])
    mkplot([r.TLine(xran[0], v, xran[1], v) for v in [-3, 3]], drawOpt='same',
           can=ppad, attr=[{'color': 12, 'width': 2, 'line': 7}])

    pulls = [pull_graph(ratio, s) for s in analytic_scens.values()]
    ppad = mkplot(pulls, drawOpt='samePE', can=ppad,
                  attr=default_attributes(open_markers=False),
                  **VAR_PLOT[variable])


    phist = ppad.pltables[0]
    phist.SetNdivisions(505, "XYZ")
    phist.SetLabelFont(42, "XYZ")
    phist.SetLabelOffset(0.007, "XYZ")
    phist.SetLabelSize(0.12, "XYZ")
    phist.GetXaxis().SetTitleOffset(0.95)
    phist.GetYaxis().SetTitleOffset(0.52)
    phist.SetTitleSize(0.14, "XYZ")
    # ppad.SetGridy()

    can.add_tobject(rpad)
    can.add_tobject(ppad)


    add_auxiliary_info(can, 2012, prelim=True)
    return can


def main(args):
    """Main"""
    set_TDR_style()

    rfile = r.TFile.Open(args.ratiofile)
    ratio = get_graph(rfile, 'r_chic2_chic1', args.variable)

    if args.variable == 'costh':
        analytic_scens = get_analytic_scenarios_costh(args.delta_lambda,
                                                      args.lambda1, args.lambda2,
                                                      args.no_extremes)
    else:
        # not implemented currently for phi direction
        analytic_scens = OrderedDict()

    if args.pull_plot:
        can = make_ratio_pull_plot(ratio, analytic_scens, args.variable)
    else:
        can = make_ratio_plot(ratio, analytic_scens, args.variable)

    cond_mkdir_file(args.outfile)
    can.SaveAs(args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make plots comparing'
                                     ' the ratio with analytical scenarios')
    parser.add_argument('ratiofile', help='File containing the ratio')
    parser.add_argument('-o', '--outfile', help='Directory into which the '
                        'produced plot is placed',
                        default='r_chic2_chic1_comp_analytic.pdf')
    parser.add_argument('-d', '--delta-lambda', help='Delta lambda to be used '
                        'for the plots', default=None, type=float)
    parser.add_argument('--lambda1', help='Lambda1 that should be used in the '
                        'plot (can be used to enter an additional NRQCD'
                        'scenario)', default=None, type=float)
    parser.add_argument('--lambda2', help='Lambda2 that should be used in the '
                        'plot (can be used to enter an additional NRQCD'
                        'scenario)', default=None, type=float)
    parser.add_argument('-n', '--no-extremes', help='Do not plot the extreme '
                        'scenarios', default=False, action='store_true')
    parser.add_argument('-p', '--pull-plot', help='Plot a pull panel',
                        action='store_true', default=False)

    var_sel = parser.add_mutually_exclusive_group()
    var_sel.add_argument('--costh', action='store_const', dest='variable',
                         const='costh', help='Input ratios are vs costh')
    var_sel.add_argument('--phi', action='store_const', dest='variable',
                         const='phi', help='Input ratios are vs phi')
    parser.set_defaults(variable='costh')

    clargs = parser.parse_args()
    main(clargs)
