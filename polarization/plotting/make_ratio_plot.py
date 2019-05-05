#!/usr/bin/env python
"""
Script to make ratio plots that can at least be used in ANs
"""

import pickle
import os

from collections import OrderedDict

from scipy.stats import chi2

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.plot_helpers import mkplot, default_attributes
from utils.misc_helpers import create_random_str, parse_binning
from utils.symbolic import lth_1, lth_2
from utils.hist_utils import project, divide, rebin_1d_binning
from utils.graph_utils import divide_hist


RATIO_NAME = 'r_chic2_chic1_v_costh' # from model
XLABELS = {'phi': '#varphi^{HX}_{fold}', 'costh': '|cos#vartheta^{HX}|'}
XRANGES = {'phi': [0, 90], 'costh': [0, 1]}
YLABEL = 'N^{#chi_{c2}} / N^{#chi_{c1}}'

# For now hardcoded
YEAR = '2012'
TRIGGER = 'HLT_Dimuon8_Jpsi'

# The file containing the fitted values necessary to compute the kappas for the
# analytic ratio funcs along the phi direction
FITRESFILE = os.path.join(os.environ['CHIB_CHIC_POLFW_DIR'], 'polarization',
                          'plotting', 'fitted_lambdas_CS.pkl')
with open(FITRESFILE, 'r') as fitf:
    GEN_FIT_RES = pickle.load(fitf)


# The scenario that has lth_1 = 0.65 and lth_2 = -0.35
# Should be usable for 12 - 30 GeV pT bin (costh)
NRQCD_SCENARIO = ('chic1_R_14o73', 'chic2_R1_0_R2_25o106')

ANA_FUNC_STRINGS = {
    'phi': '[0] * (1 + [2] * cos(x * pi / 90)) / (1 + [1] * cos(x * pi / 90))',
    'costh': '[0] * (3 + [1]) / (3 + [2]) * (1 + [2] * x*x) / (1 + [1] * x*x)'
}

POL_SCENARIOS = {
    'phi': (
        # more or less grouped according to the absolute value of Delta lambda
        ('chic1_R_18o31', 'chic2_R1_0_R2_25o58'),
        ('chic1_R_22o29', 'chic2_R1_0_R2_35o62'),

        ('chic1_R_14o33', 'chic2_R1_0_R2_5o18'),
        ('chic1_R_6o7', 'chic2_R1_0_R2_134o214'),

        ('chic1_R_2o3', 'chic2_R1_2o5_R2_2o5'), # unpolarized
    ),
    'costh': (
        ('chic1_R_2o3', 'chic2_R1_2o5_R2_2o5'), # unpolarized
        ('chic1_R_14o73', 'chic2_R1_0_R2_25o106'), # nrqcd prediction

        ('chic1_R_0', 'chic2_R1_0_R2_0'), # max negative
        ('chic1_R_1', 'chic2_R1_0_R2_1') # max positive
    )
}

ANALYSIS_BINNING = {
    'phi': parse_binning("0:90,7"),
    'costh': parse_binning("0,0.075,0.15,0.225,0.3,0.375,0.45,0.625")
}

SCEN_ATTR = default_attributes(size=0, line=7)
LINE_STYLES = [7, 9, 2, 3, 5, 1]
for iatt, att in enumerate(SCEN_ATTR):
    att['color'] = 12
    att['line'] = LINE_STYLES[iatt % len(LINE_STYLES)]

RATIO_ATTR = default_attributes(open_markers=False)

LEG_POS = {
    'phi': (0.2, 0.15, 0.9, 0.33),
    'costh': (0.2, 0.15, 0.65, 0.3)
}


def get_func(variable):
    """
    Get the TF1 belonging to the passed variable, without fixing any parameters
    yet
    """
    func_str = ANA_FUNC_STRINGS[variable]
    return r.TF1(create_random_str(), func_str, *XRANGES[variable])


def get_lth_vals(chi1_scen, chi2_scen):
    """
    Get the lth values belonging to the two passed scenarios
    """
    lth1 = lth_1(R=chi1_scen.split('_')[-1].replace('o', '/'))
    lth2 = lth_2(R1=chi2_scen.split('_')[2].replace('o', '/'),
                 R2=chi2_scen.split('_')[4].replace('o', '/'))
    return lth1, lth2


def get_kappa_vals(chi1_scen, chi2_scen):
    """
    Get the kappa values belonging to the passed scenarios
    """
    def get_kappa(scenario, costh_bound=0.625):
        """
        Get the kappa value for a given scenario in symmetric costh region
        """
        if not '2o' in scenario:
            lth = GEN_FIT_RES['CS'][scenario]['lth']
            lph = GEN_FIT_RES['CS'][scenario]['lph']
        else:
            lth, lph = 0, 0

        return (3 - costh_bound**2) / (3 + lth * costh_bound**2) * lph

    return get_kappa(chi1_scen), get_kappa(chi2_scen)


def get_analytic_shape(chi1_scen, chi2_scen, variable):
    """
    Get the analytic function for the passed scenarios and variable
    """
    func = get_func(variable)
    if variable == 'phi':
        par1, par2 = get_kappa_vals(chi1_scen, chi2_scen)
    elif variable == 'costh':
        par1, par2 = get_lth_vals(chi1_scen, chi2_scen)

    func.FixParameter(1, par1)
    func.FixParameter(2, par2)

    return func


def get_analytic_scenarios(variable):
    """
    Get the analytic scenario functions
    """
    shapes = OrderedDict()
    for chi1_scen, chi2_scen in POL_SCENARIOS[variable]:
        shapes[(chi1_scen, chi2_scen)] = get_analytic_shape(chi1_scen, chi2_scen, variable)

    return shapes


def get_corr_map(corr_f, variable):
    """
    Get the correction map for one file (either chi1 or chi2) for a given
    variable
    """
    gen_h = corr_f.Get('fold_costh_phi_JpsiPt_JpsiRap_gen_HX')
    reco_h = corr_f.Get('fold_costh_phi_JpsiPt_JpsiRap_reco_HX')

    # project the histograms onto the variable and rebin them into the analysis
    # binning
    proj_dir = 0 if variable == 'costh' else 1
    gen_h = rebin_1d_binning(project(gen_h, proj_dir),
                             ANALYSIS_BINNING[variable])
    reco_h = rebin_1d_binning(project(reco_h, proj_dir),
                             ANALYSIS_BINNING[variable])
    return divide(reco_h, gen_h)


def get_correction_map(variable, chi1_f, chi2_f):
    """
    Get the correction map for the passed variable (integrated)
    """
    chi1_map = get_corr_map(chi1_f, variable)
    chi2_map = get_corr_map(chi2_f, variable)

    return divide(chi2_map, chi1_map)


def fit_to_dist(func, dist):
    """Fit the function to the distribution and return the chi2 and ndf values"""
    fitres = dist.Fit(func, 'SEX0q0')
    return fitres.Chi2(), fitres.Ndf()


def format_fr(label, fit_results):
    """
    Format the fit results and the label in a nice way
    """
    # return '{:<8} #chi^{{2}} / ndf = {:.1f} / {}'.format(label, *fit_results)
    p_val = chi2.sf(*fit_results)
    exp = np.floor(np.log10(p_val))
    if exp > -2:
        p_val_fmt = '{:.2f}'.format(p_val)
    else:
        p_val_fmt = '{:.1f} #times 10^{{{:.0f}}}'.format(p_val * 10**(-exp), exp)

    return '{:<8}, p = {}'.format(label, p_val_fmt)


def get_legend_keys(fit_res):
    """
    Get the legend keys for all the fit results
    """
    leg_entries = []

    for scen, fitr in fit_res.iteritems():
        chi1_s, chi2_s = scen # unpack scenario
        lth1 = lth_1(R=chi1_s.split('_')[-1].replace('o', '/'))
        lth2 = lth_2(R1=chi2_s.split('_')[2].replace('o', '/'),
                     R2=chi2_s.split('_')[4].replace('o', '/'))
        dlam = float(lth2 - lth1)

        leg_entries.append(
            format_fr('#Delta#lambda_{{#theta}} = {:+.1f}'.format(dlam), fitr)
        )

    return leg_entries


def make_plot_comp_ana_shapes(graph, variable):
    """
    Make the comparison plot with the analytical shapes
    """
    ana_shapes = get_analytic_scenarios(variable)
    fit_res = OrderedDict()
    for scen, shape in ana_shapes.iteritems():
        fit_res[scen] = fit_to_dist(shape, graph)

    can = mkplot(graph, drawOpt='PE', attr=RATIO_ATTR,
                 xRange=XRANGES[variable], xLabel=XLABELS[variable],
                 yRange=[0, 0.75], yLabel=YLABEL)

    mkplot(ana_shapes.values(), can=can, drawOpt='same', attr=SCEN_ATTR,
           legPos=LEG_POS[variable], legEntries=get_legend_keys(fit_res),
           legOpt='L')

    add_auxiliary_info(can, [YEAR])
    return can


def get_graph(graphfile, corrmap, corr):
    """
    Get the graph from the file and possibly apply corrections
    """
    graph = graphfile.Get(RATIO_NAME)
    if corr:
        graph = divide_hist(graph, corrmap)
    return graph


def main(args):
    """Main"""
    set_TDR_style()

    cntfile = r.TFile.Open(args.centralfile)

    # Open correction maps here, to keep histos alive
    chi1_cf = r.TFile.Open(args.chi1corrfile)
    chi2_cf = r.TFile.Open(args.chi2corrfile)

    var = args.direction
    corrmap = get_correction_map(var, chi1_cf, chi2_cf)

    graph = get_graph(cntfile, corrmap, not args.uncorr)
    can = make_plot_comp_ana_shapes(graph, var)

    if args.compgraph is not None:
        cmpfile = r.TFile.Open(args.compgraph)
        cmpgraph = get_graph(cmpfile, corrmap, not args.uncorr)

        mkplot(cmpgraph, can=can, drawOpt='samePE', attr=RATIO_ATTR[1:])

    can.SaveAs(args.outfile)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make ratio plots '
                                     'with at least AN quality')
    parser.add_argument('centralfile', help='File containing the central result'
                        ' graph')
    parser.add_argument('chi1corrfile', help='File containing the raw chi1 '
                        'inputs for the correction maps')
    parser.add_argument('chi2corrfile', help='File containing the raw chi2 '
                        'inputs for the correction maps')
    parser.add_argument('--uncorr', default=False, action='store_true',
                        help='Do not correct the ratio')
    parser.add_argument('-o', '--outfile', help='Name of the produced pdf file',
                        default='ratio.pdf')
    parser.add_argument('-c', '--compgraph', help='Add a graph from this file '
                        'as comparison (without any rescaling or fitting)',
                        default=None)

    dir_sel = parser.add_mutually_exclusive_group()
    dir_sel.add_argument('--costh', action='store_const', dest='direction',
                         const='costh',
                         help='Do ratio plots for costh direction')
    dir_sel.add_argument('--phi', action='store_const', dest='direction',
                         const='phi', help='Do ratio plots for phi direction')
    parser.set_defaults(direction='costh')

    clargs = parser.parse_args()
    main(clargs)
