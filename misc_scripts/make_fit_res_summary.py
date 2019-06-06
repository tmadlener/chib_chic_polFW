#!/usr/bin/env python
"""
Script to generate a tex stub that can be included in a .tex file and gives
a "reasonable" appearance
"""

from __future__ import print_function

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import glob
import json
import re
from os import path

from collections import OrderedDict

from utils.reporting import open_tex_file, create_figure
from utils.plot_decoration import PLOT_LABELS_LATEX
from utils.misc_helpers import parse_binning, parse_func_var
from utils.roofit_utils import param_str

import utils.RooDoubleCB
import utils.RooErfExponential
import utils.RooPowerlawExponential


NICE_VARS = {'phi_HX_fold': r'\varphi^{HX}_{\text{fold}}', 'costh_HX_fold':
             r'|\cos\vartheta^{HX}|'}
LABEL_BASE = {'phi_HX_fold': r'${1:.0f} < {0} < {2:.0f}$',
              'costh_HX_fold': r'${1:.3f} < {0} < {2:.3f}$'}

# maximum number of figures per page
MAX_FIG_P_PAGE = 6
SORT_RGX = re.compile(r'_([0-9]+)_massfit\.pdf')

# Define an order such that the similar plots always appear next to each other
DESIRED_PLOT_ORDER = [
    'CBmass1', 'CBmass2',
    'CBsigma1', 'CBsigma2',
    'mu_bkg', 'sigma_bkg',
    'lambda_bkg', 'Nbkg',
    'Nchic0', 'r_chic0_chic1'
    'Nchic1', 'Nchic2',
]

# Parameters that should not be put into the fitted / fixed parameter table
EXCLUDE_PARS = [
    'm_psiPDG', 'm_chic0PDG', 'm_chic1PDG', 'm_chic2PDG', 'CBwidth'
]

def get_bin_idx(plotname):
    """
    Get the bin index from the plot name
    """
    return int(SORT_RGX.search(plotname).group(1))


def get_binning(config):
    """
    Read the fit config json and get the bin variable and the binning
    """
    binning = parse_binning(config['binning'][1])
    bins = zip(binning[:-1], binning[1:])

    return parse_func_var(config['bin_vars'][1])[0], bins


def get_fit_res_names(indir, bin_var, binning):
    """
    Get the list of fit_res files
    """
    plot_files = glob.glob(path.join(indir, 'JpsiPt_0_*_HX_fold_?_massfit.pdf'))
    plots = OrderedDict()
    for iplot, pname in enumerate(sorted(plot_files, key=get_bin_idx)):
        label = LABEL_BASE[bin_var].format(NICE_VARS[bin_var],
                                           *binning[iplot])
        plots[label] = pname

    return plots


def get_plot_idx_func(bin_var):
    """
    Get the function that returns the sort index for graph plots, such that they
    are matched against DESIRED_PLOT_ORDER
    """
    repl_string = '_v_{}.pdf'.format(bin_var)

    def get_plot_idx(plotname):
        """
        Closure over repl_string doing the actual work
        """
        try:
            pdfname = plotname.split('/')[-1]
            return DESIRED_PLOT_ORDER.index(pdfname.replace(repl_string, ''))
        except ValueError:
            pass
        return len(DESIRED_PLOT_ORDER)

    return get_plot_idx


def get_param_plot_names(indir, bin_var):
    """
    Get the list of param v costh graphs
    """
    param_files = glob.glob(path.join(indir, '*_v_{}.pdf'.format(bin_var)))
    plots = OrderedDict()
    for iplot, pname in enumerate(sorted(param_files,
                                         key=get_plot_idx_func(bin_var))):
        label = path.basename(pname).replace('_v_{}.pdf'.format(bin_var), '')
        if label in PLOT_LABELS_LATEX:
            label = PLOT_LABELS_LATEX[label]
        else:
            label = label.replace('_', r'\_')
        plots[label] = pname

    return plots


def create_fit_param_table(frfile, fit_config):
    """
    Write a table with all the parameters that are fitted or fixed
    """
    fitfile = r.TFile.Open(frfile)
    wsp = fitfile.Get('ws_mass_fit')

    ret_str = []

    ret_str.append(r'\begin{table}')
    ret_str.append(r'\begin{tabulary}{1.0\linewidth}{l c }')
    ret_str.append(r'parameter & value \\ \hline')
    for exp in fit_config['expression_strings']:
        par = exp.split('[')[0]
        if par in EXCLUDE_PARS or par.startswith('expr'):
            continue
        par_text = r'\texttt{{{}}}'.format(par).replace('_', r'\_')
        ret_str.append(r'{} & {} \\'.format(par_text, param_str(wsp, par)))

    ret_str.append(r'\end{tabulary}')
    ret_str.append(r'\end{table}')

    return '\n'.join(ret_str)


def main(args):
    """Main"""
    with open(args.fitconfig, 'r') as fit_conf_f:
        fit_config = json.load(fit_conf_f)

    bin_var, binning = get_binning(fit_config)
    fit_plots = get_fit_res_names(args.plotdir, bin_var, binning)

    graph_plots = get_param_plot_names(args.plotdir,
                                       bin_var.replace('_HX_fold', ''))

    with open_tex_file(args.outfile) as tex_file:
        tex_file.write(r'\paragraph{Free parameters}')
        tex_file.write(create_figure(graph_plots))

        tex_file.write('\n\\clearpage\n\\paragraph{Fitted / fixed parameters}\n')
        tex_file.write(create_fit_param_table(args.fitresfile, fit_config))

        tex_file.write('\n\\clearpage\n\\paragraph{Fit result plots}\n')
        tex_file.write(create_figure(fit_plots))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script that generates a tex '
                                     'stub that contains all the fit results')
    parser.add_argument('fitconfig', help='Fit configuration json file (for '
                        'binning information)')
    parser.add_argument('plotdir', help='Input directory that contains the '
                        'plots')
    parser.add_argument('fitresfile', help='File containing the fit results')
    parser.add_argument('-o', '--outfile', help='Ouptut file name',
                        default='fit_report.tex')

    clargs = parser.parse_args()
    main(clargs)
