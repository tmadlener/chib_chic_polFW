#!/usr/bin/env python
"""
Script to generate a tex stub that can be included in a .tex file and gives
a "reasonable" appearance
"""

from __future__ import print_function

import glob
import json
import sys
import re
from os import path, environ

from collections import OrderedDict

from utils.reporting import open_tex_file
from utils.plot_decoration import PLOT_LABELS_LATEX

NICE_VARS = {'phi_HX_fold': r'\varphi^{HX}_{\text{fold}}', 'costh_HX_fold':
             r'|\cos\vartheta^{HX}|'}
LABEL_BASE = {'phi_HX_fold': r'${1:.0f} < {0} < {2:.0f}$',
              'costh_HX_fold': r'${1:.3f} < {0} < {2:.3f}$'}

# maximum number of figures per page
MAX_FIG_P_PAGE = 6
SORT_RGX = re.compile(r'bin_([0-9]+)\.pdf')

try:
    LATEX_EXE = environ['LATEX_EXE']
except KeyError:
    LATEX_EXE = 'pdflatex'
    pass

PLOT_COMMAND = r'\includegraphics[width=0.5\linewidth]'
if LATEX_EXE == 'xelatex':
    PLOT_COMMAND = r'\rootfig{0.475}'

# Define an order such that the similar plots always appear next to each other
DESIRED_PLOT_ORDER = [
    'CBmass1', 'CBmass2',
    'CBsigma1', 'CBsigma2',
    'mu_bkg', 'sigma_bkg',
    'lambda_bkg', 'Nbkg',
    'Nchic1', 'Nchic2',
    'Nchic0',
]


def get_bin_idx(plotname):
    """
    Get the bin index from the plot name
    """
    return int(SORT_RGX.search(plotname).group(1))


def get_binning(bin_info_file):
    """
    Read the bin-info json and get the bin variable and the binning
    """
    with open(bin_info_file, 'r') as bin_info_f:
        bin_info = json.load(bin_info_f)

    return bin_info['bin_variable'], bin_info['costh_bins']


def get_fit_res_names(indir, bin_var, binning):
    """
    Get the list of fit_res files
    """
    plot_files = glob.glob(path.join(indir, 'mass_fit_config_costh_bin_?.pdf'))
    plots = OrderedDict()
    for iplot, pname in enumerate(sorted(plot_files, key=get_bin_idx)):
        label = LABEL_BASE[bin_var].format(NICE_VARS[bin_var],
                                           *binning[iplot])
        plots[label] = pname

    return plots


def get_plot_idx(plotname):
    """
    Get the sort index for the graph plots, such that they are matched against
    the DESIRED_PLOT_ORDER
    """
    try:
        pdfname = plotname.split('/')[-1]
        return DESIRED_PLOT_ORDER.index(pdfname.replace('_v_costh.pdf', ''))
    except ValueError:
        pass
    return len(DESIRED_PLOT_ORDER)


def get_param_plot_names(indir):
    """
    Get the list of param v costh graphs
    """
    param_files = glob.glob(path.join(indir, '*_v_costh.pdf'))
    plots = OrderedDict()
    for iplot, pname in enumerate(sorted(param_files, key=get_plot_idx)):
        label = path.basename(pname).replace('_v_costh.pdf', '')
        if label in PLOT_LABELS_LATEX:
            label = PLOT_LABELS_LATEX[label]
        else:
            label = label.replace('_', r'\_')
        plots[label] = pname

    return plots


def create_subfloat(plot, label):
    """
    Create a subfloat entry
    """
    return (r'\subfloat[][LABEL]{PLTCMD{PLOT}}'
            .replace('LABEL', label).replace('PLOT', plot)
            .replace('PLTCMD', PLOT_COMMAND))


def create_figure(plots):
    """
    Make the latex figure
    """
    ret_str = []

    fig_started = False
    for iplot, (label, pname) in enumerate(plots.iteritems()):
        if iplot % MAX_FIG_P_PAGE == 0:
            if fig_started:
                ret_str.append(r'\end{figure}')
                ret_str.append('')
                fig_started = False
            if not fig_started:
                ret_str.append(r'\begin{figure}')
                ret_str.append('')
                fig_started = True

        ret_str.append(create_subfloat(pname, label))

        if iplot % 2 != 0:
            ret_str.append('')

    if fig_started:
        ret_str.append(r'\end{figure}')

    return '\n'.join(ret_str)


def main(args):
    """Main"""
    bin_var, binning = get_binning(args.bininfofile)
    fit_plots = get_fit_res_names(args.inputdir, bin_var, binning)

    graph_plots = get_param_plot_names(args.inputdir)

    with open_tex_file(args.outfile) as tex_file:
        tex_file.write(r'\paragraph{Free parameters}')
        tex_file.write(create_figure(graph_plots))

        tex_file.write('\n\\clearpage\n\\paragraph{Fit result plots}\n')
        tex_file.write(create_figure(fit_plots))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script that generates a tex '
                                     'stub that contains all the fit results')
    parser.add_argument('inputdir', help='Input directory that contains the '
                        'plots')
    parser.add_argument('bininfofile', help='Bin info json file created by fit '
                        'script')
    parser.add_argument('-o', '--outfile', help='Ouptut file name',
                        default='fit_report.tex')

    clargs = parser.parse_args()
    main(clargs)
