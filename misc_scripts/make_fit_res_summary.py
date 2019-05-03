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
from os import path


NICE_VARS = {'phi_HX_fold': r'\varphi^{HX}_{\text{fold}}', 'costh_HX_fold':
             r'|\cos\vartheta^{HX}|'}
LABEL_BASE = {'phi_HX_fold': r'${1:.0f} < {0} < {2:.0f}$',
              'costh_HX_fold': r'${1:.3f} < {0} < {2:.3f}$'}

# maximum number of figures per page
MAX_FIG_P_PAGE = 6
SORT_RGX = re.compile(r'bin_([0-9]+)\.pdf')

def get_bin_idx(plotname):
    """
    Get the bin index from the plot name
    """
    return int(SORT_RGX.search(plotname).group(1))


def get_binning(indir):
    """
    Read the bin-info json and get the bin variable and the binning
    """
    json_files = glob.glob(path.join(indir, '*.json'))
    json_files.remove(path.join(indir, 'fit_model.json'))

    if len(json_files) > 1:
        print('Have more than 1 json file: ')
        for ifile, filen in enumerate(json_files):
            print('[{}] {}'.format(ifile, filen))
        while 1:
            ifile = raw_input('Chose one of the above files (q to quit): ')
            if ifile == 'q':
                sys.exit(1)
            if ifile < len(json_files):
                bin_info_file = json_files[ifile]
    else:
        bin_info_file = json_files[0]

    with open(bin_info_file, 'r') as bin_info_f:
        bin_info = json.load(bin_info_f)

    return bin_info['bin_variable'], bin_info['costh_bins']


def get_fit_res_names(indir):
    """
    Get the list of fit_res files
    """
    return glob.glob(path.join(indir, 'mass_fit_config_costh_bin_?.pdf'))


def create_subfloat(plot, label):
    """
    Create a subfloat entry
    """
    # print(r'\subfloat[][LABEL]{\includegraphics[width=0.5\linewidth]{PLOT}}'
    print(r'\subfloat[][LABEL]{\rootfig{0.45}{PLOT}}'
          .replace('LABEL', label).replace('PLOT', plot))


def create_figure(plots, bin_var, binning):
    """
    Make the latex figure
    """
    plots = sorted(plots, key=get_bin_idx)

    fig_started = False
    for iplot, pname in enumerate(plots):
        if iplot % 6 == 0:
            if fig_started:
                print(r'\end{figure}')
                print('')
                fig_started = False
            if not fig_started:
                print(r'\begin{figure}')
                fig_started = True

        label_txt = LABEL_BASE[bin_var].format(NICE_VARS[bin_var],
                                               *binning[iplot])
        create_subfloat(pname, label_txt)
        if iplot % 2 != 0:
            print('')

    if fig_started:
        print(r'\end{figure}')


def main(args):
    """Main"""
    bin_var, binning = get_binning(args.inputdir)
    plots = get_fit_res_names(args.inputdir)
    create_figure(plots, bin_var, binning)




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script that generates a tex '
                                     'stub that contains all the fit results')
    parser.add_argument('inputdir', help='Input directory that contains the '
                        'plots (and also the bin info json)')

    clargs = parser.parse_args()
    main(clargs)
