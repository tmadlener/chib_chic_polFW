#!/USSR/bin/env python
"""
Make the plots that are planned to go into the publication
"""
import os
# import logging
# logging.basicConfig(level=logging.DEBUG,
#                    format='%(levelname)s - %(funcName)s: %(message)s')

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.hist_utils import (
    draw_var_to_hist, set_hist_opts, set_bins_to_zero, set_range_hist,
    get_y_max
)
from utils.setup_plot_style import set_TDR_style, add_auxiliary_info
from utils.plot_helpers import default_colors
#################################################################
# some global variables, that should at one point be changed to #
# input parameters probably                                     #
#################################################################
# base directory where all the data is and the plots will be stored
BASE_DIR = '/'.join([os.environ['WORK'], 'NewPolMethodTests', 'data', 'Chic'])
# input data files, that will all be put onto one plot
# IINPUT_FILES = {
#     '2012': '/'.join([BASE_DIR, 'tuples', 'chic_tuple_ptmerged_2D.root']),
#     '2016': '/'.join([BASE_DIR, 'tuples', 'data_2016', 'chic_tuple_ptmerged_2D.root']),
#     '2017': '/'.join([BASE_DIR, 'tuples', 'data_2017', 'chic_tuple_pt3_2D.root'])
# }
# INPUT_FILES = {
#     '2012': '/'.join([BASE_DIR, 'tuples', 'data_2012', 'cowboy_rejection_after_fit', 'chic_tuple_pt2_2D.root']),
#     '2017': '/'.join([BASE_DIR, 'tuples', 'data_2017', 'chic_tuple_pt3_2D.root'])
# }

INPUT_FILES = {
    '2017': '../foo_test_creation_script/chic_tuple_pt3.root'
}

# entries to the legend
LEG_ENTRIES = {
    '2012': '8 TeV (2012)',
    '2016': '13 TeV (2016)',
    '2017': '13 TeV (2017)'
}
# reference frames that should be plotted
PLOT_FRAMES = ['HX', 'CS', 'PX']
# output base directory
OUTPUT_BASE_DIR = '/'.join([BASE_DIR, 'plot_tests'])
# name of the Tyree in the datafile
TREE_NAME = 'chic_tuple'
# names of the weights
WEIGHTS = {'chic1': 'wChic1', 'chic2': 'wChic2'}
# number of bins to use for the histograms
N_BINS = {'phi': 20, 'costh': 16}

##############################################################################
# global settings for plotting (that override the settings in set_TDR_style) #
##############################################################################
CANVAS_SIZE = [600, 600]
MARKER_STYLE = [20, 22, 23, 21, 33]
MARKER_STYLE_OPEN = [24, 26, 32, 25, 27]
MARKER_SIZE = 2
LINE_WIDTH = 2
COLOR = {
    '2012': 1,
    '2016': default_colors()[0],
    '2017': default_colors()[1]
}

LEG_TEXT_FONT = 42
LEG_TEXT_SIZE = 0.04
LEG_SEPARATION = 0.01 # entry separation in legend

## settings for ratio plots
RATIO_Y_RANGE = {'chic1': None, 'chic2': None, 'ratio': [0, 1]}
RATIO_X_RANGE = {'phi': [0,90], 'costh': [0, 1]}
RATIO_GRID = False
RATIO_Y_LABEL = '#chi_{c2} / #chi_{c1}'
RATIO_LEG_POS = [0.25, 0.15, 0.45, 0.30] # TODO: this has to be dynamically done depending on the number of points, etc.

## some settings for dist plots
DIST_Y_LABEL = 'events / sum of signal #chi_{c1} and #chi_{c2} events'
DIST_LEG_POS = { # TODO; this has to be done dynamically at some point
    'CS': [0.2, 0.6, 0.65, 0.75],
    'HX': [0.5, 0.65, 0.95, 0.8],
    'PX': [0.6, 0.5, 0.9, 0.65] # TODO: only CS and HX "optimized" for now
}

####################################################################
# global dict of open TFiles to avoid letting them go out of scope #
# using a dict to avoid opening the same file multiple times       #
####################################################################
_OPEN_FILES = {}


def get_dist_hists(datafile, frame):
    """
    Create the folded phi and costh hists from the passed datafile.
    """
    if not datafile in _OPEN_FILES:
        f = r.TFile.Open(datafile)
        _OPEN_FILES[datafile] = f
    else:
        f = _OPEN_FILES[datafile]

    t = f.Get(TREE_NAME)
    draw_names = {'costh': 'TMath::Abs(' + '_'.join(['costh', frame]) + ')',
                    'phi': '_'.join(['phi', frame, 'fold'])}

    # maybe need another argument to this function to avoid overwriting histograms
    hists = {
        'phi': {
            'chic1': r.TH1D('_'.join(['h', 'phi', frame, 'chic1']),
                            ''.join([';#phi^{', frame,'}_{folded}']), N_BINS['phi'], 0, 90),
            'chic2': r.TH1D('_'.join(['h', 'phi', frame, 'chic2']),
                            ''.join([';#phi^{', frame,'}_{folded}']), N_BINS['phi'], 0, 90)
        },
       'costh': {
            'chic1': r.TH1D('_'.join(['h', 'costh', frame, 'chic1']),
                            ''.join([';|cos#theta^{', frame,'}|']), N_BINS['costh'], 0, 1),
            'chic2': r.TH1D('_'.join(['h', 'costh', frame, 'chic2']),
                            ''.join([';|cos#theta^{', frame,'}|']), N_BINS['costh'], 0, 1)
        }
    }

    for state in ['chic1', 'chic2']:
        for var in ['costh', 'phi']:
            set_hist_opts(hists[var][state])
            draw_var_to_hist(t, hists[var][state], draw_names[var], '', WEIGHTS[state])

    create_ratio_hists(hists)

    return hists


def create_ratio_hists(hists, sanitize=True):
    """Create the ratios of the phi and put them into the passed dictionary"""
    if sanitize:
        for var in hists:
            for state in hists[var]:
                set_bins_to_zero(hists[var][state], 0, True)

    for var in ['costh', 'phi']:
        hists[var]['ratio'] = hists[var]['chic2'].Clone(hists[var]['chic2'].GetName().replace('chic2', 'ratio'))
        hists[var]['ratio'].Divide(hists[var]['chic1'])


def collect_hists(input_files, frames):
    """
    Collect all histograms from all files and all frames
    """
    histograms = {}
    for year in input_files:
        histograms[year] = {}
        for frame in frames:
            histograms[year][frame] = get_dist_hists(INPUT_FILES[year], frame)

    return histograms


def get_plot_name(plot, var, years, frame):
    """Get a unique canvas name that can also be used as savename"""
    return '_'.join([plot, var, frame] + years)


def get_hists_to_plot(hists, var, years, frame, plot):
    """
    Get a list of histograms to plot from the dictionary
    """
    p_hists = []
    for year in years:
        p_hists.append(hists[year][frame][var][plot])
        p_hists[-1].SetName(p_hists[-1].GetName() + '_' + year)

    return p_hists

def setup_legend(position):
    """Setup the legend"""
    leg = r.TLegend(*position)
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(LEG_TEXT_FONT)
    leg.SetTextSize(LEG_TEXT_SIZE)
    leg.SetBorderSize(0)
    leg.SetEntrySeparation(LEG_SEPARATION)

    return leg


def make_ratio_plot(hists, var, years, frame):
    """
    Make ratio plots
    """
    set_TDR_style()
    plot_name = get_plot_name('ratio', var, years, frame)
    can = r.TCanvas(plot_name, '', 50, 50, CANVAS_SIZE[0], CANVAS_SIZE[1])
    can.cd()
    if RATIO_GRID:
        can.SetGrid()

    leg = setup_legend(RATIO_LEG_POS)

    plot_hists = get_hists_to_plot(hists, var, years, frame, 'ratio')
    for i, h in enumerate(plot_hists):
        h.SetLineWidth(LINE_WIDTH)
        h.SetLineColor(COLOR[years[i]])
        h.SetMarkerSize(MARKER_SIZE)
        h.SetMarkerStyle(MARKER_STYLE[i])
        h.SetMarkerColor(COLOR[years[i]])
        h.SetYTitle(RATIO_Y_LABEL)
        set_range_hist(h, RATIO_X_RANGE[var], RATIO_Y_RANGE['ratio'])

        h.Draw('E1 same')

        leg.AddEntry(h, LEG_ENTRIES[years[i]], 'ple')


    leg.Draw()
    add_auxiliary_info(can, years)
    can.Update()
    can.SaveAs('.'.join([plot_name, 'pdf']))

def get_aux_info_pos(frame, var, plot):
    """Get the position of the auxiliary info"""
    if frame == 'CS' and plot == 'dist':
        return 'left'

    return 'right'


def make_dist_plot(hists, var, years, frame):
    """Make chic distribution plots"""
    set_TDR_style()
    plot_name = get_plot_name('dist', var, years, frame)
    can = r.TCanvas(plot_name, '', 50, 50, *CANVAS_SIZE)
    can.cd()

    chic1_hists = get_hists_to_plot(hists, var, years, frame, 'chic1')
    chic2_hists = get_hists_to_plot(hists, var, years, frame, 'chic2')

    def normalize_hists(h_chic1, h_chic2):
        """
        Normalize the corresponding histograms of each year to the
        sum of signal events in the histograms

        NOTE: assuming that the indices correspond to the same year
        in both passed lists
        """
        for i,h in enumerate(h_chic1):
            n_sig_events = h.Integral() + h_chic2[i].Integral()
            h.Scale(1 / n_sig_events)
            h_chic2[i].Scale(1 / n_sig_events)

    class Legend():
        """Small wrapper around TLegend to avoid some if-else\'s in this function"""
        def __init__(self, var, frame):
            self.var = var
            if self.var == 'costh':
                self.leg = setup_legend(DIST_LEG_POS[frame])
                self.leg.SetNColumns(2)

        def AddEntry(self, *args):
            if self.var == 'phi':
                return
            self.leg.AddEntry(*args)

        def Draw(self, *args):
            if self.var == 'phi':
                return
            self.leg.Draw(*args)


    normalize_hists(chic1_hists, chic2_hists)
    y_max = get_y_max(chic1_hists + chic2_hists) * 1.1

    leg = Legend(var, frame)

    for i, h in enumerate(chic1_hists):
        h2 = chic2_hists[i]
        h.SetLineWidth(LINE_WIDTH)
        h.SetLineColor(COLOR[years[i]])
        h.SetMarkerSize(MARKER_SIZE)
        h.SetMarkerColor(COLOR[years[i]])
        h.SetMarkerStyle(MARKER_STYLE_OPEN[i])
        h.SetYTitle(DIST_Y_LABEL)
        set_range_hist(h, RATIO_X_RANGE[var], [0, y_max])

        h2.SetLineWidth(LINE_WIDTH)
        h2.SetLineColor(COLOR[years[i]])
        h2.SetMarkerSize(MARKER_SIZE)
        h2.SetMarkerColor(COLOR[years[i]])
        h2.SetMarkerStyle(MARKER_STYLE[i])
        h2.SetYTitle(DIST_Y_LABEL)
        set_range_hist(h, RATIO_X_RANGE[var], [0, y_max])
        leg.AddEntry(h2, '#chi_{c2}', 'ple')
        leg.AddEntry(h, '#chi_{c1} ' + LEG_ENTRIES[years[i]], 'ple')

        h.Draw('E1 same')
        h2.Draw('E1 same')

    leg.Draw()
    add_auxiliary_info(can, years, get_aux_info_pos(frame, var, 'dist'))
    can.Update()
    can.SaveAs('.'.join([plot_name, 'pdf']))


def main():
    """Main"""
    histograms = collect_hists(INPUT_FILES, PLOT_FRAMES)
    # TODO: make plots

    #years_to_plot = ['2012', '2016']
    years_to_plot = ['2017']
    # years_to_plot = ['2012', '2017']

    make_ratio_plot(histograms, 'costh', years_to_plot, 'HX')
    make_ratio_plot(histograms, 'phi', years_to_plot, 'HX')
    make_ratio_plot(histograms, 'costh', years_to_plot, 'CS')
    make_ratio_plot(histograms, 'phi', years_to_plot, 'CS')

    make_dist_plot(histograms, 'costh', years_to_plot, 'HX')
    make_dist_plot(histograms, 'phi', years_to_plot, 'HX')
    make_dist_plot(histograms, 'costh', years_to_plot, 'CS')
    make_dist_plot(histograms, 'phi', years_to_plot, 'CS')


if __name__ == '__main__':
    main()
