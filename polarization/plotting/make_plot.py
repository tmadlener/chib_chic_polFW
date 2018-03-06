#!/usr/bin/env python
"""
Script for creating a plot with multiple inputs
"""

import re

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.plot_helpers import mkplot, default_colors
from utils.setup_plot_style import set_TDR_style, add_lumi_info
from utils.PlotServer import PlotServer
from utils.hist_utils import get_y_max
from utils.data_base import JsonDataBase
from utils.misc_helpers import get_full_trigger

# values at which chic1 and chic2 should be normalized to at 0 if that
# normalization option is chosen
chic1_norm_0 = 1.0
chic2_norm_0 = 0.5

# fit ranges for different variables
const_fit_range = {
    'costh': (0, 0.5),
    'phi': (0, 90), # full range
    'cosalpha': (0, 1) # full range
}

def get_lumi_text(trg_years):
    """
    Get the integrated lumi text for

    TODO: proper doc
    """
    lumi_db = JsonDataBase()
    lumi_text_frags = []
    for year, trg, dmc in trg_years:
        if dmc == 'mc': # NOTE: no lumi for MC at the moment
            continue
        trigger = get_full_trigger(trg)
        trg_lumi = lumi_db.get_int_lumi(year, trigger)
        energy = lumi_db.get_energy(year)
        logging.debug('Got int. lumi = {} and cms energy = {} for year {} and'
                      'trigger {}'.format(trg_lumi, energy, year, trigger))
        if trg_lumi > 0:
            lumi_text_frags.append('{0:.1f} fb^{{-1}} ({1} TeV)'
                                   .format(trg_lumi, energy))

    if len(lumi_text_frags) > 1:
        return ' + '.join(lumi_text_frags)
    return lumi_text_frags[0]


def get_plot_attributes(year, dmc, plot):
    """
    Get the plot attributes for a given plot
    """
    year_colors = {
        '2012': 1,
        '2016': default_colors()[0],
        '2017': default_colors()[1]
    }
    marker_style = {'2012': 20, '2016': 22, '2017': 23}
    marker_style_open = {'2012': 24, '2016': 26, '2017': 32}
    fill_style = {'chic1': 3354, 'chic2': 3345,
                  'chib1': 3354, 'chib2': 3345 }
    size = 2
    linewidth = 2

    color = year_colors[year]
    if plot == 'ratio' or plot == 'chic2' or plot == 'chib2':
        marker = marker_style[year]
    else:
        marker = marker_style_open[year]

    attr = {'color': color, 'marker': marker, 'size': size, 'width': linewidth}
    if dmc == 'mc':
        attr['size'] = 0 # do not plot markers for mc
        attr['width'] = 0
        if plot == 'ratio':
            attr['fillalpha'] = (color, 0.5)
        else:
            attr['fillstyle'] = fill_style[plot]

    return attr


def setup_legend(frame, ratio):
    """
    Setup the legend
    """
    position = {
        'CS': [0.2, 0.6, 0.65, 0.75],
        'HX': [0.5, 0.65, 0.95, 0.8],
        'PX': [0.6, 0.5, 0.9, 0.65] # TODO: only CS and HX "optimized" for now
    }
    # TODO: dynamic placement
    leg = r.TLegend(*position[frame])
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetEntrySeparation(0.01)
    leg.SetBorderSize(0)
    return leg


def normalize_histos(hchic1, hchic2, norm_mode):
    """
    Normalize the chic1 and chic2 histograms according to a certain
    normalization mode.

    Args:
        hchic1, hchic2 (ROOT.TH1): Histogram that are normalized assuming they
            belong together (i.e. are from the same year)
        norm_mode (str): Has to be one of the following:
            - 'nsignal': both histograms are normalized by the combined sum of
              their integrals
            - 'at0': both histograms are normalized according to their value at
              0. The chic1 histograms will be normalized to 1, the chic2
              histograms will be normalized to 0.5
    """
    logging.debug('Normalizing histograms {} and {} with mode {}'
                      .format(hchic1.GetName(), hchic2.GetName(), norm_mode))
    if norm_mode == 'nsignal':
        n_signal = hchic1.Integral() + hchic2.Integral()
        if n_signal <= 0:
            logging.warning('sum of integrals is 0 in normalizing histograms: '
                            '{} and {}'
                            .format(hchic1.GetName(), hchic2.GetName()))
        else:
            hchic1.Scale(1.0 / n_signal)
            hchic2.Scale(1.0 / n_signal)
    elif norm_mode == 'at0':
        get_bin_0 = lambda h: h.GetBinContent(h.FindBin(0))
        n_chic1 = get_bin_0(hchic1)
        if n_chic1 > 0:
            hchic1.Scale(chic1_norm_0 / n_chic1)
        n_chic2 = get_bin_0(hchic2)
        if n_chic2 > 0:
            hchic2.Scale(chic2_norm_0 / n_chic2)

        if n_chic1 <= 0 or n_chic2 <= 0:
            logging.warning('bin contents at 0 were {} and {} for histograms: '
                            '{} and {}'
                            .format(n_chic1, n_chic2,
                                    hchic1.GetName(), hchic2.GetName()))
    else:
        logging.warning('norm_mode was not a valid option. Did not normalize '
                         'histograms {} and {}'
                         .format(hchic1.GetName(), hchic2.GetName()))


def normalize_ratio(hratio, fit_range=(0, 0.5)):
    """ Normalize the ratio histogram by a constant obtained from a fit to a
    constant.

    TODO: proper doc
    """
    log_level = logging.getLogger().getEffectiveLevel()
    logging.debug('Fitting ratio hist {} to constant in range: [{}, {}]'
                  .format(hratio.GetName(), *fit_range))
    const_func = r.TF1('const_func', '[0]', *fit_range)
    fit_opt = 'SR0' if log_level == logging.DEBUG else 'SRq0'
    fit_rlt = hratio.Fit(const_func, fit_opt)

    if int(fit_rlt) != 0:
        logging.warning('Problem in fitting ratio histogram to constant')
        return -1, -1

    if log_level == logging.DEBUG:
        fit_rlt.Print()

    norm = fit_rlt.Parameter(0)
    hratio.Scale(1.0 / norm)

    return fit_rlt.Chi2(), fit_rlt.Ndf()


def get_plot_hists(pserver, trg_years, ratio, var, frame, pt, state='chic',
                   norm_mode='nsignal', norm_ratio=False):
    """
    Get all plots that are necessary for the desired plot

    Args:
        trg_years (list of tuples): each year needs an accompanying trigger
            and wether it is data or mc
    """
    hists = {}
    for year, trg, dmc in trg_years: # discard leg entries for now
        logging.debug('pt = {}, year = {}, dmc = {}, var = {}, frame = {}'
                              .format(pt, year, dmc, var, frame))
        if ratio:
            hist = pserver.get_hist(dmc, year, trg, pt, var, frame, 'ratio')
            if hist is not None:
                hists[(year, dmc, 'ratio')] = hist
                if norm_ratio:
                    chi2, ndf = normalize_ratio(hists[(year, dmc, 'ratio')],
                                                const_fit_range[var])
        else:
            hist = pserver.get_hist(dmc, year, trg, pt, var, frame, state + '1')
            hists[(year, dmc, state + '1')] = hist
            hist = pserver.get_hist(dmc, year, trg, pt, var, frame, state + '2')
            hists[(year, dmc, state + '2')] = hist

            normalize_histos(hists[(year, dmc, state + '1')],
                             hists[(year, dmc, state + '2')], norm_mode)

    return hists


def get_trigger_year_info(inputlist):
    """
    Parse the arguments from the --input list and return necessary info
    """
    logging.debug('Parsing input arguments: {}'.format(inputlist))
    ids = []
    leg_entries = {}

    for inp in inputlist:
        logging.debug('Parsing input: {}'.format(inp))
        args = inp.split(':')
        if len(args) == 3:
            args.append('') # append empty leg entry
        if len(args) >= 4:
            year, trg, data = args[:3]
            ids.append((year, trg, data),)
            leg_id = (year, data),
            leg_entries[leg_id] = args[3]
        else:
            logging.warning('Could not obtain all the necessary information '
                            'from {}. Skipping this input'.format(inp))

    return ids, leg_entries


def nice_leg_entries(leg_entries):
    """Replace dict strings with nicer plotting strings"""
    for i, entry in enumerate(leg_entries):
        leg_entries[i] = re.sub(r'chi(c|b)(\d)', r'#chi_{\1\2}', entry)

    return leg_entries


def get_xlabel(var, frame):
    """Get nice x-labels"""
    if var == 'phi':
        return '#phi^{{{0}}}_{{folded}}'.format(frame)
    x_label = var.replace('alpha', r'#alpha').replace('th', r'#theta')
    return '|{0}^{{{1}}}|'.format(x_label, frame)


def get_ylabel(dist, norm_at_zero, norm_ratio):
    """
    Get appropriate y-label
    """
    y_label = '#chi_{c2} / #chi_{c1}'
    if norm_ratio:
        y_label += ' (normalized by fit to const.)'
    if dist:
        if norm_at_zero:
            y_label = 'events (arb. units)'
        else:
            y_label = 'events / sum of signal #chi_{c1} and #chi_{c2} events'
    return y_label


def main(args):
    """Main"""
    # need PlotServer here, to keep the root file open
    pserver = PlotServer(args.histfile)

    norm_mode = 'at0' if args.normalize_at_zero else 'nsignal'

    year_trg_ids, leg_entries = get_trigger_year_info(args.input)
    hists = get_plot_hists(pserver, year_trg_ids, args.ratio, args.variable,
                           args.frame, args.ptbin, state=args.state,
                           norm_mode=norm_mode, norm_ratio=args.norm_ratio)


    # split into data and mc hist due to different plotting styles for the two
    data_hists = {k: hists[k] for k in hists if 'data' in k}
    mc_hists = {k: hists[k] for k in hists if 'mc' in k}

    leg = None
    data_leg = []
    mc_leg = []
    if args.legend:
        leg = setup_legend(args.frame, args.ratio )
        if args.ratio:
            data_leg = [leg_entries[(k[0], k[1]),] for k in data_hists]
            mc_leg = [leg_entries[(k[0], k[1]),] for k in mc_hists]
        else:
            data_leg = [', '.join([leg_entries[(k[0], k[1]),], k[2]])
                        for k in data_hists]
            data_leg = nice_leg_entries(data_leg)
            mc_leg = [', '.join([leg_entries[(k[0], k[1]),], k[2]])
                        for k in mc_hists]
            mc_leg = nice_leg_entries(mc_leg)


    dpl_attr = [get_plot_attributes(*k) for k in data_hists]
    mpl_attr = [get_plot_attributes(*k) for k in mc_hists]

    y_label = get_ylabel(not args.ratio, args.normalize_at_zero, args.norm_ratio)
    if args.state == 'chib':
        y_label = re.sub(r'c(\d)',r'b\1', y_label)

    x_label = get_xlabel(args.variable, args.frame)

    # to have same y-range for everything, have to do it manually, as
    # mkplot can't handle it in two different calls
    y_max = get_y_max(data_hists.values() + mc_hists.values()) * 1.1

    set_TDR_style()
    can = None # declare here to be able to switch the two calls below
    can = mkplot(mc_hists.values(), yRange = [0, y_max], drawOpt='E2',
                 attr=mpl_attr, can=can, leg=leg, legEntries=mc_leg,
                 legOpt='F', xLabel=x_label, yLabel=y_label)
    can = mkplot(data_hists.values(), yRange=[0, y_max], drawOpt='E1',
                 attr=dpl_attr, can=can, leg=leg, legEntries=data_leg,
                 legOpt='PLE', xLabel=x_label, yLabel=y_label)

    lumi_text = get_lumi_text(year_trg_ids)
    add_lumi_info(can, lumi_text)
    can.SaveAs(args.output)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for creating pub '
                                     'quality style plots')
    parser.add_argument('histfile', help='root files containing the histograms '
                        'in the structure as it is produced by the '
                        'collect_histograms script')
    parser.add_argument('-i', '--input', action='append',
                        help='Specify one given input in the format: '
                        'YYYY:TRG:DATAORMC:[:LEG], where the year, the trigger,'
                        ' data or mc and the (optional) legend entry are '
                        'separated by a colon')
    parser.add_argument('-o', '--output', default='pub_plot.pdf',
                        help='name of the output plot')
    parser.add_argument('-v', '--variable', type=str, default='costh',
                        help='variable to plot')
    parser.add_argument('-f', '--frame', type=str, default='HX',
                        help='reference frame')
    parser.add_argument('-pt', '--ptbin', type=int, default=0,
                        help='pt bin to plot')
    parser.add_argument('-l', '--legend', action='store_true', default=False,
                        help='put legend onto the plot')
    parser.add_argument('-nz', '--normalize-at-zero', default=False,
                        action='store_true', help='Normalize distribution plots'
                        ' at their value at 0. Default is to normalize by the '
                        'sum of the chic1 and chic2 Integrals')
    parser.add_argument('-nr', '--norm-ratio', help='Normalize the ratio by a '
                        'factor obtained from a fit to a constant',
                        action='store_true', default=False)


    plot_type = parser.add_mutually_exclusive_group(required=True)
    plot_type.add_argument('-r', '--ratio', action='store_true',
                           help='Make ratio plot')
    plot_type.add_argument('-d', '--dist', action='store_true',
                           help='Make dist plot')


    state_sel = parser.add_mutually_exclusive_group()
    state_sel.add_argument('--chic', action='store_const', dest='state',
                           const='chic', help='Do the plots for the chic')
    state_sel.add_argument('--chib', action='store_const', dest='state',
                           const='chib', help='Do the plots for the chib')
    parser.set_defaults(state='chic')



    clargs = parser.parse_args()
    main(clargs)
