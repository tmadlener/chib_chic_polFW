"""
Module containing helper functions for handling ROOT TH1Ds (or similar)
"""

import numpy as np

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname) - %(funcName)s: %(message)s')

from utils.misc_helpers import make_iterable

def draw_var_to_hist(tree, hist, var, cut='', weight=None):
    """
    Fill passed variable from TTree into TH1 using TTree.Draw().

    Args:
        tree (ROOT.TTree): A TTree that holds all the branches that are needed
            for drawing (including cuts and weight)
        hist (ROOT.TH1): A histogram that will be filled using the Draw()
            command.
        var (str): The name or the expression to draw into the histogram
        cut (str): A cut expression that will be applied (if not empty)
        weight (str): The name of the branch that should be used as weight for
            filling the histogram

    Returns:
        int: The return value of the TTree.Draw() command
    """
    plot_str = " >> ".join([var, hist.GetName()])
    if weight is not None:
        if cut:
            cut = "".join(["(", cut, ")"])
            cut = " * ".join([weight, cut])
        else:
            cut = weight

    return tree.Draw(plot_str, cut)


def combine_cuts(cuts, comb=' && '):
    """Concatenate list of cuts into a cut expression.

    Concatenates the list of cuts into an expression where all the cuts will be
    the AND combination of the individual cuts.

    Args:
        cuts (list of str): the individual cuts that should be combined into one
        comb (str): Logical expression to combine the cuts (e.g. && or ||)

    Returns:
        str: A cut expression string with all cuts combined via &&
    """
    safe_cuts = ("".join(["(", cut, ")"]) for cut in cuts)
    return comb.join(safe_cuts)


def set_hist_opts(hists):
    """Set some 'sane' default options to the passed histogram.

    Args:
        hists (ROOT.TH1 or list): the histogram(s) to set the options to.
    """
    hists = make_iterable(hists)
    for hist in hists:
        hist.SetStats(0) # disable stat box
        hist.Sumw2()


def set_bins_to_zero(hist, thresh=0):
    """Set bins under threshold to zero.

    Args:
        hist (ROOT.TH1): histogram for which bins should be set to zero.
        thresh (float, optional): Threshold below which bins should be set to 0.
        verbose (bool, optional): Make function print out the bins it sets to 0.
    """
    logging.debug('Checking {} for bins with entries below {}'
                  .format(hist.GetName(), thresh))
    neg_bins = [(i, b) for i, b in enumerate(hist) if b < thresh]
    for negb, cont in neg_bins:
        hist.SetBinContent(negb, 0)
        hist.SetBinError(negb, 0)
        logging.debug('Set bin {} to 0, content was {}'.format(negb, cont))


def set_labels(hist, xlabel='', ylabel=''):
    """
    Set the labels to the passed histogram.

    Args:
        hist (ROOT.TH1): Histogram for which the labels should be set
        xlabel, ylabel (str, optional): The x-, resp. y-label to set
    """
    if xlabel:
        hist.SetXTitle(xlabel)
    if ylabel:
        hist.SetYTitle(ylabel)


def get_y_max(hists):
    """
    Get the maximum y-value of all histograms

    Args:
        hists (list or ROOT.TH1): list of ROOT.TH1 for which the maximum y-value
            should be obtained

    Returns:
        float: The maximum y-value of all passed histograms
    """
    return max(h.GetBinContent(h.GetMaximumBin()) for h in make_iterable(hists))


def get_y_min(hists):
    """
    Get the minimum y-value of all histograms

    Args:
        hists (list or ROOT.TH1): list of ROOT.TH1 for which the minimum y-value
            should be obtained

    Returns:
        float: The minimum y-value of all passed histograms
    """
    return min(h.GetBinContent(h.GetMinimumBin()) for h in make_iterable(hists))


def get_x_max(hists):
    """
    Get the maximum x-value of all histograms

    Args:
        hists (list or ROOT.TH1): list of ROOT.TH1 for which the maximum x-value
            should be obtained

    Returns:
        float: The maximum x-value of all passed histograms
    """
    max_bin = lambda h: h.GetNbinsX()
    return max(h.GetBinLowEdge(max_bin(h)) + h.GetBinWidth(max_bin(h))
               for h in make_iterable(hists))


def get_x_min(hists):
    """
    Get the minimum x-value of all histograms

    Args:
        hists (list or ROOT.TH1): list of ROOT.TH1 for which the minimum x-value
            should be obtained

    Returns:
        float: The minimum x-value of all passed histograms
    """
    # 0 is underflow bin (thus starting at 1)
    return min(h.GetBinLowEdge(1) for h in make_iterable(hists))


def set_range_hist(hist, x_range=None, y_range=None):
    """
    Set the range to the histogram.

    Args:
        hist (ROOT.TH1): histogram for which the range should be set
        x_range, y_range (list or tuple, optional): x- resp. y-range to be used.
            Must be at least two numbers or None. If None this range will not be
            set for this histogram
    """
    if x_range is not None:
        hist.GetXaxis().SetRangeUser(x_range[0], x_range[1])
    if y_range is not None:
        hist.GetYaxis().SetRangeUser(y_range[0], y_range[1])


def set_range_hists(hists, x_range=None, y_range=None):
    """
    Set the range to all histograms

    Args:
        hists (list): ROOT.TH1s for which the range should be set
        x_range, y_range (list or tuple, optional): x- resp. y-range to be used.
            Must be at least two numbers or None. If None this range will not be
            set for the histograms
    """
    for hist in hists:
        set_range_hist(hist, x_range, y_range)


def set_common_range(hists, axis='xy', dscale=0.1, drange=[None, None]):
    """
    Set the same range to all histograms

    Determine the range such that all points from all histograms fit an set it
    to all histograms. It is also possible to pass in a default range that
    overrides the automatically determined ranges

    Args:
        hists (list): ROOT.TH1s for which the range should be unified
        axis (str, optional): Axis for which the range should be unified.
            Defaults to 'xy' so that the x and y axis will be set to the same
            range on all histograms
        dscale (float, optional): Scale factor to be applied to the minimum and
            maximum before setting the range.
        drange (list, optional): Default range. If any of the two values is
            not None, this value will be used instead of the automatically
            determined

    """
    if drange is None: # "normalize" None input
        drange = [None, None]

    # depending on the sign of the min and maximum values it is necessary to
    # either add or subtract from the values to have the widening in range go
    # into the right direction
    dmin = lambda val: dscale * val if val < 0 else -dscale * val
    dmax = lambda val: dscale * val if val > 0 else -dscale * val

    # define dict for mapping axis directions to corresponding functions for
    # obtaining extremal values
    rfuncs = {
        'y': {'max': get_y_max, 'min': get_y_min},
        'x': {'max': get_x_max, 'min': get_x_min}
    }
    ranges = {'x': None, 'y': None}

    # Define a lambda that gives either the default value (that is provided) as
    # an argument to this function or use the value determined from the
    # histograms. (Use a lambda instead of a full blown function)
    orideval = lambda val, dval, defv: val + dval(val) if defv is None else defv

    for axd in ['x', 'y']:
        if axd in axis:
            amin = rfuncs[axd]['min'](hists)
            amax = rfuncs[axd]['max'](hists)
            ranges[axd] = [orideval(amin, dmin, drange[0]),
                           orideval(amax, dmax, drange[1])]

    set_range_hists(hists, ranges['x'], ranges['y'])


def get_quantiles(hist, quantiles):
    """
    Get the (approximate) quantiles from the passed histogram

    Depending on the binning it is possible that the quantiles are not exactly
    at the desired values. In this case, the upper bound of the bin, that exceeds
    the quantile will be used.

    Args:
        hist (ROOT.TH1): histogram of distribution for which quantile cuts
            should be determined
        quantiles (list): quantile values between 0 and 1

    Returns:
        list: x-value at the upper bounds of the bins that exceed the quantiles
    """
    integral = hist.Integral()
    # calculate the cumulative sum of the bin entries, excluding under- and
    # overflow
    cum_sum = np.cumsum([hist[i] for i in xrange(1, hist.GetNbinsX() +1)])
    cum_sum /= integral # normalize

    # get the first indices in the normalized cumulative sum array that are
    # greater or equal to the desired quantiles
    q_bins = [np.where(cum_sum >= q)[0][0] for q in quantiles]

    x_axis = hist.GetXaxis()
    # when retrieving the upper bound adjust for the indexing in ROOT histograms
    # (starts at 1 and ends at Nbins for the "nominal" bins, 0 is underflow,
    # Nbins + 1 is overflow)
    return [x_axis.GetBinUpEdge(b + 1) for b in q_bins]
