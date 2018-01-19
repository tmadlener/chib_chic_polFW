"""
Module containing helper functions for handling ROOT TH1Ds (or similar)
"""

import numpy as np

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


def combine_cuts(cuts):
    """Concatenate list of cuts into a cut expression.

    Concatenates the list of cuts into an expression where all the cuts will be
    the AND combination of the individual cuts.

    Args:
        cuts (list of str): the individual cuts that should be combined into one

    Returns:
        str: A cut expression string with all cuts combined via &&
    """
    safe_cuts = ("".join(["(", cut, ")"]) for cut in cuts)
    return " && ".join(safe_cuts)


def set_hist_opts(hist):
    """Set some 'sane' default options to the passed histogram.

    Args:
        hist (ROOT.TH1): the histogram to set the options to.
    """
    hist.SetStats(0) # disable stat box
    hist.Sumw2()


def set_bins_to_zero(hist, thresh=0, verbose=False):
    """Set bins under threshold to zero.

    Args:
        hist (ROOT.TH1): histogram for which bins should be set to zero.
        thresh (float, optional): Threshold below which bins should be set to 0.
        verbose (bool, optional): Make function print out the bins it sets to 0.
    """
    neg_bins = [(i, b) for i, b in enumerate(hist) if b < thresh]
    for negb, cont in neg_bins:
        hist.SetBinContent(negb, 0)
        hist.SetBinError(negb, 0)

    if verbose:
        print('checked {} for bins with entries below {}'.
              format(hist.GetName(), thresh))
        for negb, cont in neg_bins:
            print('Set bin {} to 0, content was {}'.format(negb, cont))


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
        hists (list): list of ROOT.TH1 for which the maximum y-value should be
            obtained

    Returns:
        float: The maximum y-value of all passed histograms
    """
    return max([h.GetBinContent(h.GetMaximumBin()) for h in hists])

def get_x_max(hists):
    """
    Get the maximum x-value of all histograms

    Args:
        hists (list): list of ROOT.TH1 for which the maximum x-value should be
            obtained

    Returns:
        float: The maximum x-value of all passed histograms
    """
    max_bin = lambda h: h.GetNbinsX()
    return max([h.GetBinLowEdge(max_bin(h)) + h.GetBinWidth(max_bin(h))
                for h in hists])


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


def set_common_range(hists, axis='xy', dscale=0.1):
    """
    Set the same range to all histograms

    Determine the the range such that all points from all histograms fit an set
    it to all histograms.

    Args:
        hists (list): ROOT.TH1s for which the range should be unified
        axis (str, optional): Axis for which the range should be unified.
            Defaults to 'xy' so that the x and y axis will be set to the same
            range on all histograms
        dscale(float, optional): Scale factor to be applied to the minimum and
            maximum before setting the range.
    """
    x_range = None
    y_range = None

    # depending on the sign of the min and maximum values it is necessary to
    # either add or subtract from the values to have the widening in range go
    # into the right direction
    dmin = lambda val: dscale * val if val < 0 else -dscale * val
    dmax = lambda val: dscale * val if val > 0 else -dscale * val

    if 'y' in axis.lower():
        y_min = min([h.GetBinContent(h.GetMinimumBin()) for h in hists])
        y_max = get_y_max(hists)
        y_range = [y_min + dmin(y_min), y_max + dmax(y_max)]

    if 'x' in axis.lower():
        x_min = min([h.GetBinLowEdge(1) for h in hists]) # 0 is underflow bin
        x_max = get_x_max(hists)
        x_range = [x_min + dmin(x_min), x_max + dmax(x_max)]

    set_range_hists(hists, x_range, y_range)


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
