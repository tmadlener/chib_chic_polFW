"""
Module containing helper functions for handling ROOT TH1Ds (or similar)
"""

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname) - %(funcName)s: %(message)s')

from root_numpy import fill_hist

from utils.misc_helpers import (
    make_iterable, create_random_str, get_vals_from_rwbuffer
)

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


def _get_y_max_hist(hists):
    """
    Get the maximum y-value of all histograms

    Args:
        hists (list or ROOT.TH1): list of ROOT.TH1 for which the maximum y-value
            should be obtained

    Returns:
        float: The maximum y-value of all passed histograms
    """
    return max(h.GetBinContent(h.GetMaximumBin()) for h in make_iterable(hists))


def _get_y_min_hist(hists):
    """
    Get the minimum y-value of all histograms

    Args:
        hists (list or ROOT.TH1): list of ROOT.TH1 for which the minimum y-value
            should be obtained

    Returns:
        float: The minimum y-value of all passed histograms
    """
    return min(h.GetBinContent(h.GetMinimumBin()) for h in make_iterable(hists))


def _get_x_max_hist(hists):
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


def _get_x_min_hist(hists):
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


def get_binning(hist, axis='X'):
    """
    Get the binning of the passed histogram along the passed axis

    Args:
        hist (ROOT.TH1 or inheriting): Histogram for which the binning is of
            interest
        axis (char): Either 'X', 'Y' or 'Z' (depending on the passed hist not
            all are possible!). Defines for which axis the binning is desired

    Returns:
        np.array: The array with the bin eges
    """
    n_bins = getattr(hist, 'GetNbins' + axis)()
    # first check if the histogram has non-uniform binning
    # TODO: Check if there is a direct way to check if the edges buffer is null
    # instead of going through the try-except machinery
    edges = getattr(hist, 'Get' + axis + 'axis')().GetXbins().GetArray()
    try:
        return get_vals_from_rwbuffer(edges, n_bins + 1)
    except IndexError:
        pass
    return np.linspace(get_x_min(hist), get_x_max(hist), n_bins + 1)


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
        'y': {'max': _get_y_max_hist, 'min': _get_y_min_hist},
        'x': {'max': _get_x_max_hist, 'min': _get_x_min_hist}
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


def divide(num, denom, **kwargs):
    """
    Divide to histograms and return the ratio

    Args:
        num (ROOT.TH1): numerator histogram
        denom (ROOT.TH1): denominator histogram

    Keyword Args:
        name (str): If not empty this will be set as the name of the ratio TH1
        [x|y]label (str): Labels for the x and y-axis.

    Returns:
        ratio (ROOT.TH1): Ratio histogram obtained from cloning num and then
            dividing it by denom
    """
    ratio = num.Clone(kwargs.pop('name', ''))
    ratio.Divide(denom)

    set_labels(ratio, kwargs.pop('xlabel', ''), kwargs.pop('ylabel', ''))

    return ratio


def create_histogram(var, hist_sett, **kwargs):
    """
    Create a ROOT histogram from the passed variable(s)

    Args:
        var (np.array): Array with maximum of 3 columns containing the variables
            to plot.
        hist_set (tuple): Histogram settings, that are directly unpacked into
            the constructor of the ROOT histogram

    Keyword Args:
        name (str, optional): Name to be used for the histogram
        weights (np.array, optional): weight array with the same number of
             events as the var array. Each entry corresponds to the weight of
             the event
        {x,y,z}_axis (str): axis labels to be set for the histogram

    Returns:
         ROOT.TH{1,2,3}D: The histogram with the dimension corresponding to the
             number of columns of var
    """
    name = kwargs.pop('name', '')
    if not name:
        name = create_random_str()
    # use the number of dimensions from the var to determine which sort of
    # histogram to use
    ndim = var.shape
    if len(ndim) == 1:
        ndim = 1
    else:
        ndim = ndim[1]

    if ndim > 3 or ndim < 0:
        logging.error('Dimension of histogram is {}. Cannot create histogram'
                      .format(ndim))
        raise TypeError('Invalid number of dimensions in create_histograms')

    hist_type = 'TH{}D'.format(ndim)
    try:
        hist = getattr(r, hist_type)(name, '', *hist_sett)
    except TypeError as exc:
        logging.error('Could not construct TH{}D with passed hist_sett: {}'
                      .format(ndim, hist_sett))
        raise exc

    set_hist_opts(hist)

    # set axis labels
    xax, yax, zax = (kwargs.pop(a, '') for a in ['x_axis', 'y_axis', 'z_axis'])
    if xax:
        hist.SetXTitle(xax)
    if yax:
        hist.SetYTitle(yax)
    if zax:
        hist.SetZTitle(zax)

    fill_hist(hist, var, weights=kwargs.pop('weights', None))

    return hist


def _get_hist_sett(var, nbins=None, minx=None, maxx=None, hist_sett=None):
    """
    Get the "optimal" histogram settings for the passed variable or use
    some provided values.

    Returns:
        tuple: Tuple that can be directly unpacked in the constructor of TH1
    """
    if hist_sett is not None:
        return hist_sett

    if nbins is None:
        nbins = 100
    if minx is None:
        minx = np.min(var)
    if maxx is None:
        maxx = np.max(var)

    return (nbins, minx, maxx)


def hist1d(var, **kwargs):
    """
    Create a TH1D from the passed var.

    This is a convenience wrapper around create_histogram that does some
    automatic determination of the "optimal" histogram settings.

    Args:
        var (np.array): 1d array containing the variable to plot

    Keyword Args:
        hist_sett (tuple, optional): Histogram settings that are directly
            unpacked in the constructor of the ROOT histogram.
            NOTE: this overrides the nbins and min and max settings
        nbins (int, optional): Number of bins to use
        min (float, optional): Lower bound of histogram
        max (float, optional): Upper bound of histogram

    Returns:
        ROOT.TH1D: The histogram of the passed variable

    See also:
        create_histogram
    """
    hist_sett = _get_hist_sett(var, kwargs.pop('nbins', None),
                               kwargs.pop('min', None), kwargs.pop('max', None),
                               kwargs.pop('hist_sett', None))

    # use the name of the variable if it has one and nothing else is set
    if kwargs.get('x_axis', None) is None and hasattr(var, 'name'):
        kwargs['x_axis'] = var.name

    return create_histogram(var, hist_sett, **kwargs)


def hist2d(varx, vary, **kwargs):
    """
    Create a TH2D from the passed varx and vary

    This is a convenience wrapper around create_histogram that does some
    automatic determination of the "optimal" histogram settings

    Keyword Args:


    Returns:
        ROOT.TH2D: The 2d histogram of the passed variables

    See also:
        create_histogram
    """
    hist_sett = kwargs.pop('hist_sett', None)
    if hist_sett is None:
        x_sett = _get_hist_sett(varx, kwargs.pop('nbinsx', None),
                                kwargs.pop('minx', None),
                                kwargs.pop('maxx', None),
                                kwargs.pop('x_hist_sett', None))
        y_sett = _get_hist_sett(vary, kwargs.pop('nbinsy', None),
                                kwargs.pop('miny', None),
                                kwargs.pop('maxy', None),
                                kwargs.pop('y_hist_sett', None))
        hist_sett = x_sett + y_sett

    # use the name of the variables if they have one and nothing else is set
    if kwargs.get('x_axis', None) is None and hasattr(varx, 'name'):
        kwargs['x_axis'] = varx.name
    if kwargs.get('y_axis', None) is None and hasattr(vary, 'name'):
        kwargs['y_axis'] = vary.name

    return create_histogram(np.array([varx, vary]).T, hist_sett, **kwargs)
