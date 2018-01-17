"""
Module containing helper functions for handling ROOT TH1Ds (or similar)
"""

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

