"""
Module containing helper functions for plotting and related stuff

Attributes:
    _colors (list): list of index and TColor for caching default colors
    _color_indices (list): list of indices for the different default colors,
        which will be used in ROOTs SetColor methods

    Neither of the two above variables should be accessed directly, but the
    default_colors method should be used
"""

import ROOT as r

from collections import Iterable
from utils.misc_helpers import create_random_str
from utils.hist_utils import set_common_range

_colors = []
_color_indices = []

def default_colors():
    """Get a list of some nicer colors for plotting than ROOTs default.

    Creates the colors that correspond to MATLABs default palette so that they
    can be used with ROOT as well. On the first call the colors will be
    initialized and stored in module variables. Subsequent calls will only
    retrieve the cached values.

    Returns:
        list: list of integers that can be used in ROOTs SetColor methods.
    """
    from ROOT import TColor

    global _colors
    global _color_indices

    if not _colors: # only define colors if they have not been defined already
        rgbcolors = [ # rgb values of MATLABs default palette
            [0, 0.4470, 0.7410],
            [0.8500, 0.3250, 0.0980],
            [0.9290, 0.6940, 0.1250],
            [0.4940, 0.1840, 0.5560],
            [0.4660, 0.6740, 0.1880],
            [0.3010, 0.7450, 0.9330],
            [0.6350, 0.0780, 0.1840]
        ]

        # bind TColor.GetFreeColorIndex to a shorter name for less typing
        # repeatedly calling TColor.GetFreeColorIndex in the list creation is OK
        # since it returns a different index only once it has been used.
        get_idx = TColor.GetFreeColorIndex
        _colors = [(get_idx(), TColor(get_idx(),
                                      rgbcolors[i][0], rgbcolors[i][1], rgbcolors[i][2]))
                   for i in xrange(len(rgbcolors))]

        _color_indices = [col[0] for col in _colors]

    return _color_indices


def set_color(pltable, col):
    """
    Set the passed color for all attributes of the passed plottable.

    Args:
        pltable (ROOT.TObject): Plottable root object
        col (int): Color index to use for the plottable
    """
    col_attributes = ['SetLineColor', 'SetMarkerColor']
    for attr in col_attributes:
        if hasattr(pltable, attr):
            getattr(pltable, attr)(col)


def set_attributes(pltable, **kwargs):
    """
    Set attributes to passed plotable.

    Function checks which functions are available and sets the attributes
    specified by the kwargs.

    Args:
        pltable (ROOT.TObject): Plottable root object

    Keyword Args:
        color (int): TColor index specifying a color to be set
        marker (int): TMarker style index
        line (int): Towline style index
        width (int): TLine width (in pixels)
        size (float): TMarker size
    """
    def conditional_set(key, set_func):
        """
        Helper function that tries to pop key from kwargs and calls the set_func
        with the argument if it is set.

        Args:
            key (str): keyword to check for in kwargs
            set_func (function): Function taking exactly one argument that will
                be called if the key is found in the kwargs.
        """
        arg = kwargs.pop(key, None)
        if arg is not None:
            set_func(pltable, arg)

    arg_func_pairs = (
        ('color', set_color),
        ('marker', lambda p, m: p.SetMarkerStyle(m)),
        ('size', lambda p, s: p.SetMarkerSize(s)),
        ('width', lambda l, w: l.SetLineWidth(w)),
        ('line', lambda l, s: l.SetLineStyle(s))
    )

    for arg, func in arg_func_pairs:
        conditional_set(arg, func)


def plot_on_canvas(can, plots, **kwargs):
    """
    Generic plotting function, that puts all the plots on plots onto the same
    canvas.

    Args:
        can (ROOT.Canvas): Canvas to which the plots should be added
        plots (list): Plottable root objects (specifying Draw() method)

    Keyword Args:
        colors (list. optional): list of colors to be used for plotting
            (otherwise default colors will be used)
        drawOpt (Ste, optional): option that will be passed to the Draw() method
        leg (ROOT.TLegend, optional): Put the passed TLegend onto the plot
        legEntries (list): list of string (at least as long as the plot list)
            from which the keys for the legend are taken

    Returns:
        ROOT.TCanvas: TCanvas that has been passed in
    """
    can.cd()

    colors = kwargs.pop('colors', default_colors())
    get_col = lambda i: colors[ i % len(colors) ]

    leg_option = kwargs.pop('drawOpt', '')
    draw_option = ''.join(['same', leg_option])
    legend = kwargs.pop('leg', None)
    leg_entries = kwargs.pop('legEntries', [h.GetName() for h in plots])

    # cant make pop above default to 'ple', since that alters the
    # draw option defaults
    # If passed draw option is 'H' only, specify the 'ple' option for the
    # legend
    if not leg_option or leg_option == 'H':
        leg_option = 'ple'

    for i in range(len(plots)):
        set_color(plots[i], get_col(i))
        plots[i].Draw(draw_option)

        if legend is not None:
            legend.AddEntry(plots[i], leg_entries[i], leg_option)

    if legend is not None:
        legend.Draw()

    can.Update()

    return can


def mkplot(pltables, **kwargs):
    """
    Plot all pltables onto a canvas and return the canvas

    Args:
        pltables: single plotable ROOT object or list of plotable ROOT objects

    Keyword Args:
        colors (list. optional): list of colors to be used for plotting
            (otherwise default colors will be used)
        drawOpt (Ste, optional): option that will be passed to the Draw() method
        leg (ROOT.TLegend, optional): Put the passed TLegend onto the plot
        legEntries (list): list of string (at least as long as the plot list)
            from which the keys for the legend are taken
        autoRange (bool): Automatically determine the y-range from the passed
            histograms and set it such that all points should fit (not taking)
            into account error bars

    Returns:
        ROOT.TCanvas: Canvas with all the plotables drawn onto it
    """
    # allow single plots to be handled the same as a list of plots
    if not isinstance(pltables, Iterable):
        pltables = [pltables]

    can_name = create_random_str()
    can = r.TCanvas(can_name, '', 50, 50, 600, 600)

    if kwargs.pop('autoRange', False) and all(p.InheritsFrom('TH1')
                                              for p in pltables):
        set_common_range(pltables, axis='y', dscale=0.05)

    plot_on_canvas(can, pltables, **kwargs)

    return can
