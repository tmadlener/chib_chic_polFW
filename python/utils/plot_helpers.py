"""
Module containing helper functions for plotting and related stuff

Attributes:
    _colors (list): list of index and TColor for caching default colors
    _color_indices (list): list of indices for the different default colors,
        which will be used in ROOTs SetColor methods

    Neither of the two above variables should be accessed directly, but the
    default_colors method should be used
"""

import numpy as np

import ROOT as r

from itertools import product
from utils.misc_helpers import create_random_str, make_iterable
from utils.hist_utils import set_common_range, set_labels, divide

_colors = []
_color_indices = []

class TCanvasWrapper(object):
    """
    Class that wraps a TCanvas and exposes all its functionality externally
    but also puts attached plots and other TObjects to the members so that they
    stay alive as long as the TCanvas itself
    """
    def __init__(self, canvas):
        """
        Initialize from an already existing TCanvas
        """
        self.wrapped_canvas = canvas
        self.pltables = []
        self.attached_tobjects = []

    def __getattr__(self, attr):
        """
        Forward all calls to the canvas
        """
        # intercept the adding function calls
        if attr.startswith('add'):
            self.__getattribute__(attr)

        # remove all attached things when the canvas itself is cleared
        if attr == 'Clear' or attr == 'Reset':
            self.pltables = []
            self.attached_tobjects = []

        return self.wrapped_canvas.__getattribute__(attr)

    def add_pltables(self, pltables):
        """
        Add plotable objects
        """
        for plt in make_iterable(pltables):
            # do not store duplicates
            if not plt in self.pltables:
                self.pltables.append(plt)

    def add_tobject(self, tobject):
        """
        Add any TObject
        """
        # do not store duplicates
        if not tobject in self.attached_tobjects:
            self.attached_tobjects.append(tobject)


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
            [0, 0.4470, 0.7410],      # darkish blue
            [0.8500, 0.3250, 0.0980], # "orangeish" red
            [0.9290, 0.6940, 0.1250], # yellow
            [0.4940, 0.1840, 0.5560], # purple
            [0.4660, 0.6740, 0.1880], # green
            [0.3010, 0.7450, 0.9330], # lightish blue
            [0.6350, 0.0780, 0.1840]  # wine red
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


def default_attributes(diff_markers=True, size=2, linewidth=2,
                       open_markers=False, unique_colors=True):
    """
    Get a list of some sensible default attributes

    Args:
        diff_markers (bool, optional): Use different markers (True)
        size (float, optional): Size of the markers (2)
        linewidth (int, optional): Linewidth of TLines (2)
        open_markers (bool, optional): Use open markers instead of filled ones
            (False)
        unique_colors (bool): Match each marker with an individual color instead
            of using all available combinations of color and marker (True)

    Returns:
        list: List containing a list of dictionaries as taken by the 'attr'
            argument of plot_on_canvas
    """
    dcols = default_colors()
    if diff_markers:
        if open_markers:
            markers = [24, 26, 32, 25, 27]
        else:
            markers = [20, 22, 23, 21, 33]
    else:
        markers = [6]

    # first combine different colors with different markers
    if unique_colors: # each marker gets its own color
        mark_col = zip(markers, dcols)
    else: # make all possible combinations, iterating colors first
        mark_col = product(markers, dcols)

    attributes = []
    for marker, color in mark_col:
        attributes.append({
            'color': color, 'marker': marker, 'size': size, 'width': linewidth
        })

    return attributes


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
        fill (int): FillColor
        fillalpha (tuple): TColor index specifying fill color and float between
            0 and 1 specifying fill alpha
        fillstyle (int): A fillstyle to be set
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
        ('line', lambda l, s: l.SetLineStyle(s)),
        ('fill', lambda p, c: p.SetFillColor(c)),
        ('fillalpha', lambda p, fa: p.SetFillColorAlpha(fa[0], fa[1])),
        ('fillstyle', lambda p, s: p.SetFillStyle(s))
    )

    for arg, func in arg_func_pairs:
        conditional_set(arg, func)


def setup_legend(*position):
    """
    Setup and return a ROOT.TLegend object

    Args:
        position: Coordinates of the TLegend. Forwarded to the constructor

    Returns:
        ROOT.TLegend
    """
    leg = r.TLegend(*position)
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetEntrySeparation(0.01)
    leg.SetBorderSize(0)

    return leg


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
        drawOpt (str, optional): option that will be passed to the Draw() method
        leg (ROOT.TLegend, optional): Put the passed TLegend onto the plot
        legEntries (list): list of string (at least as long as the plot list)
            from which the keys for the legend are taken
        attr (list of dicts, optional): list of attributes to pass to the
            set_attributes function. Overrides color argument!

    Returns:
        ROOT.TCanvas: TCanvas that has been passed in
    """
    can.cd()

    attributes = kwargs.pop('attr', None)
    if attributes is None:
        colors = kwargs.pop('colors', default_colors())
        attributes = []
        for col in colors:
            attributes.append({'color': col})

    get_att = lambda i: attributes[ i % len(attributes) ]

    leg_option = kwargs.pop('legOpt', '')
    draw_option = kwargs.pop('drawOpt', '')
    if not leg_option:
        leg_option = draw_option
    # draw_option = ''.join(['same', draw_option])

    legend = kwargs.pop('leg', None)
    leg_entries = kwargs.pop('legEntries', [h.GetName() for h in plots])

    # cant make pop above default to 'ple', since that alters the
    # draw option defaults
    # If passed draw option is 'H' only, specify the 'ple' option for the
    # legend
    if not leg_option or leg_option == 'H':
        leg_option = 'ple'

    for i, plot in enumerate(plots):
        set_attributes(plot, **get_att(i))
        if i == 0:
            plot.Draw(draw_option)
        else:
            plot.Draw(draw_option + 'same')

        if legend is not None:
            legend.AddEntry(plot, leg_entries[i], leg_option)

    if legend is not None:
        legend.Draw()

    can.Update()

    return can


def mkplot(pltables, **kwargs):
    """
    Plot all pltables onto a canvas and return the canvas

    Args:
        pltables: single plotable ROOT object or list of plotable ROOT objects.
            If this is an empty list an empty canvas (or the passed in canvas)
            will be returned

    Keyword Args:
        colors (list. optional): list of colors to be used for plotting
            (otherwise default colors will be used)
        drawOpt (str, optional): option that will be passed to the Draw() method
        leg (ROOT.TLegend, optional): Put the passed TLegend onto the plot
        legEntries (list): list of string (at least as long as the plot list)
            from which the keys for the legend are taken
        yRange, xRange (list): Two element list defining the range for the given
            axis. Valid options are: two floats or any mixture of one float and
            a None value or even two None values. For any passed value that
            value will be used as axis range, for any None value an appropriate
            value will be determined from the passed pltables
        can (ROOT.TCanvas): Do not create new canvas but use passed canvas to
            plot on
        log[xyz] (boolean, optional): Set the [xyz] axis to log scale

    Returns:
        TCanvasWrapper: Transparent wrapper class around a TCanvas that forwards
            all calls to the TCanvas but keeps the plotted objects alive as long
            as the TCanvas is alive

    See also:
        plot_on_canvas
    """
    # allow single plots to be handled the same as a list of plots
    pltables = make_iterable(pltables)

    can = kwargs.pop('can', None)
    if can is None:
        can_name = create_random_str()
        can = r.TCanvas(can_name, '', 50, 50, 600, 600)

    if not isinstance(can, TCanvasWrapper):
        can = TCanvasWrapper(can)

    # Check if at least one pltable has been passed and return the canvas if not
    # NOTE: Can't return earlier with None, since a canvas is expected to be
    # returned, and only here is it certain, that we return the right canvas
    if len(pltables) < 1:
        return can

    only_hists = all(p.InheritsFrom('TH1') for p in pltables)

    x_range = kwargs.pop('xRange', None)
    if x_range is not None and only_hists:
        set_common_range(pltables, axis='x', dscale=0.05, drange=x_range)
    y_range = kwargs.pop('yRange', None)
    if y_range is not None and only_hists:
        set_common_range(pltables, axis='y', dscale=0.05, drange=y_range)

    y_label = kwargs.pop('yLabel', '')
    x_label = kwargs.pop('xLabel', '')
    if (y_label or x_label) and only_hists:
        for h in pltables:
            set_labels(h, x_label, y_label)

    leg = kwargs.pop('leg', None)
    if leg is None:
        legPos = kwargs.pop('legPos', None)
        if legPos is not None and len(legPos) == 4:
            leg = setup_legend(*legPos)

    plot_on_canvas(can, pltables, leg=leg, **kwargs)

    can.add_pltables(pltables)
    if leg is not None:
        can.add_tobject(leg)

    if kwargs.pop('logy', False):
        can.SetLogy()
    if kwargs.pop('logx', False):
        can.SetLogx()
    if kwargs.pop('logz', False):
        can.SetLogz()

    return can


def setup_latex():
    """
    Setup a TLatex for putting text onto the plot

    Returns:
        ROOT.TLatex: TLatex object with basic settings
    """
    latex = r.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)

    return latex


def put_on_latex(latex, text_info):
    """
    Put all the text onto the passed latex

    Args:
        latex (ROOT.TLatex): The TLatex that will be used for drawing
        text_info (list of tuples): List containing the position and the
            text to put onto the plot in the format (leftpos, toppos, text)
    """
    for left, top, text in text_info:
        latex.DrawLatex(left, top, text)


def _set_ratio_properties(hist):
    """
    Set the properties of the histogram so that they look somewhat consistent
    for the ratio pad of the baseline plot
    """
    hist.GetYaxis().SetTitleSize(0.08)
    hist.GetXaxis().SetTitleSize(0.08)
    hist.GetYaxis().SetLabelSize(0.08)
    hist.GetXaxis().SetLabelSize(0.08)
    hist.GetYaxis().SetTitleOffset(0.5)


def baseline_plot(baseline, compplots, **kwargs):
    """
    Make a plot and compare the compplots with the baseline plots
    """
    comp_attr = kwargs.pop('attr', None)
    if comp_attr is None:
        comp_attr = default_attributes(open_markers=False, size=1.0)

    # the baseline will always be black. Try to match the size of the markers to
    # the one that were requested by the user
    base_attr = {'color': 1, 'marker': 20, 'size': 1.5}
    sizes = np.unique([a['size'] for a in comp_attr if 'size' in a])
    if len(sizes) == 1:
        base_attr['size'] = sizes[0]

    attr = [base_attr] + comp_attr

    # use the xLabel only for the lower plot
    xLabel = kwargs.pop('xLabel', None)

    # add the baseline name to the legend entries (if any)
    legEntries = kwargs.pop('legEntries', None)
    base_name = kwargs.pop('basename', 'baseline')
    if legEntries is not None:
        legEntries = [base_name] + legEntries

    # setup canvas
    can = kwargs.pop('can', None)
    if can is None:
        can_name = create_random_str()
        can = r.TCanvas(can_name, '', 50, 50, 600, 600)
    can.cd()

    ppad = r.TPad('_'.join([can.GetName(), 'plotpad']), 'plotpad', 0, 0.3, 1, 1)
    r.SetOwnership(ppad, False)
    ppad.Draw()

    # plot the comparison plot
    ppad = mkplot([baseline] + make_iterable(compplots),
                  attr=attr, can=ppad,
                  legEntries=legEntries,
                  **kwargs)

    can.cd()
    rpad = r.TPad('_'.join([can.GetName(), 'ratiopad']), 'rpad', 0, 0, 1, 0.33)
    rpad.SetBottomMargin(0.2)
    r.SetOwnership(rpad, False)
    rpad.Draw()

    # remove some kwargs
    for kwarg in ['yLabel', 'legPos', 'leg', 'legEntries', 'yRange', 'logy']:
        kwargs.pop(kwarg, None)

    ratios = [divide(p, baseline) for p in make_iterable(compplots)]
    for ratio in ratios:
        _set_ratio_properties(ratio)

    # determine the ratios and plot them
    rpad = mkplot(ratios,
                  attr=comp_attr, can=rpad, xLabel=xLabel,
                  yLabel=' / '.join([kwargs.pop('compname', 'distribution(s)'),
                                     base_name]),
                  yRange=kwargs.pop('yRangeRatio', [None, None]), **kwargs)

    # attach all plots to the returned canvas
    if not isinstance(can, TCanvasWrapper):
        can = TCanvasWrapper(can)

    # can.add_pltables(ppad.pltables + rpad.pltables)
    # for obj in ppad.attached_tobjects:
    #     can.add_tobject(obj)
    can.add_tobject(ppad)
    can.add_tobject(rpad)

    return can
