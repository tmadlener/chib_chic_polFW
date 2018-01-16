"""
Module containing helper functions for plotting and related stuff

Attributes:
    _colors (list): list of index and TColor for caching default colors
    _color_indices (list): list of indices for the different default colors,
        which will be used in ROOTs SetColor methods

    Neither of the two above variables should be accessed directly, but the
    default_colors method should be used
"""

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
    pltable.SetLineColor(col)
    pltable.SetMarkerColor(col)


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


