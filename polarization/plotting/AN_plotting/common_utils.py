#!/usr/bin/env python
"""
Module containing commonly used functionality of the scripts in this directory
"""

from utils.graph_utils import scale_graph, fit_to_graph

def get_graph(rfile, param, var='costh'):
    """Get the ratio graph from the file"""
    return rfile.Get('{}_v_{}_HX_fold_bin_0'.format(param, var))


def rescale_graph(graph, temp_graph):
    """
    Rescale the graph to the template graph (and also shift values along the x-axis)
    """
    norm, _ = fit_to_graph(graph, temp_graph)

    return scale_graph(graph, 1. / norm)

