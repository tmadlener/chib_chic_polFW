#!/usr/bin/env python
"""
Script to generate the hepdata yaml tables from the chic2 / chic1 yield ratio
graphs
"""

import numpy as np
import yaml
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True


from utils.graph_utils import get_errors


RATIO_NAME = 'r_chic2_chic1_v_{}_HX_fold_bin_0'

def get_graph(gfile, direction):
    """Get the graph from the passed file"""
    return gfile.Get(RATIO_NAME.format(direction))


def convert_graph(graph, direction):
    """
    Convert the TGraphAsymmErrors to the hepdata expected format

    See:
    https://hepdata-submission.readthedocs.io/en/latest/data_yaml.html
    """
    x_vals  = np.array(graph.GetX())
    x_lo, x_hi, y_lo, y_hi = get_errors(graph)

    # number of digits after comma for bin-borders and bin-centers depends
    # on the angle under consideration
    fmt_bin = lambda x: int(round(x))
    fmt_cent = lambda x: round(x, 2)
    fmt_dep = lambda x: round(x, 4)
    if direction == 'costh':
        fmt_bin = lambda x: round(x, 4)
        fmt_cent = fmt_bin

    # Create the independent variable table in the format that hepdata expectes
    table = {}
    table['independent_variables'] = [
        {'header': {'name': direction, 'units': ''},
         'values': []}
    ]
    for cent, lo, hi in zip(x_vals, x_lo, x_hi):
        table['independent_variables'][0]['values'].append(
            {'low': fmt_bin(cent - lo), 'high': fmt_bin(cent + hi),
             'value': fmt_cent(cent)}
        )

    # Create the dependent variable table
    table['dependent_variables'] = [
        {'header': {'name': 'chi_c2 over chi_c1 ratio'},
         'values': []}
    ]
    y_vals = np.array(graph.GetY())

    for cent, lo, hi in zip(y_vals, y_lo, y_hi):
        table['dependent_variables'][0]['values'].append(
            {'value': fmt_dep(cent),
             'errors': [
                 # hepdata wants the minus sign in front of the lower errors
                 {'asymerror': {'minus': -fmt_dep(lo), 'plus': fmt_dep(hi)},
                  'label': 'stat'}
             ]
            }
        )

    return table


def main(args):
    """Main"""
    rfile = r.TFile.Open(args.graphfile)
    graph = get_graph(rfile, args.direction)

    table = convert_graph(graph, args.direction)
    with open(args.yamlfile, 'w') as yfile:
        yaml.dump(table, yfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to generate a yaml file'
                                     ' in the hepdata expected format.')
    parser.add_argument('graphfile', help='file containing the chic2 / chic1 '
                        'yield ratio graph')
    parser.add_argument('yamlfile', help='output yaml file')

    dir_sel = parser.add_mutually_exclusive_group()
    dir_sel.add_argument('--costh', action='store_const', dest='direction',
                         const='costh',
                         help='Assume that the ratio is vs costh')
    dir_sel.add_argument('--phi', action='store_const', dest='direction',
                         const='phi', help='Assume that the ratio is vs phi')
    parser.set_defaults(direction='costh')


    clargs = parser.parse_args()
    main(clargs)
