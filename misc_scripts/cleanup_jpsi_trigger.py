#!/usr/bin/env python
"""
Script to "compactify" the trigger information in Jpsi 2012 data tuples
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')


from utils.data_handling import (
    add_branch, get_dataframe, get_treename, check_branch_available
)
from utils.selection_functions import TRIGGER_BIT_MAP


def get_trigger_bit_event(data, trigger):
    """
    Get the trigger bit for the passed trigger
    """
    bit = np.zeros(data.shape[0], dtype=int)
    trig_bit = TRIGGER_BIT_MAP[trigger] # simply error out if we don't find it

    for trig in data.columns:
        bit |= trig_bit * (data.loc[:, trig] == 1)

    return bit


def get_load_vars(infile, trigger):
    """
    Get the list of variables to load
    """
    triggers = [trigger + '*']
    rfile = r.TFile.Open(infile)
    tree = rfile.Get(get_treename(infile))

    # Updating a branch is not possible.
    if check_branch_available(tree, 'trigger', True):
        logging.warn('The branch \'trigger\' is already in the TTree of the '
                     'input file. In place updating of the Branch is not '
                     'possible!')
    return triggers


def main(args):
    """Main"""
    in_data = get_dataframe(args.inputfile, columns=get_load_vars(args.inputfile, args.trigger))
    triggers = get_trigger_bit_event(in_data, args.trigger)
    add_branch(triggers, 'trigger', args.inputfile, get_treename(args.inputfile))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to summarize all '
                                     'trigger information into a bitfield for '
                                     'jpsi data files')
    parser.add_argument('inputfile', help='The file to process')
    parser.add_argument('trigger', help='trigger that should be processed')


    clargs = parser.parse_args()
    main(clargs)
