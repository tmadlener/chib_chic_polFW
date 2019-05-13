#!/usr/bin/env python
"""
Module to define functionality to load root libs
"""

import os

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.ERROR,
                    format='%(levelname)s - %(funcName)s: %(message)s')

def try_load_lib(path, name=''):
    """
    Try and load the .so specified by path using the ROOT module loading
    functionality.

    Args:
        path (str): Path to the .so file to load

    Raises:
        ImportError: In case ROOT is not able to load the library
    """
    if r.gSystem.Load(path) != 0:
        # For some reason only loading the shared object below fails with an
        # undefined symbol error sometimes. This can normally be fixed by
        # loading libRooFitCore first. In principle this could be done
        # independently of first trying to load the shared object, but loading
        # libRooFitCore displays the splash text from RooFit which is something
        # I do not like entirely, so we do it only when really necessary
        logging.debug('Cannot load the {} shared object from {}. Trying after '
                      'loading libRooFitCore first'.format(name, path))
        r.gSystem.Load('libRooFitCore')

        if r.gSystem.Load(path) != 0:
            # If we still fail to load the library then bail out. This is a
            # bigger problem than we thought
            raise ImportError('Cannot load the {} shared object from {}.\nMake'
                              ' sure that it is compiled and linked against the'
                              ' proper libraries.\nCheck {} for this and run '
                              '\'make clean; make\' there and see if this fixes'
                              ' things.'.format(name, path, os.path.dirname(path)))
