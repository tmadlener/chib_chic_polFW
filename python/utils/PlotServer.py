#!/usr/bin/env python
"""
Module that provides a plot server that retrieves plot from a root file
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s: %(message)s')

from decorator import decorator

from utils.misc_helpers import get_full_trigger

def try_file_access(proto_err_msg):
    """
    Decorator for trying to access the root file 'safely'
    """
    @decorator
    def access_wrapper(func, *args):
        """Decorator that will be returned"""
        def try_access(*args):
            """Function that tries the access and returns obj or None"""
            obj_or_dir = func(*args)
            if not obj_or_dir:
                logging.warning('plot server access: ' +
                                proto_err_msg.format(*args))
                return None
            return obj_or_dir
        return try_access(*args)

    return access_wrapper


class PlotServer(object):
    """
    Simple plot server that retrieves and returns TObjects from a TFile with an
    internal TDirectory structure.
    """
    def __init__(self, filename):
        """
        Args:
            filename (str): path to the rootfile from which objects should be
                obtained
        """
        self.hfile = r.TFile.Open(filename)


    @try_file_access('Error in retrieving hist: year={2}, trigger={3}, pt={4}, '
                     'var={5}, frame={6}, plot={7}, dmc={1}')
    def get_hist(self, dmc, year, trigger, pt, var, frame, plot):
        """
        Get a histogram from a file where the internal directory structure uses
        the data/mc, year, trigger and pt information and the plots are
        identified by the variable name, the reference frame and a plot
        identifier.

        Args:
            dmc (str): 'data' or 'mc' identifier
            year (int or str):
            trigger (str): minimum part of a trigger that uniquely identifies it
            pt (int or str): pt bin
            var (str): the variable (e.g. 'costh' or 'phi')
            frame (str): reference frame (in capitals)
            plot (str): which plot (e.g. 'chic1' or 'ratio')

        Returns:
            ROOT.TH1: TH1D found at location specified by arguments or None
        """
        logging.debug('Trying to get plot from file {0}, Parameters are: '
                      'year={1}, trigger={2}, pt={3}, var={4}, frame={5}, '
                      'plot={6}, dmc={7}'.format(self.filename(), year, trigger,
                                                 pt, var, frame, plot, dmc))
        full_trigger = get_full_trigger(trigger)
        subdir = '/'.join([dmc, str(year), full_trigger, str(pt)])
        histname = '_'.join([plot, var, frame])
        return self._get_by_str('/'.join([subdir, histname]))


    @try_file_access('Directory {1}/{2}/{3} not found')
    def get_dir_contents(self, dmc, year, trigger, pt):
        """
        TODO: doc
        """
        logging.debug('Trying to get directory, Parameters are: year={}, '
                      'trigger={}, pt={}, dmc={}'
                      .format(year, trigger, pt, dmc))
        return self._get_by_str('/'.join([dmc, year, trigger, pt]))


    def filename(self):
        """
        Get the filename from which the objects are served

        Returns:
            str: filename
        """
        return self.hfile.GetName()


    def _get_by_str(self, obj_str):
        """
        Get an object from the histo file via a string describing the path to it

        Args:
            obj_str (str): String with path to desired object

        Returns:
            ROOT.TObject found at the place pointed to by obj_str
        """
        logging.debug('Trying to get object {} from file {}'
                      .format(obj_str, self.filename()))
        return self.hfile.Get(obj_str)
