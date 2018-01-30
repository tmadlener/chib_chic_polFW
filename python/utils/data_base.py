#!/usr/bin/env python
"""
Module for providing all sorts of common data that is necessary for
all sorts of scripts.

Provides only access to 'database', which is set up as a bunch of json files
in the ./database directory

Attributes:
    json_data_cache (dict): Cached contents of database/data.json
"""

import json
import os
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s: %(message)s')

from decorator import decorator

def default_value_access(def_val, proto_err_msg):
    """
    Decorator for accessing the database 'safely' and getting a warning message
    and a default value in case of failure

    Args:
        def_val: default value that should be returned in case access fails
        proto_err_msg (str): String that calls .format(*args) to log failure

    Returns:
        Element from dictionary or default value
    """
    @decorator
    def access_wrapper(func, *args):
        """The decorator that will be returned"""
        def try_db_access(*args):
            """The function that is actually called and tries the access"""
            try:
                return func(*args)
            except KeyError:
                logging.warning('database access: ' +
                                proto_err_msg.format(*args))
                return def_val

        return try_db_access(*args)

    return access_wrapper


class JsonDataBase(object):
    """
    Simple json file data base access class for getting some commonly used
    information through one place.
    """
    def __init__(self, db_file=''):
        """
        Args:
            db_file (str, optional): path to the database file that should be
                used
        """
        if db_file:
            self.db_path = db_file
        else:
            self.db_path = '/'.join([os.environ['CHIB_CHIC_POLFW_DIR'],
                                     'database', 'data.json'])
        logging.debug('Initializing json data base from {}'
                      .format(self.db_path))
        with open(self.db_path) as dataf:
            self.data = json.load(dataf)


    @default_value_access(-1, 'Trigger {2} not found for year {1}')
    def get_int_lumi(self, year, trigger):
        """
        Get integrated luminosity for a trigger and a given year

        Args:
            year (int or str)
            trigger (str): full name of trigger path

        Returns:
            float: The integrated luminosity in fb^-1 or -1 in case of access
                fail
        """
        logging.debug('Trying to get integrated luminosity for year {} and '
                      'trigger {}'.format(year, trigger))
        return self._get_from_db('int_lumi', lambda d, y, t: d[y][t],
                                 str(year), trigger)


    @default_value_access(-1, 'Year {1} has no entry for energy')
    def get_energy(self, year):
        """"
        Get cms energy for a given year

        Args:
            year (int or str)

        Returns:
            int: cms energy in TeV
        """
        logging.debug('Trying to get cms energy for year {}'.format(year))
        return self._get_from_db('cms_energy', lambda d, y: d[y], str(year))


    def _get_from_db(self, main_key, sub_key_func, *sub_keys):
        """
        Main entry point to get something from the database

        Args:
            main_key (str): the main key in the data base (i.e. top level). This
                sub dictionary will be passed to the sub_key_func
            sub_key_func (func): Function taking a dictionary and a the
                *sub_keys as arguments and returns from the passed dictionary
            *sub_keys: The list of argument that will be passed to the
                sub_key_func

        Returns:
            The return value of the sub_key_func
        """
        main_dict = self.data[main_key]
        return sub_key_func(main_dict, *sub_keys)

