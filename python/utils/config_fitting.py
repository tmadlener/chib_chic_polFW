#!/usr/bin/env python
"""
Module that can be initialized by a json file containing the fit model and all
settings for fitting
"""

import json

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname) - %(funcName)s: %(message)s')

from utils.FitModel import FitModel
from utils.roofit_utils import fix_params
import utils.RooDoubleCB # Make the double sided CB available via the factory

def _full_model_expr(model_type, name, sub_models, event_yields):
    """
    Construct the full model expression

    NOTE: Might fail at runtime when this string is used as factory expression
    """
    sub_expr = ', '.join([
        ' * '.join(['{}'.format(y), '{}'.format(sub_models[i])])
        for i, y in enumerate(event_yields)
    ])

    return '{}::{}({})'.format(model_type, name, sub_expr)


def _try_factory(wsp, expr):
    """Try to run the expression through the RooWorkspace.factory"""
    obj = wsp.factory(expr)
    if not obj:
        logging.error('\'{}\' failed to create the an object in the workspace'
                      .format(expr))
        return False
    return True


class ConfigFitModel(FitModel):
    """
    Fit model initialized and constructed using a json config file
    """
    def __init__(self, configfile):
        """
        Args:
            configfile (str): Name of the json config file from which the model
                should be constructed
        """
        with open(configfile, 'r') as conf:
            logging.info('Reading config file \'{}\' to define model'
                         .format(configfile))
            config = json.load(conf)

        self.mname = config['fit_variable']
        self.full_model = config['full_model']['name']

        self.nevent_yields = []
        self.components = []
        self.model_expressions = []
        self.expression_strings = []

        sub_models = config["sub_models"]
        for model in sub_models:
            name = model['name']

            # check if the event yield has to be created in the workspace or if
            # it should already be present. If there are no brackets in the
            # expression, then assume that it is already present
            try:
                ev_yield, y_range = model['event_yield'].split('[')
                self.expression_strings.append(model['event_yield'])
            except ValueError:
                ev_yield = model['event_yield']

            self.nevent_yields.append(ev_yield)

            plot_config = model['plot_config']

            self.components.append([name] + plot_config)
            self.model_expressions.append(model['expression']
                                          .format(name, self.mname))

        self.full_model_expr = _full_model_expr(config['full_model']['type'],
                                                self.full_model,
                                                [c[0] for c in self.components],
                                                self.nevent_yields)

        if 'expression_strings' in config:
            self.expression_strings += config['expression_strings']

        self.fix_vars = []
        if 'fix_vars' in config:
            for var in config['fix_vars']:
                self.fix_vars.append((var.keys()[0], var.values()[0]))

        self.floating_costh = []
        if 'floating_costh' in config:
            self.floating_costh = config['floating_costh']

        self.legpos = config['plot_config']['legpos']


    def define_model(self, wsp):
        """
        Build the model from the configuration file into the passed workspace

        Args:
            wsp (ROOT.RooWorkspace): Workspace into which the model should be
                constructed

        Returns:
            success (boolean): True if all of the expressions were succesfully
                imported, False if at least one of them failed
        """
        # first make sure to define all variables in the expression strings,
        # since the models might depend on them. Afterwards build the submodels
        # and finally construct the fullmodel from there
        success = []
        for expr in self.expression_strings:
            success.append(_try_factory(wsp, expr))

        for model in self.model_expressions:
            success.append(_try_factory(wsp, model))
        success.append(_try_factory(wsp, self.full_model_expr))

        fix_params(wsp, self.fix_vars)

        if not all(success):
            logging.error('Could not succesfully define the model in the '
                          'workspace')
            return False

        return True
