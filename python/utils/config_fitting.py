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
            self.expression_strings.append(model['event_yield'])
            # strip the range from the variable name
            self.nevent_yields.append(model['event_yield'].split('[')[0])

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
        """
        # first make sure to define all variables in the expression strings,
        # since the models might depend on them. Afterwards build the submodels
        # and finally construct the fullmodel from there
        for expr in self.expression_strings:
            wsp.factory(expr)
        for model in self.model_expressions:
            wsp.factory(model)
        wsp.factory(self.full_model_expr)

        self.fix_params(wsp, self.fix_vars)
        wsp.Print()
