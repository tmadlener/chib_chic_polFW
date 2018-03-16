#!/usr/bin/env python

from utils.FitModel import FitModel

"""
Module containing the basic setup for the chib fitting
"""

import json

class ChibMassModel(FitModel):
    """
    Chib mass model
    """
    def __init__(self, config_file):
        """
        Args:
            config_file: Name of the json file containing the model definition
        """
        self.config_file = config_file

        with open(config_file, 'r') as f:
            data = json.load(f)
            self.mname = data["chi_fitvar"]["name"]
            self.full_model = data["chi_modelname"]
            self.model_strings = data["chi_model"]

        self.components = (
            ('chib1', 7, 417, '#chi_{b1}'),
            ('chib2', 7, 632, '#chi_{b2}'),
            ('background', 7, 1, 'background')
            )
        self.legpos = (0.18, 0.5, 0.38, 0.7)


    def define_model(self, ws):
        """
        Build the chib mass model in the passed workspace

        Args:
            ws (ROOT.RooWorkspace): Workspace into which the mass model is
                constructed
        """  

        for fac in self.model_strings:
            ws.factory(fac)

