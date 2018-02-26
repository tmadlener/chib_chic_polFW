#!/usr/bin/env python
"""
Module containing the basic setup for the jpsi mass fitting
"""
import ROOT as r
import ROOT.RooFit as rf

from utils.roofit_utils import ws_import, get_var
from utils.FitModel import FitModel

class JpsiMassModel(FitModel):
    """
    Jpsi mass model
    """
    def __init__(self, mname):
        """
        Args:
            mname (str): Name of the jpsi mass variable in the TTree
        """
        self.signal = 'M_jpsi_signal'
        self.bkg_model = 'M_bkg_jpsi'
        self.full_model = 'M_fullModel_jpsi'
        self.mname = mname

        self.components = (
            (self.signal, 7, 417),
            (self.bkg_model, 7, 1)
        )

        # self._define_model(wsp)


    def define_model(self, wsp):
        """
        Define the jpsi mass model in the passed workspace

        Args:
            wsp (ROOT.RooWorkspace): workspace into which the fit model is imported
        """
        wsp.factory('RooCBShape::{}({}, CBmass_jpsi[3.1, 2.95, 3.24],'
                    'CBsigma_jpsi[0.0035, 0, 0.5], CBalpha_jpsi[0.6, 0.2, 2.5], '
                    'CBn_jpsi[2.5, 1.8, 6])'.format(self.signal, self.mname))

        get_var(wsp, 'CBalpha_jpsi').setConstant(True)
        get_var(wsp, 'CBn_jpsi').setConstant(True)

        wsp.factory('RooExponential::{}({}, lambda_jpsi[0, -10, 10])'
                    .format(self.bkg_model, self.mname))

        wsp.factory('Njpsi[1e6, 0, 5e7]')
        wsp.factory('Nbkg_jpsi[1e5, 0, 5e6]')

        wsp.factory('SUM::{}(Njpsi * {},  Nbkg_jpsi * {})'
                    .format(self.full_model, self.signal, self.bkg_model))
