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
        self.legpos = (0.7, 0.8, 0.9, 0.9)
        self.nevent_vars = ['Njpsi', 'Nbkg_jpsi']
        self.components = (
            (self.signal, 7, 632, 'J/#psi'),
            (self.bkg_model, 7, 1, 'BG')
        )

        # self._define_model(wsp)


    def define_model(self, wsp):
        """
        Define the jpsi mass model in the passed workspace

        Args:
            wsp (ROOT.RooWorkspace): workspace into which the fit model is imported
        """
        # rapidity dependent sigma of J/psi
        wsp.factory('CBsigma_p0_jpsi[0.02, 0.015, 0.025]')
        wsp.factory('CBsigma_p1_jpsi[0, -0.01, 0.01]')
        wsp.factory('CBsigma_p2_jpsi[0.01, 0.005, 0.025]')
        CBsigma_jpsi = r.RooFormulaVar('CBsigma_jpsi',
                                       '@0 + @1 * abs(@3) + @2 * abs(@3) * abs(@3)',
                                       r.RooArgList(get_var(wsp, 'CBsigma_p0_jpsi'),
                                                    get_var(wsp, 'CBsigma_p1_jpsi'),
                                                    get_var(wsp, 'CBsigma_p2_jpsi'),
                                                    get_var(wsp, 'JpsiRap')))

        wsp.factory('CBalpha_p0_jpsi[1.729, 1.2, 2.5]')
        wsp.factory('CBalpha_p1_jpsi[0.191, 0, 0.5]')
        CBalpha_jpsi = r.RooFormulaVar('CBalpha_jpsi', '@0 + @1 * abs(@2)',
                                       r.RooArgList(get_var(wsp, 'CBalpha_p0_jpsi'),
                                                    get_var(wsp, 'CBalpha_p1_jpsi'),
                                                    get_var(wsp, 'JpsiRap')))

        wsp.factory('CBmass_p0_jpsi[3.094, 3.086, 3.098]')
        wsp.factory('CBmass_p1_jpsi[0.001, -0.002, 0.002]')
        wsp.factory('CBmass_p2_jpsi[-0.003, -0.005, 0.001]')
        CBmass_jpsi = r.RooFormulaVar('CBmass_jpsi',
                                      '@0 + @1 * abs(@3) + @2 * abs(@3) * abs(@3)',
                                      r.RooArgList(get_var(wsp, 'CBmass_p0_jpsi'),
                                                   get_var(wsp, 'CBmass_p1_jpsi'),
                                                   get_var(wsp, 'CBmass_p2_jpsi'),
                                                   get_var(wsp, 'JpsiRap')))

        # Fraction to ensure that sigma 2 is smaller than sigma 1
        wsp.factory('r_sigma2_jpsi[0.5, 0, 0.9]')
        CBsigma2_jpsi = r.RooFormulaVar('CBsigma2_jpsi', '@0 * @1',
                                        r.RooArgList(get_var(wsp, 'r_sigma2_jpsi'),
                                                     CBsigma_jpsi))

        ws_import(wsp, r.RooArgSet(CBalpha_jpsi, CBmass_jpsi, CBsigma2_jpsi))

        wsp.factory('RooCBShape::{}({}, CBmass_jpsi, CBsigma_jpsi, CBalpha_jpsi,'
                    'CBn_jpsi[2.5, 1.8, 6])'.format('CB_shape1_jpsi', self.mname))
        # wsp.factory('RooCBShape::{}({}, CBmass_jpsi, CBsigma2_jpsi, CBalpha_jpsi,'
        #             ' CBn_jpsi)'.format('CB_shape2_jpsi', self.mname))
        wsp.factory('RooGaussian::{}({}, CBmass_jpsi, CBsigma2_g_jpsi[0.01, 0.0, 0.025])'
                    .format('CB_shape2_jpsi', self.mname))


        wsp.factory('SUM::{}(frac_CB1[0.5, 0, 1] * CB_shape1_jpsi, CB_shape2_jpsi)'
                    .format(self.signal))

        wsp.factory('RooExponential::{}({}, lambda_jpsi[0, -10, 10])'
                    .format(self.bkg_model, self.mname))

        # fix some parameters (as was done previously, values also from there)
        par_fix_vals = [
            # ('CBsigma_p1_jpsi', 0),
            ('CBmass_p1_jpsi', 0),
            ('CBmass_p2_jpsi', 0),
            # ('CBalpha_p1_jpsi', 0),
            # ('CBsigma_p2_jpsi', 0.0125),
            ('CBn_jpsi', 2.5)
        ]

        self.fix_params(wsp, par_fix_vals)

        wsp.factory('{}[1e6, 0, 5e7]'.format(self.nevent_vars[0])) # Njpsi
        wsp.factory('{}[1e5, 0, 5e6]'.format(self.nevent_vars[-1])) # Nbkg_jpsi

        wsp.factory('SUM::{}({} * {},  {} * {})'
                    .format(self.full_model,
                            self.nevent_vars[0], self.signal,
                            self.nevent_vars[-1], self.bkg_model))
