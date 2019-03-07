#!/usr/bin/env python
"""
Module containing the basic setup for the chic fitting
"""
import ROOT as r
import ROOT.RooFit as rf

from utils.roofit_utils import ws_import, get_var, fix_params
from utils.constants import (
    m_psiPDG, m_chic0PDG, m_chic1PDG, m_chic2PDG
)
from utils.FitModel import FitModel


class ChicMassModel(FitModel):
    """
    Chic mass model
    """
    def __init__(self, mname):
        """
        Args:
            mname (str): Name of the chic mass variable in the TTree
        """
        self.mname = mname
        self.full_model = 'M_fullModel'
        self.chic1 = 'M_chic1'
        self.chic2 = 'M_chic2'
        self.chic0 = 'M_chic0'
        self.signal = 'M_signal'
        self.bkg_model = 'M_background'
        self.legpos = (0.18, 0.5, 0.38, 0.7)
        self.nevent_yields = ['Nchic0', 'Nchic1', 'Nchic2', 'Nbkg']
        self.floating_costh = []

        self.components = (
            (self.chic0, 7, 901, '#chi_{c0}'),
            (self.chic1, 7, 417, '#chi_{c1}'),
            (self.chic2, 7, 632, '#chi_{c2}'),
            (self.bkg_model, 7, 1, 'BG'),
        )


    def define_model(self, ws):
        """
        Build the chic mass model in the passed workspace

        Args:
            ws (ROOT.RooWorkspace): Workspace into which the mass model is
                constructed
        """
        ## build mass model
        ws.factory('RooCBShape::{}({}, CBmass1[3.5,3.45,3.54],'
                   'CBsigma1[0.008, 0.003, 0.02], CBalpha1[0.6, 0.2, 1.1], '
                   'CBn[2.5, 1.8, 3.2])'.format(self.chic1, self.mname))

        pes = r.RooFormulaVar('PES', '(@0 - {0}) / {1}'.format(m_psiPDG, m_chic1PDG - m_psiPDG),
                              r.RooArgList(get_var(ws, 'CBmass1')))
        ws_import(ws, pes)
        CBmass0 = r.RooFormulaVar('CBmass0', '(@0 * {0}) + {1}'.format(m_chic0PDG - m_psiPDG, m_psiPDG),
                                  r.RooArgList(get_var(ws, 'PES')))
        CBmass2 = r.RooFormulaVar('CBmass2', '(@0 * {0}) + {1}'.format(m_chic2PDG - m_psiPDG, m_psiPDG),
                                  r.RooArgList(get_var(ws, 'PES')))
        ws_import(ws, r.RooArgSet(CBmass0, CBmass2))

        # ws.factory('CBmass2[3.54, 3.49, 3.58]')
        # ws.factory('CBmass0[3.42, 3.395, 3.45]')

        CBsigma0 = r.RooFormulaVar('CBsigma0', '@0 * {0}'.format((m_chic0PDG - m_psiPDG) / (m_chic1PDG - m_psiPDG)),
                                   r.RooArgList(get_var(ws, 'CBsigma1')))
        CBsigma2 = r.RooFormulaVar('CBsigma2', '@0 * {0}'.format((m_chic2PDG - m_psiPDG) / (m_chic1PDG - m_psiPDG)),
                                   r.RooArgList(get_var(ws, 'CBsigma1')))
        ws_import(ws, r.RooArgSet(CBsigma0, CBsigma2))

        # ws.factory('CBsigma2[0.008, 0.003, 0.02]')
        # ws.factory('CBsigma0[0.008, 0.003, 0.02]')

        ws.factory('RooCBShape::{}({}, CBmass2, CBsigma2, CBalpha2[0.6, 0.2, 1.1], CBn)'
                   .format(self.chic2, self.mname))
        # ws.factory('RooCBShape::M_chic2(chicMass, CBmass2, CBsigma2, CBalpha1, CBn)')
        CBalpha0 = r.RooFormulaVar('CBalpha0', '(@0 + @1) / 2.0',
                                   r.RooArgList(get_var(ws, 'CBalpha1'),
                                                get_var(ws, 'CBalpha2')))
        # r.RooArgList(ws.var('CBalpha1'), ws.var('CBalpha1')))
        ws_import(ws, r.RooArgSet(CBalpha0))

        ws.factory('RooVoigtian::{}({}, CBmass0, CBsigma0, CBwitdh[0.0104])'
                   .format(self.chic0, self.mname))

        ## background
        # ws.factory('q01S[3.1, 3.0, 3.2]')
        # ws.factory('alpha1[0.6, 0, 5.0]')
        # ws.factory('beta1[-2.5, -10, 0]')

        # a1 = r.RooFormulaVar('a1', 'TMath::Abs(@0 - @1)',
        #                      r.RooArgList(get_var(ws, 'chicMass'),
        #                                   get_var(ws, 'q01S')))
        # b1 = r.RooFormulaVar('b1', '@0 * (@1 - @2)',
        #                      r.RooArgList(get_var(ws, 'beta1'),
        #                                   get_var(ws, 'chicMass'),
        #                                   get_var(ws, 'q01S')))
        # signum1 = r.RooFormulaVar('signum1',
        #                           '(TMath::Sign(-1., @0 - @1) + 1) / 2.',
        #                           r.RooArgList(get_var(ws, self.mname),
        #                                        get_var(ws, 'q01S')))
        # ws_import(ws, r.RooArgSet(a1, b1, signum1))

        # M_background = r.RooGenericPdf(self.bkg_model, 'signum1 * pow (a1, alpha1) * exp(b1)',
        #                                r.RooArgList(ws.function('a1'), ws.function('signum1'), ws.function('b1'),
        #                                             ws.var('alpha1')))

        # polynomial background
        ws.factory('BK_p1[0, -1, 1]')
        ws.factory('BK_p2[0, -1, 1]')
        M_background = r.RooPolynomial(self.bkg_model, self.bkg_model, ws.var('chicMass'),
                                      r.RooArgList(ws.var('BK_p1'), ws.var('BK_p2')))

        # M_background = r.RooExponential(self.bkg_model, self.bkg_model, ws.var('chicMass'), ws.var('BK_p1'))


        ws_import(ws, M_background)

        fix_params(ws, [
            ('CBn', 2.75),
            ('BK_p2', 1e-10),
            # ('q01S', 3.1)
        ])
        ws.var('CBmass1').setVal(3.510)
        # ws.var('alpha1').setConstant(True)
        # ws.var('q01S').setConstant(True)
        # ws.var('beta1').setConstant(True)

        ## combine model
        ws.factory('{}[10000, 0, 200000]'.format(self.nevent_vars[1])) # Nchic1
        ws.factory('{}[300, 0, 10000]'.format(self.nevent_vars[0])) # Nchic0
        ws.factory('{}[30000, 0, 200000]'.format(self.nevent_vars[-1])) # Nbkg

        ## leave Nchic2 free and make ratio a formula var
        ws.factory('{}[10000, 0, 200000]'.format(self.nevent_vars[2])) # Nchic2
        r_c2_c1 = r.RooFormulaVar('r_chic2_chic1', '@0 / @1',
                                  r.RooArgList(get_var(ws, self.nevent_vars[2]),
                                               get_var(ws, self.nevent_vars[1])))
        ws_import(ws, r_c2_c1)

        ## make ratio free and Nchic2 depend on it
        # ws.factory('r_chic2_chic1[0.5, 0, 1]')
        # Nchic2 = r.RooFormulaVar(self.nevent_vars[2], '@0 * @1',
        #                          r.RooArgList(get_var(ws, 'r_chic2_chic1'),
        #                                       get_var(ws, self.nevent_vars[1])))
        # ws_import(ws, Nchic2)


        Nsignal = r.RooFormulaVar('Nsignal', '@0 + @1 + @2',
                                  r.RooArgList(get_var(ws, self.nevent_vars[1]),
                                               get_var(ws, self.nevent_vars[2]),
                                               get_var(ws, self.nevent_vars[0])))
        ws_import(ws, Nsignal)

        ws.factory('SUM::{}({} * {}, {} * {}, {} * {})'
                   .format(self.signal,
                           self.nevent_vars[1], self.chic1,
                           self.nevent_vars[2], self.chic2,
                           self.nevent_vars[0], self.chic0))
        ws.factory('SUM::{}(Nsignal * {}, {} *  {})'
                   .format(self.full_model, self.signal,
                           self.nevent_vars[-1], self.bkg_model))
