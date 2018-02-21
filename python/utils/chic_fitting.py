#!/usr/bin/env python
"""
Module containing the basic setup for the chic fitting
"""

import ROOT as r
import ROOT.RooFit as rf

from utils.roofit_utils import ws_import

def chi_mass_model(ws, mname='chicMass'):
    """
    Build the chic mass model in the passed workspace

    Args:
        ws (ROOT.RooWorkspace): Workspace into which the mass model is constructed
        mname (str, optional): Name of the chic mass variable in the passed
            workspace
    """
    # Constant values used in the fit
    m_psiPDG = 3.096916
    m_chic0PDG = 3.41475
    m_chic1PDG = 3.51066
    m_chic2PDG = 3.55620

    ## build mass model
    ws.factory('RooCBShape::M_chic1({}, CBmass1[3.5,3.45,3.54],'
               'CBsigma1[0.008, 0.003, 0.02], CBalpha1[0.6, 0.2, 1.1], '
               'CBn[2.5, 1.8, 3.2])'.format(mname))

    pes = r.RooFormulaVar('PES', '(@0 - {0}) / {1}'.format(m_psiPDG, m_chic1PDG - m_psiPDG),
                          r.RooArgList(ws.var('CBmass1')))
    ws_import(ws, pes)
    CBmass0 = r.RooFormulaVar('CBmass0', '(@0 * {0}) + {1}'.format(m_chic0PDG - m_psiPDG, m_psiPDG),
                              r.RooArgList(ws.function('PES')))
    CBmass2 = r.RooFormulaVar('CBmass2', '(@0 * {0}) + {1}'.format(m_chic2PDG - m_psiPDG, m_psiPDG),
                              r.RooArgList(ws.function('PES')))
    ws_import(ws, r.RooArgSet(CBmass0, CBmass2))

    # ws.factory('CBmass2[3.54, 3.49, 3.58]')
    # ws.factory('CBmass0[3.42, 3.395, 3.45]')

    CBsigma0 = r.RooFormulaVar('CBsigma0', '@0 * {0}'.format((m_chic0PDG - m_psiPDG) / (m_chic1PDG - m_psiPDG)),
                               r.RooArgList(ws.var('CBsigma1')))
    CBsigma2 = r.RooFormulaVar('CBsigma2', '@0 * {0}'.format((m_chic2PDG - m_psiPDG) / (m_chic1PDG - m_psiPDG)),
                               r.RooArgList(ws.var('CBsigma1')))
    ws_import(ws, r.RooArgSet(CBsigma0, CBsigma2))

    # ws.factory('CBsigma2[0.008, 0.003, 0.02]')
    # ws.factory('CBsigma0[0.008, 0.003, 0.02]')

    ws.factory('RooCBShape::M_chic2(chicMass, CBmass2, CBsigma2, CBalpha2[0.6, 0.2, 1.1], CBn)')
    # ws.factory('RooCBShape::M_chic2(chicMass, CBmass2, CBsigma2, CBalpha1, CBn)')
    CBalpha0 = r.RooFormulaVar('CBalpha0', '(@0 + @1) / 2.0',
                               r.RooArgList(ws.var('CBalpha1'), ws.var('CBalpha2')))
    # r.RooArgList(ws.var('CBalpha1'), ws.var('CBalpha1')))
    ws_import(ws, r.RooArgSet(CBalpha0))

    ws.factory('RooVoigtian::M_chic0(chicMass, CBmass0, CBsigma0, CBwitdh[0.0104])')

    ## background
    ws.factory('q01S[3.1, 3.0, 3.2]')
    ws.factory('alpha1[0.6, 0, 5.0]')
    ws.factory('beta1[-2.5, -10, 0]')

    a1 = r.RooFormulaVar('a1', 'TMath::Abs(@0 - @1)', r.RooArgList(ws.var('chicMass'), ws.var('q01S')))
    b1 = r.RooFormulaVar('b1', '@0 * (@1 - @2)', r.RooArgList(ws.var('beta1'), ws.var('chicMass'), ws.var('q01S')))
    signum1 = r.RooFormulaVar('signum1', '(TMath::Sign(-1., @0 - @1) + 1) / 2.',
                              r.RooArgList(ws.var('chicMass'), ws.var('q01S')))
    ws_import(ws, r.RooArgSet(a1, b1, signum1))

    # ws.factory('BK_p1[0, -1, 1]')
    # ws.factory('BK_p2[0, -1, 1]')
    # M_background = r.RooPolynomial('M_background', 'M_background', ws.var('chicMass'),
    #                               r.RooArgList(ws.var('BK_p1'), ws.var('BK_p2')))

    M_background = r.RooGenericPdf('M_background', 'signum1 * pow (a1, alpha1) * exp(b1)',
                                   r.RooArgList(ws.function('a1'), ws.function('signum1'), ws.function('b1'),
                                                ws.var('alpha1')))

    ws_import(ws, M_background)

    # ws.var('BK_p2').setVal(1e-10)
    # ws.var('BK_p2').setConstant(True)
    ws.var('CBn').setVal(2.75)
    ws.var('CBn').setConstant(True)
    ws.var('CBmass1').setVal(3.510)
    # ws.var('alpha1').setConstant(True)
    # ws.var('q01S').setConstant(True)
    # ws.var('beta1').setConstant(True)

    ## combine model
    ws.factory('Nchic1[10000, 0, 200000]')
    ws.factory('r_chic2_chic1[0.5, 0, 1]')
    ws.factory('Nchic0[300, 0, 10000]')
    ws.factory('Nbkg[30000, 0, 200000]')

    Nchic2 = r.RooFormulaVar('Nchic2', '@0 * @1', r.RooArgList(ws.var('r_chic2_chic1'), ws.var('Nchic1')))
    ws_import(ws, Nchic2)
    Nsignal = r.RooFormulaVar('Nsignal', '@0 + @1 + @2',
                              r.RooArgList(ws.var('Nchic1'), ws.function('Nchic2'), ws.var('Nchic0')))
    ws_import(ws, Nsignal)

    ws.factory('SUM::M_signal(Nchic1 * M_chic1, Nchic2 * M_chic2, Nchic0 * M_chic0)')
    ws.factory('SUM::M_fullModel(Nsignal * M_signal, Nbkg *  M_background)')


def make_mass_fit_plot(ws, pdfname, mname='chicMass', add_cut=''):
    """
    Make plot of mass fit

    Args:
        ws (ROOT.RooWorkspace): workspace containing the data as well as the fit
           pdf and the fit results
        pdfname (str): Name of the pdffile that will be created and stored
        mname (str, optional): Name of the chic mass variable in the workspace
        add_cut (str, optional): Any additional cuts that should be applied to
            the full dataset contained in the workspace before plotting
    """
    plot_data = ws.data('full_data').reduce(add_cut)
    frame = ws.var(mname).frame(rf.Bins(80))
    full_pdf = ws.pdf('M_fullModel')

    plot_data.plotOn(frame, rf.MarkerSize(0.8))
    full_pdf.plotOn(frame, rf.LineWidth(2))
    full_pdf.plotOn(frame, rf.Components('M_chic0'), rf.LineStyle(7),
                   rf.LineColor(901), rf.LineWidth(2))
    full_pdf.plotOn(frame, rf.Components('M_chic1'), rf.LineStyle(7),
                   rf.LineColor(417), rf.LineWidth(2))
    full_pdf.plotOn(frame, rf.Components('M_chic2'), rf.LineStyle(7),
                   rf.LineColor(632), rf.LineWidth(2))
    full_pdf.plotOn(frame, rf.Components('M_background'), rf.LineStyle(7),
                   rf.LineColor(1), rf.LineWidth(2))

    can = r.TCanvas('can', 'can', 50, 50, 600, 600)
    can.cd()
    frame.Draw()

    can.SaveAs(pdfname)


def do_mass_fit(ws, savename, add_cut=''):
    """
    Run the mass fit (extended maximum likelihood)

    Args:
        ws (ROOT.RooWorkspace): workspace containing the model to fit and the
            data to fit to
        savename (str): name under which the snapshot of the variables after the
            fit as well as the fit results will be stored into the workspace
            ('snap_{savename}' and 'fit_res_{savename}')
        add_cut (str, optional): Any additional cuts that should be applied to
            the full dataset contained in the workspace before doing the fit
    """
    fit_data = ws.data('full_data')
    if add_cut:
        fit_data = fit_data.reduce(add_cut) #NOTE: possibly redundant

    fit_results = ws.pdf('M_fullModel').fitTo(fit_data, rf.Minos(True),
                                              rf.NumCPU(4), rf.Save(True),
                                              rf.Extended(True), rf.Offset(False))

    fit_results.Print()
    print fit_results.status(), fit_results.covQual()

    ws.saveSnapshot('snap_{}'.format(savename), ws.allVars())
    fit_results.SetName('fit_res_{}'.format(savename))
    ws_import(ws, fit_results)
