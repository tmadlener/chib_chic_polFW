#!/usr/bin/env python
"""
Module to facilitate some handling of RooFit objects
"""
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import pandas as pd

from utils.misc_helpers import get_np_from_tmatrix


def ws_import(wsp, *args):
    """
    Import to workspace.

    Since the 'import' keyword is reserved in python, it is not possible to
    call workspace.import directly. This is just a wrapper to save some typing

    Args:
        wsp (ROOT.RooWorkspace): workspace into which everything should be
            imported
        args: Objects that can be passed to ROOT.RooWorkspace.import
    """
    getattr(wsp, 'import')(*args)


def get_var(wsp, varname):
    """
    Get the variable or function from the workspace

    Args:
        wsp (ROOT.RooWorkspace): workspace containing the desired variable or
            function
        varname (str): name of the variable or function

    Returns:
        ROOT.RooRealVar or ROOT.RooFormulaVar: Variable or function found under
            the passed name in the workspace
    """
    var = wsp.var(varname)
    if var:
        return var
    var = wsp.function(varname)
    if var:
        return var


def get_var_err(wsp, varname, fit_res=None):
    """
    Get the value and error of a variable from the workspace

    Args:
        wsp (ROOT.RooWorkspace): workspace containing the variable
        varname (str): name of the variable
        fit_res (ROOT.RooFitResult, optional): If the desired variable is a
            ROOT.RooFormulaVar, the error depends on other values and they are
            taken from a fit result. If no fit result is passed in such a case
            the error will be set to -1

    Returns:
        tuple: First element is the value of the variable, second element is the
            error on the variable
    """
    var = get_var(wsp, varname)
    if hasattr(var, 'getError'):
        return var.getVal(), var.getError()
    if hasattr(var, 'getPropagatedError') and fit_res is not None:
        return var.getVal(), var.getPropagatedError(fit_res)
    return var.getVal(), -1


def get_chi2_ndf(fit_res, frame, pdfname, histname):
    """
    Get the chi2 value und the number degrees of freedom between a pdf and a
    histogram for a given fit result

    Args:
        fit_res (ROOT.RooFitResult): Fit result corresponding to the pdf.
            Necessary to obtain the number of floating parameters in the fit
        frame (ROOT.RooPlot): Plot containing the pdf and the histogram for
            which the chi2 should be calculated
        pdfname (str): Name of the pdf as it can be found on the frame
        histname (str): name of the datahist as it can be found on the frame

    Returns:
        tuple: chi2 and ndf as floats
    """
    n_float_pars = fit_res.floatParsFinal().getSize()
    n_bins = frame.GetNbinsX()
    ndf = n_bins - n_float_pars

    # reduced chi square between pdf and histogram assuming n_float_pars free
    # parameters in the fit
    chi2_ndf = frame.chiSquare(pdfname, histname, n_float_pars)
    chi2 = chi2_ndf * ndf

    logging.debug('fit results have {} floating parameters and histogram uses '
                  '{} bins. Reduced chi2 from frame is {:.2f}'
                  .format(n_float_pars, n_bins, chi2_ndf))

    return chi2, ndf


def get_args(rooarglist):
    """
    Get all elements from a ROOT.RooArgList.

    In principle it is possible to iterate over a RooArgList, however, it does
    not stop at the last element, resulting in a possible null-pointer
    dereference. This is just a helper function to circumvent this

    Args:
        rooarglist (ROOT.RooArgList)

    Returns:
        a generator expression yielding all elements of the passed RooArgList
    """
    n_elem = len(rooarglist) # Why ever this works
    return (rooarglist[i] for i in xrange(n_elem))


def get_corr_matrix(fit_result):
    """
    Get the correlation matrix from the passed fit_result as a numpy 2D array

    Args:
        fit_result (ROOT.RooFitResult): Fit result for which the correlation
            matrix of the free parameters should be obtained

    Returns:
        pandas.DataFrame: A dataframe containing the correlation matrix and the
            names of the variables as columns and indices (to facilitate
            plotting via seaborn)
    """
    # corr matrix of free params
    corr_matrix = get_np_from_tmatrix(fit_result.correlationMatrix())
    float_pars = [p.GetName() for p in get_args(fit_result.floatParsFinal())]

    return pd.DataFrame(corr_matrix, index=float_pars, columns=float_pars)
