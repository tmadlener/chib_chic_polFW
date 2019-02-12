#!/usr/bin/env python
"""
Module to facilitate some handling of RooFit objects
"""
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')
import ROOT as r
import numpy as np
import pandas as pd

from utils.misc_helpers import get_np_from_tmatrix, get_bin_centers


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


def _get_var_vals(wsp, var, snapshots):
    """
    Get the variable (and its uncertainties for all snapshots)
    """
    central = []
    err_low = []
    err_high = []
    wvar = get_var(wsp, var)

    for snap in snapshots:
        wsp.loadSnapshot(snap)

        central.append(wvar.getVal())
        err_low.append(-wvar.getErrorLo())
        err_high.append(wvar.getErrorHi())

    return np.array(central), np.array(err_low), np.array(err_high)


def get_var_graph(wsp, snap_base, var, n_bins, binning=None, bin_means=None,
                  fit_res_base=''):
    """
    Get the graph of a variable for all snapshots results matching fit_res_base
    """
    snapshots = [snap_base.format(i) for i in xrange(n_bins)]
    dependent = isinstance(get_var(wsp, var), r.RooFormulaVar)
    if dependent and not fit_res_base:
        logging.error('{} is a dependent variable but no base name for the fit '
                      'results is passed. Cannot calculate uncertainties'
                      .format(var))
        return None

    if dependent:
        fitresults = [fit_res_base.format(i) for i in xrange(n_bins)]
        vals = []
        errs = []
        for ibin, fitres in enumerate(fitresults):
            wsp.loadSnapshot(snapshots[ibin])
            val_err = get_var_err(wsp, var, wsp.genobj(fitres))
            vals.append(val_err[0])
            errs.append(val_err[1])

        vals = np.array(vals)
        err_lo = np.array(errs)
        err_hi = np.array(errs)
    else:
        vals, err_lo, err_hi = _get_var_vals(wsp, var, snapshots)

    # bin_means are necessary and are either determined from the binning or
    # simply chosen as the indices
    if bin_means is None:
        if binning is not None:
            bin_means = get_bin_centers(binning)
        else:
            bin_means = np.linspace(0, n_bins - 1, n_bins)

    if binning is not None:
        x_lo = np.array([bin_means[i] - binning[i] for i in xrange(n_bins)])
        x_hi = np.array([binning[i+1] - bin_means[i] for i in xrange(n_bins)])
    else:
        x_lo, x_hi = np.zeros(n_bins), np.zeros(n_bins)

    return r.TGraphAsymmErrors(n_bins, bin_means, vals, x_lo, x_hi,
                               err_lo, err_hi)
