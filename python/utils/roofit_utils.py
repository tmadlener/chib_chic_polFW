#!/usr/bin/env python
"""
Module to facilitate some handling of RooFit objects
"""
import logging
logging.basicConfig(level=logging.WARNING,
                    format='%(levelname)s - %(funcName)s: %(message)s')
import ROOT as r
import numpy as np
import pandas as pd

from utils.misc_helpers import get_np_from_tmatrix, get_bin_centers, fmt_float


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


def all_vals(var):
    """
    Get the value, minimum value and maximum value the passed variable

    Args:
        var (ROOT.RooRealVar): The variable

    Returns:
        tuple of float: the current value, the minimum value and the maximum
            value
    """
    return var.getVal(), var.getMin(), var.getMax()


def get_chi2_ndf(frame, pdfname, histname, fit_res=None):
    """
    Get the chi2 value und the number degrees of freedom between a pdf and a
    histogram for a given fit result

    Args:
        frame (ROOT.RooPlot): Plot containing the pdf and the histogram for
            which the chi2 should be calculated
        pdfname (str): Name of the pdf as it can be found on the frame
        histname (str): name of the datahist as it can be found on the frame
        fit_res (ROOT.RooFitResult, optional): Fit result corresponding to the
            pdf. Necessary to obtain the number of floating parameters in the
            fit. If None (default) only the chi2 will be calculated and ndf will
            be set to the number of bins
    Returns:
        tuple: chi2 and ndf as floats
    """
    n_bins = frame.GetNbinsX()
    if fit_res is not None:
        n_float_pars = fit_res.floatParsFinal().getSize()
        logging.debug('fit results have {} floating parameters'
                      .format(n_float_pars))
    else:
        n_float_pars = 0
        logging.debug('No fit result passed')

    ndf = n_bins - n_float_pars

    # reduced chi square between pdf and histogram assuming n_float_pars free
    # parameters in the fit
    chi2_ndf = frame.chiSquare(pdfname, histname, n_float_pars)
    chi2 = chi2_ndf * ndf

    logging.debug('Reduced chi2 from frame is {:.2f}. Have {} bins in histogram '
                  'and {} floating parameters -> chi2 / ndf = {:.2f} / {}'
                  .format(chi2_ndf, n_bins, n_float_pars, chi2, ndf))

    return chi2, ndf


def calc_info_crit(fit_res, min_nll, n_events=None):
    """
    Calculate the BIC (Bayesian Information Criterion) or the AIC (Akaike
    Information Criterion) from the passed fit result.

    Args:
        fit_res (ROOT.RooFitResult): Fit result of a maximum likelihood fit
        min_nll (float): The value of the negative log likelihood at its
            minimum. If None is passed, the minimum stored in the fit result
            will be used.
        n_events (int or None): If an int is passed it is assumed that this is
            the number of events in the sample and the BIC is calculated, if
            None is passed the AIC is calculated

    See also:
        https://en.wikipedia.org/wiki/Bayesian_information_criterion
        https://en.wikipedia.org/wiki/Akaike_information_criterion
    """
    if min_nll is None:
        min_nll = fit_res.minNll()
    n_pars = fit_res.floatParsFinal().getSize()
    if n_events is not None:
        logging.debug('Calculating BIC with -log(L) = {:.0f}, k = {}, n = {}'
                      .format(min_nll, n_pars, n_events))
        return np.log(n_events) * n_pars + 2 * min_nll

    logging.debug('Calculating AIC with -log(L) = {:.0f}, k = {}'
                  .format(min_nll, n_pars))
    return 2 * n_pars + 2 * min_nll


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
    # symmetric (HESSE) uncertainties, just in case we need them later
    err_sym = []

    wvar = get_var(wsp, var)

    for snap in snapshots:
        wsp.loadSnapshot(snap)

        central.append(wvar.getVal())
        err_low.append(-wvar.getErrorLo())
        err_high.append(wvar.getErrorHi())
        err_sym.append(wvar.getError())

    # check if all the asymmetric uncertainties are present
    if np.any(np.array(err_low) == 0) or np.any(np.array(err_high) == 0):
        logging.warning('Some of the asymmetric (MINOS) uncertainties for {} are'
                        ' not properly determined, switching to symmetric '
                        '(HESSE) uncertainties'.format(var))
        err_low = err_high = err_sym

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


def set_var(wsp, varname, val, create=False, err=None):
    """
    Set the passed value to the variable in the workspace. If the variable is
    not already present, create it first.
    """
    var = get_var(wsp, varname)
    if var is None:
        if not create:
            logging.error('Variable \'{}\' is not present in workspace'
                          .format(varname))
            return
        else:
            logging.debug('Variable \'{}\' not already present in workspace. '
                          'Creating it.'.format(varname))
            wsp.factory("{}[{}]".format(varname, val))
            var = get_var(wsp, varname)
    else:
        var.setVal(val)

    # At this point we definitely have the variable
    if err is not None:
        var.setError(err)


def fix_params(wsp, param_vals):
    """
    Fix the parameters in the workspace to given or current values

    Args:
        wsp (ROOT.RooWorkspace): workspace containing all the variables
        param_vals (list of tuples): tuples where first element is the name
            of the parameter and the second is the value to which it should
            be fixed. If the second parameter is None, it will be fixed to
            the current value in the workspace
    """
    for par, val in param_vals:
        if val is None:
            val = get_var(wsp, par).getVal()
            debug_msg = 'Fixing {} to current value in workspace: {}'
        else:
            debug_msg = 'Fixing {} to {}'

        logging.info(debug_msg.format(par, val))
        get_var(wsp, par).setVal(val)
        get_var(wsp, par).setConstant(True)


def release_params(wsp, param_names):
    """
    Release the parameters in the workspace

    Args:
        wsp (ROOT.RooWorkspace): workspace containing all the variables
        param_names (list of strings): Names of the parameters in the
            workspace to release
    """
    for par in param_names:
        logging.debug('Releasing variable \'{}\''.format(par))
        get_var(wsp, par).setConstant(False)


def try_factory(wsp, expr):
    """Try to run the expression through the RooWorkspace.factory"""
    logging.debug(expr)
    obj = wsp.factory(expr)
    if not obj:
        logging.error('\'{}\' failed to create the an object in the workspace'
                      .format(expr))
        return False
    return True


def param_str(wsp, param):
    """
    Get a string representation of the values of the parameter in the workspace.
    Uses misc_utils.fmt_float for formatting.

    Args:
        wsp (ROOT.RooWorkspace): The workspace containing the parameter
        param (str): The name of the parameter

    Returns:
        str (string that can be used in latex math environments)
    """
    var = get_var(wsp, param)
    if var.isConstant():
        return r'${}$ (fixed)'.format(fmt_float(var.getVal()))
    val = var.getVal()
    var_exp = np.floor(np.log10(val)) if val > 0 else np.floor(np.log10(-val))

    return r'${} \pm {}$'.format(fmt_float(val),
                                 fmt_float(var.getError(), var_exp.astype(int)))


def eval_pdf(pdf, var, values):
    """
    Evaluate a RooAbsPdf at all passed values

    NOTE:
    This has to be done in a bit a cumbersome fashion, since evaluate is not
    directly exposed and even if it was it does not take an argument at which
    the pdf should be evaluated.

    This function does not take into account any normalization! This is a
    RooFit "oddity" that actually makes sense (although I never can keep the
    reasoning for it in memory for longer than a few hours)

    In order to be able to compare the numbers that are returned from this
    function with the return value for another RooAbsPdf they have to be
    multiplied by a factor that relates the pdfs such that they are correctly
    normalized against each other. The easiest way for the current purpose is
    to use the respective yields

    Args:
        pdf (ROOT.RooAbsPdf): The pdf function that should be evaluated
        var (ROOT.RooAbsReal): The variable for which the pdf should be
            evaluated at the passed values
        values (np.array of floats): All values for which the pdf should be
            evaluated.

    Returns:
        np.array: The values of the passed pdf at the passed variable values
    """
    logging.debug('Evaluating pdf {} along variable {} at {} values'.
                  format(pdf.GetName(), var, len(values)))
    vals = np.ones_like(values)
    # Necessary to have the resulting pdf normalized to 1 on the range of the
    # variable. Only with this do the relative normalizations using the yields
    # make sense and give the expected results. Keeping it outside the loop to
    # avoid creating it everytime
    var_set = r.RooArgSet(var)
    for idx, val in enumerate(values):
        var.setVal(val)
        vals[idx] = pdf.getVal(var_set)

    return vals
