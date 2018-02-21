#!/usr/bin/env python
"""
Module to facilitate some handling of RooFit objects
"""

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
