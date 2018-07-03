#!/usr/bin/env python
"""
Module to do some symbolic calculations
"""

import sympy as sp

def _cov_matrix(params):
    """
    Get a symbolic covariance matrix in terms of the correlation coefficients
    and the parameter uncertainties
    """
    n_params = len(params)
    sigma = lambda par: 'sigma_' + par.name
    rho = lambda par1, par2: 'rho_' + par1.name + par2.name

    def fill_cov(irow, icol):
        """
        Helper function that creates the covariance symbols on the way
        """
        if irow == icol:
            return sp.symbols(sigma(params[irow]))**2
        else:
            # Since this will be symmetric, use corresponding notation
            # and switch indices to not create too much symbols
            if irow > icol:
                irow, icol = icol, irow
            rho12 = sp.symbols(rho(params[irow], params[icol]))
            sigma1 = sp.symbols(sigma(params[irow]))
            sigma2 = sp.symbols(sigma(params[icol]))
            return rho12 * sigma1 * sigma2

    return sp.Matrix(n_params, n_params, fill_cov)


def func_cov(func, free_params):
    """
    Calculate the the covariance matrix of the function assuming free_params
    have uncertainties

    Args:
        func (sympy expression): Function containing all the free parameters
        free_params (list of sympy symbols): The free parameters

    Returns:
        sympy expression: The covariance expression assuming known correlation
            coefficients and uncertainties for the free_params
    """
    func_m = sp.Matrix([func])
    jaco = func_m.jacobian(sp.Matrix(free_params))
    cov = _cov_matrix(free_params)

    return jaco.T.dot(cov.dot(jaco))


