#!/USSR/bin/env python
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

    return jaco * cov * jaco.T


def _lambda_theta(j, simplified=True):
    """
    Get the lambda_theta formula for the given J state,
    See PRD 83, 096001 (2011) for details
    """
    h2, g2, g3 = 0, 0, 0 # E1-approximation
    delta_1 = sp.S(2) / 5 * h2**2
    delta_2 = 2 * g2**2 + sp.S(5) / 7 * g3**2


    if simplified:
        R, R2 = sp.symbols('R, R2')
    else:
        bm1, bp1, bm2, bp2 = sp.symbols('b_0, b_{-1}, b_{+1}, b_{-2}, b_{+2}')
        R = abs(bp1)**2 + abs(bm1)**2
        R2 = abs(bp2)**2 + abs(bm2)**2

    b02 = abs(1 - R - R2)

    D1 = (2 * (1 + delta_1) * b02 + (3 - delta_1) * R) / (1 - 3*delta_1)
    D2 = (2 * (5 - delta_2) * b02 + (9 - delta_2) * R + 2 * (3 + delta_2) * R2) /\
         (1 - delta_2)

    if j == 1:
        # For j = 1, there is no Jz = 2 component -> remove it
        return (1 / D1 * (2 * b02 - R)).subs({R2: 0})

    if j == 2:
        return - 3 / D2 * (2 * b02 + R - 2 * R2)


def lth_2(**kwargs):
    """
    Get lambda theta for the J=2 state

    Keyword Args:
        R1 (float, or sympy.Symbol): The fraction of |Jz|=1,
        R2 (flaot, or sympy.Symbol): The fraction of |Jz|=2
        simplified (Boolean, optional, default=True): If false return the
            formula using actual b_m parameters instead of simplified version
            with Rs

    Returns:
        A sympy expression or a number depending on the inputs
    """
    lth = _lambda_theta(2, kwargs.pop('simplified', True))
    R1 = kwargs.pop('R1', None)
    R2 = kwargs.pop('R2', None)

    if R1 is not None:
        lth = lth.subs({'R': sp.S(R1)})

    if R2 is not None:
        lth = lth.subs({'R2': sp.S(R2)})

    return lth


def lth_1(**kwargs):
    """
    Get lambda theta for the J=1 state

    Keyword Args:
        R (float, or sympy.Symbol): The fraction of |Jz|=1,
        simplified (Boolean, optional, default=True): If false return the
            formula using actual b_m parameters instead of simplified version
            with Rs

    Returns:
        A sympy expression or a number depending on the inputs
    """
    lth = _lambda_theta(1, kwargs.pop('simplified', True))
    R = kwargs.pop('R', None)

    if R is not None:
        lth = lth.subs({'R': sp.S(R)})
    return lth
