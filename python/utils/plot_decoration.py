#!/usr/bin/env python
"""
Module that mainly holds constants values for nicer looking plots
"""

FIX_RANGES = {
    'CBmass1': [3.505, 3.5075],
    'CBmass2': [3.546, 3.556],
    'CBsigma1': [0.004, 0.006],
    'CBsigma2': [0.0055, 0.009],
    'lambda_bkg': [-3.0, -1.7],
    'mu_bkg': [3.2, 3.25],
    'sigma_bkg': [0.05, 0.13],
}

YLABELS = {
    'CBmass1': '#mu^{#chi_{c1}}',
    'CBmass2': '#mu^{#chi_{c2}}',
    'CBsigma1': '#sigma^{#chi_{c1}}',
    'CBsigma2': '#sigma^{#chi_{c2}}',
    'mu_bkg': '#mu^{Bkg}',
    'sigma_bkg': '#sigma^{Bkg}',
    'lambda_bkg': '#lambda^{Bkg}',
    'Nchic1': 'N^{#chi_{c1}}',
    'Nchic2': 'N^{#chi_{c2}}',
    'Nchic0': 'N^{#chi_{c0}}',
    'Nbkg': 'N^{Bkg}'
}

PLOT_LABELS_LATEX = {
    'CBmass1': r'$\mu^{\chi_{c1}}$',
    'CBmass2': r'$\mu^{\chi_{c2}}$',
    'CBsigma1': r'$\sigma^{\chi_{c1}}$',
    'CBsigma2': r'$\sigma^{\chi_{c2}}$',
    'mu_bkg': r'$\mu^{\text{Bkg}}$',
    'sigma_bkg': r'$\sigma^{\text{Bkg}}$',
    'lambda_bkg': r'$\lambda^{\text{Bkg}}$',
    'Nchic1': r'$N^{\chi_{c1}}$',
    'Nchic2': r'$N^{\chi_{c2}}$',
    'Nchic0': r'$N^{\chi_{c0}}$',
    'Nbkg': r'$N^{\text{Bkg}}$'
}
