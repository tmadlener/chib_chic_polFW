#!/usr/bin/env python
"""
Module that mainly holds constants values for nicer looking plots
"""

FIX_RANGES = {
    'CBmass1': [3.505, 3.5075],
    'CBmass2': [3.55, 3.5525],
    'CBsigma1': [0.004, 0.006],
    'CBsigma2': [0.0055, 0.009],
    # 'lambda_bkg': [-3.0, -1.7],
    'mu_bkg': [3.2, 3.25],
    'sigma_bkg': [0.05, 0.13],
    'r_chic0_chic1': [0.01, 0.06],
}

YLABELS = {
    'CBmass1': '#mu^{#chi_{c1}} / GeV',
    'CBmass2': '#mu^{#chi_{c2}} / GeV',
    'CBsigma1': '#sigma^{#chi_{c1}} / GeV',
    'CBsigma2': '#sigma^{#chi_{c2}} / GeV',
    'mu_bkg': '#mu^{Bkg} / GeV',
    'sigma_bkg': '#sigma^{Bkg} / GeV',
    'lambda_bkg': '#minus 1 / #lambda^{Bkg} / GeV^{-1}',
    'Nchic1': 'N^{#chi_{c1}}',
    'Nchic2': 'N^{#chi_{c2}}',
    'Nchic0': 'N^{#chi_{c0}}',
    'Nbkg': 'N^{Bkg}',
    'r_chic0_chic1': 'N^{#chi_{c0}} / N^{#chi_{c1}}',
    'r_chic2_chic1': 'N^{#chi_{c2}} / N^{#chi_{c1}}',
    'dlth': '#Delta#lambda_{#vartheta}',
    'lth': '#lambda_{#vartheta}(#chi_{c1})',
    'lph': '#lambda_{#varphi}(#chi_{c1})',
    'lth2': '#lambda_{#vartheta}(#chi_{c2})',
    'lph2': '#lambda_{#varphi}(#chi_{c2})',
    'dlph': '#Delta#lambda_{#varphi}',
    'ltilde': '#tilde{#lambda}(#chi_{c1})',
    'ltilde2': '#tilde{#lambda}(#chi_{c2})',
    'dltilde': '#Delta#tilde{#lambda}',
    'dkappa': '#Delta#kappa',
    'chicMass': 'M(J/#psi#gamma) / GeV',
    'JpsiPt': 'p_{T}^{J/#psi}',
    'costh_HX_fold': '|cos#vartheta^{HX}|',
    'phi_HX_fold': '#varphi^{HX}'
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
    'Nbkg': r'$N^{\text{Bkg}}$',
    'r_chic0_chic1': r'$N^{\chi_{c0}} / N^{\chi_{c1}}$',
    'dlth': r'$\Delta\lambda_{\vartheta}$',
    'lth': r'$\lambda_{\vartheta}(\chi_{c1})$'
}


CTH_LAB = '|cos#vartheta^{HX}|'
CTH_RAN = [0, 1]
PHI_LAB = '#varphi^{HX}_{fold}'
PHI_RAN = [0, 90]

CTH_PLOT = {'xRange': CTH_RAN, 'xLabel': CTH_LAB}
PHI_PLOT = {'xRange': PHI_RAN, 'xLabel': PHI_LAB}

VAR_PLOT = {'phi': PHI_PLOT, 'costh': CTH_PLOT}
