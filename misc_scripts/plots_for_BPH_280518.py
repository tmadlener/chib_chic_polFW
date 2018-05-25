#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import ROOT as r

from utils.data_handling import (
    get_dataframe, apply_selections, create_histogram
)
from utils.selection_functions import (
    trigger_sel, photon_sel, fiducial_muon_sel, flat_pt, jpsi_kin_sel,
    get_n_events, chic_state_sel
)
from utils.plot_helpers import mkplot, default_attributes

def costh_phi(df, frame, gen=False):
    """
    Get Abs(costh) and folded phi from dataframe for a given reference frame
    """
    add_gen = 'gen_' if gen else ''
    costh = '{}costh_{}'.format(add_gen, frame)
    phi = '{}phi_{}_fold'.format(add_gen, frame)

    return np.array([df[costh].abs(), df[phi]]).T


costh_phi_settings = {
    'costh': (8, 0, 1),
    'phi': (10, 0, 90),
    'costh_phi': (8, 0, 1, 10, 0, 90)
}


def get_var_histo(df, var, frame, gen=False):
    """
    Get the histogram for a given state and a given variable in a given frame
    """
    hist_sett = costh_phi_settings[var]
    ctp = costh_phi(df, frame, gen) # need costh-phi distribution in any case
    if var == 'costh_phi':
        return create_histogram(ctp, costh_phi_settings[var],
                                x_axis='|cos#theta^{{{}}}|'.format(frame),
                                y_axis='#phi^{{{}}}_{fold}'.format(frame))
    elif var == 'phi':
        return create_histogram(ctp[:,1], costh_phi_settings[var],
                                x_axis='#phi^{{{}}}_fold'.format(frame))
    elif var == 'costh':
        return create_histogram(ctp[:,0], costh_phi_settings[var],
                                x_axis='|cos#theta^{{{}}}|'.format(frame))
    else:
        logging.error('var has to be one of \'costh_phi\', \'costh\' or \'phi\' but was %s', var)


def get_chic2_chic1_ratio(df, var, frame, gen=False):
    """
    Get the chic2 / chic1 ratio for a given variable in a given frame
    """
    chi1_data = apply_selections(df, lambda df: chic_state_sel(df, 'chic1') )
    chi2_data = apply_selections(df, lambda df: chic_state_sel(df, 'chic2') )

    chi1_hist = get_var_histo(chi1_data, var, frame, gen)
    chi2_hist = get_var_histo(chi2_data, var, frame, gen)

    chi2_hist.Divide(chi1_hist)
    chi2_hist.SetYTitle('#chi_{c2} / #chi_{c1}')

    return chi2_hist


gen_data = get_dataframe('/afs/hephy.at/data/tmadlener01/ChicPol/Chic2012/InputFiles/NewMCGen2012/flat_tuples/pGunRunI-chic-gen_flat_new_tupling.root')
reco_data = get_dataframe('/afs/hephy.at/data/tmadlener01/ChicPol/Chic2012/InputFiles/NewMCGen2012/flat_tuples/pGunRunI-chic-lowPt_merged_new_tupling.root')

gen_selection = (
    lambda df: photon_sel(df, flat_pt(0, 1.5), gen=True),
    lambda df: fiducial_muon_sel(df, gen=True),
    lambda df: jpsi_kin_sel(df, 8, 20, 1.2, gen=True)
)

# reco_selection = gen_selection + (trigger_sel,)
reco_selection = gen_selection


print('Number of events after selection: gen: {}, reco: {}'
      .format(get_n_events(gen_data, gen_selection),
              get_n_events(reco_data, reco_selection)))

gen_ratio = get_chic2_chic1_ratio(apply_selections(gen_data, gen_selection),'costh', 'HX', True)
reco_ratio = get_chic2_chic1_ratio(apply_selections(reco_data, reco_selection),'costh', 'HX', True)

leg = r.TLegend(0.15, 0.15, 0.25, 0.25)
leg.SetBorderSize(0)

can = r.TCanvas('c', 'c', 600, 600)
can.cd()
ratio_pad = r.TPad('ratio_pad', '', 0, 0.3, 1, 1)
r.SetOwnership(ratio_pad, False)
ratio_pad.Draw()
ratio_pad.cd()

ratio_pad = mkplot([gen_ratio, reco_ratio], yRange=[0, None], drawOpt='E1', legOpt='PLE',
                   attr=default_attributes(size=1.0, linewidth=1), can=ratio_pad,
                   leg=leg, legEntries=['gen', 'reco'], yLabel='N_{#chi_{c2}} / N_{#chi_{c1}}')

can.cd()


reco_gen_dr = reco_ratio.Clone()
reco_gen_dr.Divide(gen_ratio)


dr_pad = r.TPad('double_ratio_pad', '', 0, 0, 1, 0.3)
r.SetOwnership(dr_pad, False)
dr_pad.Draw()

reco_gen_dr.GetYaxis().SetTitleSize(0.08)
reco_gen_dr.GetXaxis().SetTitleSize(0.08)
reco_gen_dr.GetYaxis().SetLabelSize(0.08)
reco_gen_dr.GetXaxis().SetLabelSize(0.08)
reco_gen_dr.GetYaxis().SetTitleOffset(0.5)
reco_gen_dr.SetXTitle('')

dr_pad = mkplot(reco_gen_dr, yRange=[0.8, 1.2], drawOpt='E1', can=dr_pad,
                attr=[{'color': 1, 'marker': 20, 'size': 1.0}], yLabel='reco / gen')

can.SaveAs('gen_reco_noTrigger_MC_costh_ratio_with_double_ratio.pdf')
