#!/usr/bin/env python
"""
Script to make plots showing acceptance effects
"""

import os
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.data_handling import get_dataframe, apply_selections
from utils.selection_functions import single_muon_sel, select_bin, flat_pt
from utils.hist_utils import hist2d, hist1d
from utils.plot_helpers import mkplot, setup_legend, default_attributes
from utils.setup_plot_style import set_TDR_style
from utils.plot_decoration import YLABELS, CTH_LAB
from utils.misc_helpers import cond_mkdir, _get_var


RDIR = os.path.join(os.environ['MY_DATA_DIR'], 'ChicPol', 'Chic2012', 'Results')
DATADIR = os.path.join(RDIR, 'input_files')
INFILE = os.path.join(DATADIR, 'presel_8_30_single_mu_3p5_vtx_prob_0p01_lt_sig_2p5.root')
OUTDIR = os.path.join(RDIR, 'thesis_plots')

SELECTIONS = (
    select_bin('chicMass', 3.2, 3.75),
    select_bin('JpsiPt', 12, 18),
    select_bin('abs(JpsiRap)', 0, 1.2),
    single_muon_sel(flat_pt(3.5, 1.6)),
)

def make_costh_phi_plot(data, frame):
    """Make costh-phi histogram"""
    hist = hist2d(_get_var(data, 'costh_{}_fold'.format(frame), np.abs),
                  _get_var(data, 'phi_{}_fold'.format(frame)),
                  nbinsx=24, minx=0, maxx=1,
                  nbinsy=32, miny=0, maxy=90)

    cth_lab = '|cos#vartheta^{{{}}}|'.format(frame)
    phi_lab = '#varphi^{{{}}}_{{fold}}'.format(frame)

    can = mkplot(hist, drawOpt='colz', xRange=[0, 1], xLabel=cth_lab,
                 yRange=[0, 90], yLabel=phi_lab)
    can.pltables[0].SetNdivisions(505, 'X')
    can.Update()
    return can


def make_costh_mu_pt_plot(data):
    """Make a comparison of costh distributions for different muon pt cuts"""
    pt_cuts = [3.5, 4, 4.5, 5, 5.5]
    mu_pt_labels = ['p_{{T}}^{{#mu}} > {:.1f} GeV'.format(p) for p in pt_cuts]
    hists = [
        hist1d(apply_selections(data, single_muon_sel(flat_pt(p, 1.6))).costh_HX_fold.abs(),
               nbins=24, min=0, max=1) for p in pt_cuts
    ]
    [h.Scale(1.0 / h.GetBinContent(1)) for h in hists]

    can = mkplot(hists, drawOpt='PE', xRange=[0, 1], xLabel=CTH_LAB,
                 attr=default_attributes(size=1.0, width=2, open_markers=False),
                 legPos=(0.65, 0.50, 0.8, 0.92), legEntries=mu_pt_labels,
                 yLabel='normalized to {} = 0 [a.u.]'.format(CTH_LAB))
    can.pltables[0].SetNdivisions(505, 'X')
    can.attached_tobjects[0].SetTextSize(0.05)

    return can


def main():
    """Main"""
    set_TDR_style()
    cond_mkdir(OUTDIR)

    data = apply_selections(get_dataframe(INFILE), SELECTIONS)

    r.gStyle.SetPadRightMargin(0.129)
    r.gStyle.SetPadLeftMargin(r.gStyle.GetPadLeftMargin() - 0.007)
    can = make_costh_phi_plot(data, 'HX')
    can.SaveAs(os.path.join(OUTDIR, 'costh_phi_fold_HX_pt_12_18_all.pdf'))

    can = make_costh_phi_plot(data, 'CS')
    can.SaveAs(os.path.join(OUTDIR, 'costh_phi_fold_CS_pt_12_18_all.pdf'))

    set_TDR_style()
    can = make_costh_mu_pt_plot(data)
    can.SaveAs(os.path.join(OUTDIR, 'costh_HX_comp_mu_pt_jpsipt_12_18_all.pdf'))


if __name__ == '__main__':
    main()
