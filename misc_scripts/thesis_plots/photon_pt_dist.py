#!/usr/bin/env python
"""
Script to make a plot of the photon pt distribution
"""

import os
import numpy

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.data_handling import get_dataframe, apply_selections
from utils.selection_functions import (
    single_muon_sel, select_bin, flat_pt, trigger_sel_, vtx_prob_sel_,
    prompt_sel_, collect_requirements, photon_sel_
)
from utils.hist_utils import hist1d
from utils.plot_helpers import mkplot
from utils.setup_plot_style import set_TDR_style

RDIR = os.path.join(os.environ['MY_DATA_DIR'], 'ChicPol', 'Chic2012', 'Results')
DATADIR = os.path.join(os.environ['MY_DATA_DIR'], 'ChicPol', 'Chic2012',
                       'InputFiles', 'Reprocessing2012Data', 'flat_tuples')
INFILE = os.path.join(DATADIR,
                      'rootuple-chic-8TeV-BCD-Reprocessing_tupling_w_mmugammma_w_lt_sig.root')
OUTDIR = os.path.join(RDIR, 'thesis_plots')

SELECTIONS = (
    select_bin('chicMass', 3.2, 3.75),
    select_bin('JpsiPt', 18, 30),
    select_bin('abs(JpsiRap)', 0, 1.2),
    single_muon_sel(flat_pt(3.5, 1.6)),
    trigger_sel_('Dimuon8_Jpsi'),
    vtx_prob_sel_(),
    prompt_sel_(),
    photon_sel_(flat_pt(0, 1.5))
)


def make_photon_pt_dist_plot(data):
    """Make photon pt distribution plot"""
    hist = hist1d(data.photonPt, min=0., max=10, nbins=50, log=False)

    can = mkplot(hist, drawOpt='PE',
                 # logx=True,
                 logy=True,
                 xRange=[0.0, 10], yRange=[0.5, 0.7e4],
                 xLabel='p_{T}^{#gamma} [GeV]',
                 yLabel='Events / 0.2 GeV',
                 attr=[{'marker': 20, 'size': 1.0, 'color': 1}])

    return can


def main():
    """Main"""
    set_TDR_style()
    r.gStyle.SetPadRightMargin(r.gStyle.GetPadRightMargin() + 0.01)

    data = apply_selections(get_dataframe(INFILE,
                                          columns=collect_requirements(SELECTIONS)),
                            SELECTIONS)

    can = make_photon_pt_dist_plot(data)
    can.SaveAs(os.path.join(OUTDIR, 'photon_pt_dist.pdf'))



if __name__ == '__main__':
    main()
