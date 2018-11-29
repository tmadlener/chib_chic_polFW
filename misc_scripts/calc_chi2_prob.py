#!/usr/bin/env python
"""
Script to calculate the chi2 probability of data to Toy MC predictions
"""
import re
import json
import numpy as np
import sympy as sp

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()
r.gSystem.Load('$CHIB_CHIC_POLFW_DIR/general/bin/CustomRooFitFunctions')

from scipy.stats import norm, chi2
from scipy.optimize import minimize
from itertools import product

from utils.symbolic import lth_1, lth_2
from utils.roofit_utils import get_var_err, get_var
from utils.hist_utils import divide
from utils.graph_utils import get_errors
from utils.plot_helpers import mkplot, setup_legend, default_colors
from utils.misc_helpers import unique_w_key


import imp
ru = imp.load_source('root_utils', "/afs/hephy.at/user/j/jnecker/code/root_makros/python/root_utils.py")



def get_fit_graph(wsp, costh_bins, costh_means):
    """
    Get the graph of the ratio from the workspace
    """
    ratio_central = []
    ratio_err = []

    n_bins = len(costh_bins)

    for i in xrange(n_bins):
        wsp.loadSnapshot('snap_costh_bin_{}'.format(i))
        central, err = get_var_err(wsp, 'ratio_chib2_chib1',
                                   wsp.genobj('fit_res_costh_bin_{}'.format(i)))
        ratio_central.append(central)
        ratio_err.append(err)

    ratio = np.array(ratio_central)
    err = np.array(ratio_err)

    c_lo = np.array([costh_means[i] - costh_bins[i][0] for i in xrange(n_bins)])
    c_hi = np.array([costh_bins[i][1] - costh_means[i] for i in xrange(n_bins)])

    return r.TGraphAsymmErrors(n_bins, np.array(costh_means), ratio,
                               c_lo, c_hi, err, err)


def calc_chi2(pred_hist, data_graph):
    """
    Calculate the chi2 between the predictor histogram and the data graph
    """
    # Need to select the bins without over- and underflow bin here
    pred = np.array([b for b in pred_hist][1:pred_hist.GetNbinsX() + 1])
    data_central = np.array(data_graph.GetY())
    _, _, data_err, _ = get_errors(data_graph)

    def chisquare((norm,)):
        return np.sum((( norm * pred - data_central) / data_err)**2)

    min_res = minimize(chisquare, np.array(data_central[0]), method='nelder-mead')
    fit_norm = min_res.x[0]

    return np.sum(((fit_norm * pred - data_central) / data_err)**2), fit_norm


def get_ratio_combs(chi1_hists, chi2_hists):
    """
    Get all the possible ratio combinations
    """
    return product(chi1_hists, chi2_hists)


def get_delta_lambda(chi1_name, chi2_name):
    """
    Get the delta lambda from the chi1 and the chi2 name
    """
    ratio_rgx = r'([0-9]o?[0-9]?)'
    chi1_rgx = r'.*_R_' + ratio_rgx
    chi2_rgx = r'.*_R1_' + ratio_rgx + '_R2_' + ratio_rgx

    lth1 = None
    lth2 = None

    m = re.search(chi1_rgx, chi1_name)
    if m:
        R = sp.S(m.group(1).replace('o', '/'))
        print R
        lth1 = lth_1(R=R)

    m = re.search(chi2_rgx, chi2_name)
    if m:
        R1 = sp.S(m.group(1).replace('o', '/'))
        R2 = sp.S(m.group(2).replace('o', '/'))
        print R1,R2
        lth2 = lth_2(R1=R1, R2=R2)

    return lth2 - lth1


def main(args):
    """Main"""
    
    with open(args.bininfo, 'r') as finfo:
        costh_info = json.load(finfo)

    fitfile = r.TFile.Open(args.fitres)
    wsp = fitfile.Get('ws_mass_fit')

    fit_graph=args.tgraph
    if fit_graph is None: fit_graph = get_fit_graph(wsp, costh_info['costh_bins'],
                              costh_info['costh_means'])
    print fit_graph

    frame = costh_info.get('frame','HX')

    histfile = r.TFile.Open(args.histfile)
    histlist = list(set(b.GetName() for b in histfile.GetListOfKeys()))

    ratio_combs = get_ratio_combs([h for h in histlist if 'chic1' in h or 'chib1' in h],
                                  [h for h in histlist if 'chic2' in h or 'chib2' in h])

    ratios = []
    for chi1_hist, chi2_hist in ratio_combs:
        h_chi1 = histfile.Get(chi1_hist)
        h_chi1.Scale(1 / h_chi1.Integral())
        h_chi2 = histfile.Get(chi2_hist)
        h_chi2.Scale(1 / h_chi2.Integral())
        ratio = divide(h_chi2, h_chi1)

        chisqu, fit_norm = calc_chi2(ratio, fit_graph)
        ratios.append((get_delta_lambda(chi1_hist, chi2_hist),
                       ratio, chisqu, fit_norm))
        
        ratio.Scale(fit_norm)

    ratios.sort(key=lambda rr: sp.N(rr[0]))
    # print before removing "duplicates"
    f=open(args.saveas+'_chi2.txt','w')
    for ratio in ratios:
        tmp_str='Delta lambda = {}, N = {:.3f}, chi2 = {:.2f}, p = {:.3e}, sigma = {:.1f}'.format(
            ratio[0], ratio[3], ratio[2],
            chi2.sf(ratio[2], fit_graph.GetN() - 1),
            norm.isf(chi2.sf(ratio[2], fit_graph.GetN() - 1)))
        print(tmp_str)
        f.write(tmp_str+'\n')
        


    # only retain one of the possible delta lambda = 0 results
    ratios = unique_w_key(ratios, lambda rr: rr[0])

    leg = setup_legend(0.12, 0.12, 0.5, 0.4)
    leg.SetNColumns(2)

    plot_toy = [sp.S(x) for x in ['-8/5', '-1', '-1/3', '0', '4/3']]

    can = mkplot([rr[1] for rr in ratios if rr[0] in plot_toy],
                 yRange=[0, 1], yLabel='#chi_{b2} / #chi_{b1}',
                 xLabel='|cos#vartheta^{'+frame+'}|',
                 leg=leg,
                 legEntries=['#Delta#lambda_{#theta} = ' + str(rr[0])
                             for rr in ratios if rr[0] in plot_toy])

    can = mkplot(fit_graph, can=can, leg=can.attached_tobjects[0],
                 legEntries=['data'],
                 attr=[{'color': 1, 'marker': 20, 'size': 2}], drawOpt='PE')
    
    #ru.add_cms_lumi(can, args.lumi)

    can.SaveAs(args.saveas+'.pdf')


    histos2plot=[]
    def_colors=default_colors()
    def_colors.append(r.kMagenta -1)
    j=0
    for i,rr in enumerate(ratios):
         if rr[0] in plot_toy:
            histos2plot.append(
                { 'histogram': rr[1],
                  'legend':'#Delta#lambda_{#vartheta} = ' + str(rr[0]),
                  'color': def_colors[j%len(def_colors)], 
                  'marker': 24+j, 
                  'markersize': 1.3,
                  'linewidth': 1})            
            j+=1
            if j==1:j+=1
    fit_graph.SetFillColor(0)
    fit_graph.SetFillStyle(0)

    histos2plot+=[{'histogram':fit_graph,
                  'legend': 'Data',
                  'color': 1, 
                  'marker': 25, 
                  'markersize': 2,
                  'linewidth': 2
                  }]

   
    ru.draw_histos1_overlayed( histos2plot,
        grid = False,
        drawlegend = True,
        legendposition='tl',
        legendcolumns=2,
        title = ';|cos#vartheta^{'+frame+'}|;Ratio #chi_{b2} / #chi_{b1}',
        saveas= args.saveas+'.pdf',
        width = 1000,
        height = 900,
        lumi=args.lumi,
        ymin=0,
        ymax=1.05,
        xmax=1,
        xmin=0,
        legendoption='LP',
        transparentlegend=True
    )


    saveto = r.TFile.Open(args.saveas+'.root','recreate')
    for (i,rr) in enumerate(ratios):
        eins=rr[0].p
        zwei='o{}'.format(rr[0].q) if rr[0].q !=1 else ''
        name='deltalambda_{}{}_h{}'.format(eins,zwei,i)
        print name
        h=rr[1]
        h.SetTitle('#Delta#lambda_{#theta} = ' + str(rr[0]))
        saveto.cd()
        h.Write(name, r.TObject.kWriteDelete)
    saveto.cd()
    fit_graph.Write('costhbinned_fitresult_graph', r.TObject.kWriteDelete)
    can.Write('complete_canvas')
    saveto.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to calculate the chi2 '
                                     'probability of data fit results to Toy MC'
                                     ' predictions')
    parser.add_argument('fitres', help='File containing the results of costh '
                        'binned mass fits')
    parser.add_argument('histfile', help='File containing the selected and '
                        'efficiency weighted costh histograms from Toy MC')
    parser.add_argument('bininfo', help='File containing the costh binning info '
                        '(as produced e.g. by the costh binnned mass fits)')
    parser.add_argument('--saveas', help='Output file name without ending, '
                        'used for pdf and root file containing hists',
                        default='simple_chi2_overview_toy_data_comp')
    parser.add_argument('--tgraph', help='TGraph for the chi2 calculation.',
                        default=None)

    clargs = parser.parse_args()
    main(clargs)

