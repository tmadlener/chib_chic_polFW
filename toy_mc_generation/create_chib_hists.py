# Plot comparing different relative polarization scenarios to costh binned fit

import ROOT as r
import numpy as np
import re, json

from utils.data_handling import get_dataframe, apply_selections
import utils.selection_functions as sf
from utils.hist_utils import create_histogram, divide
from utils.plot_helpers import mkplot
import imp
ru = imp.load_source('root_utils', "/afs/hephy.at/user/j/jnecker/code/root_makros/python/root_utils.py")


if __name__=='__main__':
    import argparse, inspect
    
    parser = argparse.ArgumentParser(description='Script to run data to MC comparison for chib')
    parser.add_argument('fitres', help='File containing the results of costh '
                        'binned mass fits')
    parser.add_argument('bininfo', help='JSON file containing the costh binning info '
                        '(as produced e.g. by the costh binnned mass fits)')
    parser.add_argument('outfileprefix', help='Output file WITHOUT extension')
    parser.add_argument('configfile',help='chib fitting config file')
    parser.add_argument('--year', help='Data taking period', default=2016)
    parser.add_argument('--maxphotoneta', default=1.2)

    args = parser.parse_args()
    print "Called with arguments:"
    print (args)
    
    
    chi2_outfile = args.outfileprefix #'chib_fastmc_histograms_costhHX.root'
    jsonbinning = args.bininfo #'/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/4setting2017/pT_10_50_h/costh_binned_fit_10_50/fitresults-dimuon_mass_8p7_11p1-dimuon_pt_10_50_costhbins.json'
    costhbinnedfit_results = args.fitres #'/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/4setting2017/pT_10_50_h/costh_binned_fit_10_50/fitresults-dimuon_mass_8p7_11p1-dimuon_pt_10_50_costhbinned_fit.root'
    #outpath = '/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/4setting2017/pT_10_50_h/plots'
    #inpath = '/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/4setting2017/pT_10_50_h/costh_binned_fit_10_50'
    #jsonbinning = inpath+'/fitresults-dimuon_mass_8p7_11p1-dimuon_pt_10_50_costhbins.json'
    #costhbinnedfit_results = inpath+'/fitresults-dimuon_mass_8p7_11p1-dimuon_pt_10_50_costhbinned_fit.root'
    #chi2_outfile = outpath+'/data_fastmc_comparison'
    year=int(args.year)
    frame='HX'
    with open(jsonbinning, 'r') as jsonf:
        frame = json.load(jsonf)['frame']
    lumi = '34.88' if year == 2016 else "37.15"

    maxphotoneta=1.5
    f=open(args.configfile,'r')
    cuts=json.load(f).get('cut_variables')
    print cuts
    if cuts:
        for c in cuts: 
            if c['name']=='photon_eta':maxphotoneta=float(c['max'])
    else: print 'Could not read "cut_variables" from config file'
    print '|eta(photon)| < '+str(maxphotoneta)

    outfile_name = chi2_outfile+'fastmc_histograms_costh{}.root'.format(frame)

    folder='20180901_WITHselection'
    chib1_pattern = '/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/toymc/{}/toymc_chib-R_{}_WITHselection.root'
    chib2_pattern = '/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/toymc/{}/toymc_chib-R1_{}_R2_{}_WITHselection.root'

    infiles = [chib1_pattern.format(folder, R) for R in (
        '2o3',
        '0',
        '1'
       )]
    infiles += [chib2_pattern.format(folder, R1, R2) for R1,R2 in (
        ('2o5','2o5'),
        ('0','0'),
        ('1','0'),
        ('0','1')
        )]

    # Most of the cuts are already done at generation
    costhvar='costh_{}'.format(frame)
    cols = [ costhvar, 
            'Jpsi{Pt,Rap}',
            '{photon,mu{N,P}}{Pt,Eta}']
    selections = (lambda d: sf.photon_sel(d, sf.flat_pt(0.4, maxphotoneta)))
    if year == 2017: selections=(
        lambda d: sf.photon_sel(d, sf.flat_pt(0.4, maxphotoneta)),
        lambda d: sf.jpsi_kin_sel(d, 10, 50, 1.2),
        lambda d: sf.single_muon_sel(d, sf.flat_pt(4.5,1.4))
        )

    print "Selections:"
    for sel in selections: print inspect.getsource(sel)
    print 'maxphotoneta={}'.format(maxphotoneta)


    chib1_regex = "R_(\d+o?\d*)"
    chib2_regex = "R1_(\d+o?\d*)_R2_(\d+o?\d*)"
    def get_chistate_dict(filename):
        R1='-99'
        R2='-99'
        state='-99'
        name='ERROR'
        m=re.search(chib1_regex, filename)
        if (m and (len(m.groups())>0)):
            R1 = m.group(1)
            name = "1_R_{}".format(R1)
        m=re.search(chib2_regex, filename)
        if(m and (len(m.groups())>1)):
                R1=m.group(1)
                R2=m.group(2)
                name = "2_R1_{}_R2_{}".format(R1,R2)
        return {"file":filename, 
                "state":state,
                "histoname":"chib{}".format(name)}

    binning=None
    with open(jsonbinning, 'r') as jsonf:
        tmpbins = json.load(jsonf)['costh_bins']
        binning = [b[0] for b in tmpbins]
        binning.append(tmpbins[-1][-1])
        binning = np.array(binning).astype(np.float)
    print 'Binning:', binning


    file_infos=[get_chistate_dict(f) for f in infiles]

    print 'Generating histograms'
    outfile = r.TFile(outfile_name,"recreate")
    for i in file_infos:
        df = get_dataframe(i['file'], treename='toy_mc', columns=cols)
        sel_df = apply_selections(df, selections)
        h_sel_costh = create_histogram(sel_df[costhvar].abs(), (len(binning)-1,binning))
        outfile.cd()
        h_sel_costh.Write(i['histoname'],r.TObject.kWriteDelete)
        print '    Written {0}'.format(i['histoname'])
    outfile.Close()

    import calc_chi2_prob
    calc_chi2_prob_args = argparse.Namespace(
        fitres=costhbinnedfit_results, 
        histfile=outfile_name, 
        bininfo=jsonbinning,
        saveas=chi2_outfile,
        lumi=lumi)
    print 'Calculating chi2 probability...'
    calc_chi2_prob.main(calc_chi2_prob_args)
    