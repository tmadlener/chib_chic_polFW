
infilename = '/afs/hephy.at/work/j/jnecker/data/chib_results/bin_8_50/fitresults-dimuon_mass_8p7_11p1-dimuon_pt_8_50.root'
workspacename = 'ws-dimuon_mass_8p7_11p1-dimuon_pt_8_50_chib1P1S'
varname = 'chi_mass_rf1S'
pdfs2plot = [{'name':'model', 'legend':'Chi Mass Model', 'color':1}] #first is model

from ROOT import TFile, TCanvas
from utils.chib_fitting import ChibMassModel # Needed for user defined functions

if __name__=='__main__':
    infile = TFile.Open(infilename)
    ws = infile.Get(workspacename)
    var = ws.var(varname)
    model = ws.pdf(pdfs2plot[0]['name'])
    frame = var.frame()
    c = TCanvas('c','c',1000, 800)
    model.plotOn(frame)
    frame.Draw()
    raw_input('...')