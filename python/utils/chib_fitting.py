#!/usr/bin/env python

from utils.FitModel import FitModel

from ROOT import gSystem
gSystem.Load('$CHIB_CHIC_POLFW_DIR/general/bin/CustomRooFitFunctions')

"""
Module containing the basic setup for the chib fitting
"""

import json

class ChibMassModel(FitModel):
    """
    Chib mass model
    """
    def __init__(self, config_file, outfolder='.'):
        """
        Args:
            config_file: Name of the json file containing the model definition
        """
        self.config_file = config_file
        from ROOT import RooMsgService
        RooMsgService.instance().setSilentMode(True)
        RooMsgService.instance().setGlobalKillBelow(4)

        with open(config_file, 'r') as f:
            data = json.load(f)
            self.mname = data["chi_fitvar"]["name"]
            self.full_model = data["chi_modelname"]
            self.model_strings = data["chi_model"]
            self.fitvarmin = data["chi_fitvar"]["min"]
            self.fitvarmax = data["chi_fitvar"]["max"]
            self.nevent_vars = data['sweight_yields']
        
        # Hard coded:
        self.background='background'
        self.components = (
            ('chib1', 7, 417, '#chi_{b1}'),
            ('chib2', 7, 632, '#chi_{b2}'),
            ('background', 7, 1, 'background')
            )
        self.legpos = (0.18, 0.5, 0.38, 0.7)        
        self.max_iterations=10
        self.outfolder=outfolder+'/'


    def define_model(self, ws):
        """
        Build the chib mass model in the passed workspace

        Args:
            ws (ROOT.RooWorkspace): Workspace into which the mass model is
                constructed
        """  

        for fac in self.model_strings:
            ws.factory(fac)
        ws.saveSnapshot('start_parameter', ws.allVars())

    def fit(self, wsp, savename, add_cut=''):
        import ROOT.RooFit as rf
        import ROOT as r
        import re
        from utils.roofit_utils import ws_import
        import logging
        logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')
        """
        Fit the defined model to the data present in the workspace (under name
        'full_data').

        Args:
            wsp (ROOT.RooWorkspace): workspace containing the model as well as
                the data to which it should be fitted
            savename(str): basename for storing the snapshot of the variables
                after the fit as well as the fit result pointer into the
                workspace. Snapshots prepended with 'snap_' and fit results
                prepended with 'fit_res' for later retrieval.
            add_cut(str, optional): Additional cut to apply to the dataset
                before fitting. NOTE: all used variables have to be present in
                the dataset.
        """
        wsp.loadSnapshot('start_parameter')
        fit_data = wsp.data('full_data').reduce(add_cut)
        x=wsp.var(self.mname)
        plot=x.frame()

        # Start variables from config file are for costh integrated,
        # therefore divide N1, N2 and Nbg by number of bins
        # Everything hard coded
        for v in ('N1','N2','Nbg'):
            n=wsp.var(v)
            n.setVal(n.getValV()/5.)
            print 'Startvalue of {} set to {}.'.format(v,n.getValV())

            
        # fit first the background
        wsp.var(self.mname).setRange("LOWSIDE", self.fitvarmin, 9.8)
        wsp.var(self.mname).setRange("HIGHSIDE", 9.98, self.fitvarmax)
        wsp.pdf(self.background).fitTo(fit_data,
                                       rf.Minos(True),
                                       rf.NumCPU(8),
                                       rf.Range("LOWSIDE,HIGHSIDE"),
                                       rf.Offset(True),
                                       rf.Minimizer('Minuit2', 'migrad')
                                       );  
       
        #then fix the background and fit again
        wsp.var('c0').setConstant(True)
        wsp.var('c1').setConstant(True)
        wsp.pdf(self.full_model).fitTo(fit_data,
                                        rf.Minos(True),
                                        rf.NumCPU(8),
                                        rf.Extended(True),
                                        rf.Offset(True),
                                        rf.Minimizer('Minuit2', 'migrad'))
        
        #then fix all but N and fit again
        for vn in ['c0','c1','ah','al','sigma1','mu1']:
            wsp.var(vn).setConstant(True)
        wsp.pdf(self.full_model).fitTo(fit_data,
                                        rf.Minos(True),
                                        rf.NumCPU(8),
                                        rf.Extended(True),
                                        rf.Offset(True),
                                        rf.Minimizer('Minuit2', 'migrad'))
        
        r.gROOT.SetBatch(True)
        c= r.TCanvas('c','c',1000,800)
        plot=plot.emptyClone('fixed_bg_plot')
        fit_data.plotOn(plot)
        wsp.pdf(self.full_model).plotOn(plot)
        plot.Draw()
        c.SaveAs(self.outfolder+'fixed_bg_plot{}.pdf'.format(savename))

        # then release all and start the fit
        for vn in ['c0','c1','ah','al','sigma1','mu1']:
            wsp.var(vn).setConstant(False)
        

# HARDCODED #       
        #if '0' in savename  and '2017' in self.outfolder:
            #wsp.var('N1').setVal(470)
            #wsp.var('N2').setVal(269)
            #wsp.var('Nbg').setVal(3880)
            #wsp.var('ah').setVal(2.24)
            #wsp.var('al').setVal(0.72)
            #wsp.var('sigma1').setVal(0.00437)
            #wsp.var('c0').setVal(0.103)
            #wsp.var('c1').setVal(-0.0284)
        #if '4' in savename and '2017' in self.outfolder:
            #wsp.var('N1').setVal(350)
            #wsp.var('N2').setVal(200)
            ##wsp.var('ah').setVal(2)
            ##wsp.var('al').setVal(0.78)
            #wsp.var('sigma1').setVal(0.0045)
# ENDE HARDCODED # 


        for i in xrange(self.max_iterations):

            print 'ITERATION {}, {}'.format(i,savename)
            vars=wsp.allVars()
            it=vars.createIterator()    
            for j in xrange(vars.getSize()):
                it.Next().Print()
            fit_results = wsp.pdf(self.full_model).fitTo(fit_data,
                                                         rf.Minos(True),
                                                         rf.NumCPU(8),
                                                         rf.Save(True),
                                                         rf.Extended(True),
                                                         rf.Offset(True),
                                                         rf.Minimizer('Minuit2', 'minimize'))

            fit_results.Print()
            logging.info('Fit status = {}, covQual = {}'.format(fit_results.status(), fit_results.covQual()))
            r.gROOT.SetBatch(True)
            c= r.TCanvas('c','c',1000,800)
            plot=plot.emptyClone('plot{}'.format(i))
            fit_data.plotOn(plot)
            wsp.pdf(self.full_model).plotOn(plot)
            plot.Draw()
            c.SaveAs(self.outfolder+'{}_iteration{}.pdf'.format(savename,i))

            if fit_results.status()==0 and i!=0:
                print "CONVERGED in {} iterations, {}".format(i,savename)
                break

            if i==self.max_iterations-1:
                print "MAXIMUM NUMBER OF ITERATIONS ({}) REACHED, {}".format(self.max_iterations,savename)


        wsp.saveSnapshot('snap_{}'.format(savename), wsp.allVars())
        fit_results.SetName('fit_res_{}'.format(savename))
        ws_import(wsp, fit_results)
        if add_cut:
            fit_data.SetName('data_{}'.format(savename))
            ws_import(wsp, fit_data)

