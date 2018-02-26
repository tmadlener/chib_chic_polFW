#!/usr/bin/env python

import ROOT as r
import ROOT.RooFit as rf

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.roofit_utils import ws_import, get_var
from utils.misc_helpers import create_random_str

class FitModel(object):
    """
    Base class for different models for fitting
    """
    def __init__(self, mname):
        """
        Args:
            mname (str): name of the mass variable (as it is stored in the
                TTree and subsequently in the workspace)
        """
        self.full_model = None
        self.components = None
        self.mname = None
        raise NotImplementedError('__init__ has to be defined by derived class')

    def define_model(self, wsp):
        """
        Define the fit model in the workspace by importing everything into it.

        Args:
            wsp (ROOT.RooWorkspace): Workspace into which the fit model is
                imported
        """
        raise NotImplementedError('_define_model has to be defined by derived '
                                  'class')

    def fit(self, wsp, savename, add_cut=''):
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
        fit_data = wsp.data('full_data').reduce(add_cut)

        fit_results = wsp.pdf(self.full_model).fitTo(fit_data,
                                                     rf.Minos(True),
                                                     rf.NumCPU(4),
                                                     rf.Save(True),
                                                     rf.Extended(True),
                                                     rf.Offset(False))

        fit_results.Print()
        logging.info('Fit status = {}, covQual = {}'
                     .format(fit_results.status(), fit_results.covQual()))

        wsp.saveSnapshot('snap_{}'.format(savename), wsp.allVars())
        fit_results.SetName('fit_res_{}'.format(savename))
        ws_import(wsp, fit_results)


    def plot(self, wsp, pdfname, snapname='', add_cut=''):
        """
        Make a plot of the model and the fitted data and save as pdf.

        Args:
            wsp (ROOT.RooWorkspace): workspace where all the information
                (including data and model) is stored
            pdfname (str): Name of the created pdf under which the plot will be
                stored.
            snapname (str, optional): Name of snapshot that will be loaded
                before plotting
            add_cut (str, optional): Additional cut to apply to the dataset
                before plotting. NOTE: all used variables have to be present in
                the dataset
        """
        if snapname:
            wsp.loadSnapshot(snapname)

        plot_data = wsp.data('full_data').reduce(add_cut)
        frame = get_var(wsp, self.mname).frame(rf.Bins(80))
        full_pdf = wsp.pdf(self.full_model)

        plot_data.plotOn(frame, rf.MarkerSize(0.8))
        full_pdf.plotOn(frame, rf.LineWidth(2))

        for name, lst, lcol in self.components:
            full_pdf.plotOn(frame, rf.Components(name), rf.LineStyle(lst),
                            rf.LineColor(lcol), rf.LineWidth(2))

        can = r.TCanvas(create_random_str(16), '', 50, 50, 600, 600)
        can.cd()
        frame.Draw()

        can.SaveAs(pdfname)


    def plot_fit_params(self, wsp, pdfname, snapname=''):
        """
        Plot all free fit parameters onto canvas and save as a pdf.

        Args:
            wsp (ROOT.RooWorkspace): workspace where the model is stored
            pdfname (str): Name of the created pdf under which the plot will be
                stored
            snapname (str, optional): Name of snapshot that will be loaded
                before plotting
        """
        if snapname:
            wsp.loadSnapshot(snapname)
        frame = get_var(wsp, self.mname).frame(rf.Title('Fit Results'))
        full_pdf = wsp.pdf(self.full_model)

        full_pdf.paramOn(frame, rf.Layout(0.1, 0.9, 0.9),
                         rf.Format('NEU', rf.AutoPrecision(2)))

        can = r.TCanvas(create_random_str(32), 'rcan', 600, 600)
        can.cd()
        frame.findObject('{}_paramBox'.format(full_pdf.GetName())).Draw()
        can.SaveAs(pdfname)
