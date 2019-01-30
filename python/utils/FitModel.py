#!/usr/bin/env python

import ROOT as r
import ROOT.RooFit as rf

import seaborn as sb
import matplotlib.pyplot as plt

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.roofit_utils import (
    ws_import, get_var, get_chi2_ndf, get_corr_matrix
)
from utils.misc_helpers import create_random_str
from utils.plot_helpers import setup_latex, put_on_latex


def set_var(wsp, varname, val):
    """
    Set the passed value to the variable in the workspace. If the variable is
    not already present, create it first.
    """
    var = get_var(wsp, varname)
    if var is None:
        logging.debug('Variable \'{}\' not already present in workspace'
                      .format(varname))
        wsp.factory("{}[{}]".format(varname, val))
    else:
        var.setVal(val)


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
        self.legpos = None
        self.nevent_yields = None
        self.floating_costh = []
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


    def fit(self, wsp, savename, add_cut='', weighted_fit=False):
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
            weighted_fit (boolean, optional): Use the weight variable (that has)
                to be already present in the passed RooDataSet in the workspace
                as weight in the fit.
        """
        fit_data = wsp.data('full_data')
        if add_cut:
            fit_data = fit_data.reduce(add_cut)

        fit_args = (
            rf.Minos(True),
            rf.NumCPU(4),
            rf.Save(True),
            rf.Extended(True),
            rf.Offset(False)
        )

        if weighted_fit:
            # Choosing the SumW2Error(True) option here should give errors on
            # the fitted parameters that reflect the available statistics of the
            # sample (at least according to the RooFit documentation)
            fit_args += (rf.SumW2Error(True),)

        fit_results = wsp.pdf(self.full_model).fitTo(fit_data, *fit_args)

        fit_results.Print()
        logging.info('Fit status = {}, covQual = {}'
                     .format(fit_results.status(), fit_results.covQual()))

        set_var(wsp, '__fit_status__', fit_results.status())
        set_var(wsp, '__cov_qual__', fit_results.covQual())

        wsp.saveSnapshot('snap_{}'.format(savename), wsp.allVars())
        fit_results.SetName('fit_res_{}'.format(savename))
        ws_import(wsp, fit_results)


    def plot(self, wsp, pdfname, snapname='', add_cut='', **kwargs):
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

        Keyword Args:
            logy (bool): Set log on y scale of distribution plot
            weighted_fit (bool): Assume that the fit has done using weights and
                calculate the data uncertainties using these weights.
        """
        fitresname = None
        if snapname:
            wsp.loadSnapshot(snapname)
            fitresname = snapname.replace('snap', 'fit_res')

        # Starting from ROOT v6.12 it is possible to set this on individual axis
        r.TGaxis.SetMaxDigits(3)

        mvar = get_var(wsp, self.mname)
        plot_data = wsp.data('full_data')
        if add_cut:
            plot_data = plot_data.reduce(add_cut)
        frame = mvar.frame(rf.Bins(80))
        frame.SetTitle("")
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.GetYaxis().SetTitle('Events / {:.1f} MeV'
                                  .format((mvar.getMax() - mvar.getMin()) / 80 * 1000))

        full_pdf = wsp.pdf(self.full_model)

        leg = self._setup_legend()

        if kwargs.pop('weighted_fit', True):
            plot_data.plotOn(frame, rf.MarkerSize(0.8), rf.Name('data_hist'),
                             rf.DataError(r.RooAbsData.SumW2))
        else:
            plot_data.plotOn(frame, rf.MarkerSize(0.8), rf.Name('data_hist'))

        full_pdf.plotOn(frame, rf.LineWidth(2),
                        # rf.ProjWData(plot_data),
                        rf.Name('full_pdf_curve'))

        leg.AddEntry(frame.getCurve('full_pdf_curve'), 'sum', 'l')

        for name, lst, lcol, legentry in self.components:
            full_pdf.plotOn(frame, rf.Components(name), rf.LineStyle(lst),
                            rf.LineColor(lcol), rf.LineWidth(2), rf.Name(name))
            leg.AddEntry(frame.getCurve(name), legentry, 'l')


        info_text = []
        if fitresname is not None:
            fit_res = wsp.genobj(fitresname)
            chi2, ndf = get_chi2_ndf(fit_res, frame, 'full_pdf_curve',
                                     'data_hist')
            info_text.append((0.15, 0.875,
                              '#chi^{{2}}/ndf = {:.1f} / {}'.format(chi2, ndf)))

        pull_frame = mvar.frame(rf.Bins(80))
        hpull = frame.pullHist('data_hist', 'full_pdf_curve', True)
        hpull.SetMarkerSize(0.8)
        pull_frame.addPlotable(hpull, 'P')

        pull_frame.SetTitle("")
        pull_frame.GetYaxis().SetTitle("pull")
        pull_frame.GetXaxis().SetTitleSize(0.08)
        pull_frame.GetYaxis().SetTitleSize(0.08)
        pull_frame.GetXaxis().SetLabelSize(0.08)
        pull_frame.GetYaxis().SetLabelSize(0.08)
        pull_frame.GetYaxis().SetTitleOffset(0.4)
        pull_frame.GetYaxis().SetRangeUser(-5.99, 5.99)

        latex = setup_latex()
        can = r.TCanvas(create_random_str(16), '', 50, 50, 600, 600)
        can.cd()

        pad = r.TPad('mass_pad', 'mass_pad', 0, 0.3, 1, 1)
        r.SetOwnership(pad, False)
        pad.Draw()
        pad.cd()
        if kwargs.pop('logy', False):
            pad.SetLogy()
            pdfname = pdfname.replace('mass_fit', 'mass_fit_log')
        frame.Draw()
        leg.Draw()
        put_on_latex(latex, info_text, ndc=True)


        can.cd()
        pull_pad = r.TPad('pull_pad', 'pull_pad', 0, 0, 1, 0.3)
        r.SetOwnership(pull_pad, False)
        pull_pad.Draw()
        pull_pad.SetGridy()
        pull_pad.SetBottomMargin(0.2)
        pull_pad.cd()
        pull_frame.Draw()

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


    def print_fit_params(self, wsp, fitresname):
        """
        Print all free fit parameters

        Args:
            wsp (ROOT.RooWorkspace): workspace where the model is stored
            fitresname (str): Name of fit results that will be printed
        """
        logging.debug('printing free parameters of %s', fitresname)
        fit_res = wsp.genobj(fitresname)
        free_params = fit_res.floatParsFinal()
        for i in xrange(free_params.getSize()):
            var = free_params[i]
            print('{} = {:.4f} +/- {:.4f}'.format(var.GetName(), var.getVal(),
                                                  var.getError()))


    def plot_corr_matrix(self, wsp, fitresname, pdfname):
        """
        Print the correlation matrix of a given fitresult
        """
        fit_res = wsp.genobj(fitresname)
        corr_matrix = get_corr_matrix(fit_res)
        fig = plt.figure(figsize=corr_matrix.shape)

        sb.heatmap(corr_matrix, annot=True, fmt='.2f', square=True)
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)

        fig.savefig(pdfname, bbox_inches='tight')


    def fix_params(self, wsp, param_vals):
        """
        Fix the parameters in the workspace to given or current values

        Args:
            wsp (ROOT.RooWorkspace): workspace containing all the variables
            param_vals (list of tuples): tuples where first element is the name
                of the parameter and the second is the value to which it should
                be fixed. If the second parameter is None, it will be fixed to
                the current value in the workspace
        """
        for par, val in param_vals:
            if val is None:
                val = get_var(wsp, par).getVal()
                debug_msg = 'Fixing {} to current value in workspace: {}'
            else:
                debug_msg = 'Fixing {} to {}'

            logging.info(debug_msg.format(par, val))
            get_var(wsp, par).setVal(val)
            get_var(wsp, par).setConstant(True)


    def release_params(self, wsp, param_names):
        """
        Release the parameters in the workspace

        Args:
            wsp (ROOT.RooWorkspace): workspace containing all the variables
            param_names (list of strings): Names of the parameters in the
                workspace to release
        """
        for par in param_names:
            get_var(wsp, par).setConstant(False)


    def _setup_legend(self):
        """
        Setup the legend and return it

        Returns:
            ROOT.TLegend: legend with basic settings
        """
        leg = r.TLegend(*self.legpos)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.SetBorderSize(0)

        return leg
