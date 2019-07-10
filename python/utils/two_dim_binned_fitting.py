#!/usr/bin/env python
"""
TODO
"""

import re

import ROOT as r
import ROOT.RooFit as rf
import numpy as np

from collections import OrderedDict

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.roofit_utils import (
    get_var, try_factory, set_var, ws_import, fix_params, get_chi2_ndf,
    get_var_err, calc_info_crit
)
from utils.misc_helpers import (
    parse_binning, get_bin_cut_root, combine_cuts, create_random_str,
    replace_all, parse_func_var
)
from utils.plot_helpers import (
    setup_legend, _setup_canvas, setup_latex, put_on_latex
)
from utils.graph_utils import assign_x
from utils.plot_decoration import YLABELS

import utils.RooDoubleCB
import utils.RooErfExponential
import utils.RooPowerlawExponential

BIN_WIDTH = 0.00625 # GeV

def get_bins(binning1, binning2, binvar1, binvar2):
    """
    Get a dict with bin borders and identifiers
    """
    bins1 = zip(binning1[:-1], binning1[1:])
    bins2 = zip(binning2[:-1], binning2[1:])

    bins = OrderedDict()

    for ibin1, bin1 in enumerate(bins1):
        for ibin2, bin2 in enumerate(bins2):
            name = '_'.join([binvar1, str(ibin1), binvar2, str(ibin2)])
            bins[name] = (bin1, bin2)

    return bins


def get_bin_cut(bin_vars, bin_borders):
    """
    Get the bin cut
    """
    cuts = []
    for i, bv in enumerate(bin_vars):
        cuts.append(get_bin_cut_root(bv, *bin_borders[i]))

    return combine_cuts(cuts)


def _full_model_expr(model_type, name, sub_models, event_yields):
    """
    Construct the full model expression

    NOTE: Might fail at runtime when this string is used as factory expression
    """
    sub_expr = ', '.join([
        ' * '.join(['{}'.format(y), '{}'.format(sub_models[i])])
        for i, y in enumerate(event_yields)
    ])

    return '{}::{}({})'.format(model_type, name, sub_expr)


def get_bin_info_text(bin_name, bins):
    """Get text that can be put into a TLatex specifying which bin is plotted"""
    var_rgx = re.compile(r'(.*)_[0-9]+_(.*)_[0-9]+')
    match = var_rgx.match(bin_name)

    var1 = YLABELS.get(match.group(1), match.group(1))
    var2 = YLABELS.get(match.group(2), match.group(2))

    def fmt_val(value):
        """Format a value somewhat nicely"""
        return '{:.3f}'.format(value).rstrip('0').rstrip('.')

    str1 = '{} < {} < {}'.format(fmt_val(bins[bin_name][0][0]), var1, fmt_val(bins[bin_name][0][1]))
    str2 = '{} < {} < {}'.format(fmt_val(bins[bin_name][1][0]), var2, fmt_val(bins[bin_name][1][1]))

    if var1 == 'JpsiPt':
        str1 += ' GeV'
    if var2 == 'JpsiPt':
        str2 += ' GeV'

    return [
        (0.22, 0.575, str1),
        (0.22, 0.5, str2)
    ]

    print bin_name
    return []


class BinnedFitModel(object):
    """
    Fit model for doing a simultaneous fit in multiple bins
    """
    def __init__(self, config):
        """
        Read all the info from the config and do some computation which are then
        needed throughout
        """
        # TODO: read the configuration from a json file and make that the
        # argument to the function
        self.binning = np.array([
            parse_binning(config['binning'][0]),
            parse_binning(config['binning'][1])
        ])
        # variables as they are defined and used for deciding the binning
        # inside ROOT with string cuts
        self.bin_cut_vars = config['bin_vars']

        # For the internal use and especially for naming the variables it is not
        # possible to have functional variable definitions. Simply use whatever
        # is the argument of func(var) (i.e. var)
        self.bin_vars = [
            parse_func_var(v)[0] for v in config['bin_vars']
        ]

        self.bins = get_bins(self.binning[0], self.binning[1],
                             self.bin_vars[0], self.bin_vars[1])

        self.fit_var = config['fit_variable']
        self.fit_range = config['fit_range']
        self.expression_strings = config['expression_strings']
        self.proto_params = config['proto_parameters']

        self.sub_models = config['sub_models']
        self.nevent_yields = [m['event_yield'] for m in self.sub_models]

        self.components = [(m['name'], m['plot']) for m in self.sub_models]

        self.full_model = config['full_model']['name']

        # prototype full model expression, where the bin model can be achieved
        # with simple replacements
        self.proto_fexpr = _full_model_expr(config['full_model']['type'],
                                            self.full_model,
                                            [c[0] for c in self.components],
                                            self.nevent_yields)

        self.fix_vars = []
        if 'fix_vars' in config:
            for var in config['fix_vars']:
                self.fix_vars.append((var.keys()[0], var.values()[0]))

        self.plot_config = config['plot_config']

        # self.config = config # for development / debugging


    def get_load_vars(self):
        """
        Get the variable names and bounds that are necessary to define the
        variables in a workspace such that data can be imported from a TTree.

        Returns:
            load_vars, bounds: List of tuples of strings, where first element is
                the name of the variable and the second is the boundaries as
                comma separated list
        """
        variables = [
            (self.fit_var, ",".join([str(v) for v in self.fit_range])),
        ]
        for ivar, var in enumerate(self.bin_cut_vars):
            var_name, var_func = parse_func_var(var)
            bounds_up = np.max(self.binning[ivar])

            # For now only handle the 'abs(var)' case separately
            if var_func is not None and 'abs' in var_func.__name__:
                bounds_low = -bounds_up
            else:
                bounds_low = np.min(self.binning[ivar])

            variables.append((var_name, "{}, {}".format(bounds_low, bounds_up)))

        return variables


    def define_model(self, wsp):
        """
        Build the model into the passed workspace

        Args:
            wsp (ROOT.RooWorkspace): Workspace into which the model should be
                constructed. The data has to be already present!

        Returns:
           success (boolean): True if all of the expressions were successfully
               handled by the workspace factory, False if at least one of them
               failed
        """
        success = []
        for expr in self.expression_strings:
            success.append(try_factory(wsp, expr))

        for bin_name in self.bins:
            for param, var_str in self.proto_params:
                # Check if we have a functional dependency or not
                if isinstance(var_str, (tuple, list)):
                    expr = self._get_func_expr(wsp, bin_name, param, var_str)
                    for dep_param, dep_var_str in self.proto_params:
                        if dep_param is not param:
                            expr = expr.replace(dep_param, dep_param + '_' + bin_name)
                # For a non-functional dependency we want to have a prototype
                # variable in the workspace
                else:
                    expr = self._get_par_expr(wsp, bin_name, param, var_str)

                success.append(try_factory(wsp, expr))

            for model in self.sub_models:
                name = '{}_{}'.format(model['name'], bin_name)
                expr = model['expression'].format(name, self.fit_var)

                for param, var_str in self.proto_params:
                    expr = expr.replace(param, param + '_' + bin_name)

                success.append(try_factory(wsp, expr))

            full_expr = self.proto_fexpr
            for model in self.sub_models:
                full_expr = full_expr.replace(model['name'],
                                              model['name'] + '_' + bin_name)
            for param, var_str in self.proto_params:
                full_expr = full_expr.replace(param, param + '_' + bin_name)

            full_expr = full_expr.replace(self.full_model,
                                          self.full_model + '_' + bin_name)
            success.append(try_factory(wsp, full_expr))

        fix_params(wsp, self.fix_vars)

        return all(success)


    def fit(self, wsp, verbosity=0):
        """
        Import the passed data into the workspace, construct the combined
        negative log-likelihood (NLL) and run the simultaneous fit.

        Args:
            wsp (ROOT.RooWorkspace): The workspace that holds the model and into
                which the data will be imported
            data (ROOT.TTree): The (unbinned) data to which the model should be
                fitted. The binning will be done in this function
            verbosity (int): Either -1, 0, 1 to steer how much output is
                generated by the fit. Follows the TMinuit convention:
                -1 -> quite, 0 -> normal, 1 -> verbose
        """
        nll_args = ( # TODO: define meaning full args
            rf.NumCPU(2),
            rf.Extended(True),
            rf.Offset(True)
        )
        sim_nll = self._create_nll(wsp, nll_args)

        minimizer = self._setup_minimizer(sim_nll, verbosity)

        fit_status = 0
        cov_qual = 0

        for i in range(3):
            logging.info('starting migrad')
            minimizer.migrad()
            logging.info('starting minos')
            minimizer.minos()

            fit_results = minimizer.save()

            fit_results.Print()
            logging.info('Fit status = {}, covQual = {}'
                         .format(fit_results.status(), fit_results.covQual()))

            fit_status *= 10
            fit_status += fit_results.status()
            cov_qual *= 10
            cov_qual += fit_results.covQual()


            if fit_results.status() == 0 and fit_results.covQual() == 3:
                break

        set_var(wsp, '__fit_status__', fit_status, create=True)
        set_var(wsp, '__cov_qual__', cov_qual, create=True)

        wsp.saveSnapshot('snap_two_dim', wsp.allVars())
        fit_results.SetName('fit_res_two_dim')
        ws_import(wsp, fit_results)
        sim_nll.SetName('sim_nll')
        ws_import(wsp, sim_nll)


    def plot(self, wsp, verbose=False, publication=False, preliminary=False):
        """
        Plot the fit results from the passed workspace

        Args:
            wsp (ROOT.RooWorkspace): The workspace that holds the model and the
                data as well as the fit results
            verbose (boolean): Print some more information about the fit status
                onto the plots (mainly for debugging)
            publication (boolean, optional): Make the plots in publication
                quality (negates the verbose argument automatically)
            preliminary (boolean, optional): Add the 'preliminary' label
        """
        data = wsp.data('full_data')

        cans = OrderedDict()
        g_chi2, g_ndf = 0, 0

        fvar = get_var(wsp, self.fit_var)
        n_bins = int(round(((fvar.getMax() - fvar.getMin()) / BIN_WIDTH)))

        wsp.loadSnapshot('snap_two_dim')

        for bin_name, bin_borders in self.bins.iteritems():
            frame, leg = self._plot_bin(wsp, data, bin_name, bin_borders, n_bins,
                                        publication)
            b_chi2, b_ndf = get_chi2_ndf(frame, 'full_pdf_curve', 'data_hist')
            logging.debug('bin: {}: chi2 / number bins = {:.2f} / {}'
                          .format(bin_name, b_chi2, b_ndf))
            g_chi2 += b_chi2
            g_ndf += b_ndf

            can = _setup_canvas(None) # returns a TCanvasWrapper
            can.cd()

            pad = r.TPad('result_pad', 'result pad', 0, 0.3, 1, 0.98)
            r.SetOwnership(pad, False)

            if publication:
                pad.SetBottomMargin(0)

            pad.Draw()
            pad.cd()
            frame.Draw()
            can.add_tobject(pad)
            leg.Draw()
            can.add_tobject(leg)

            # pulls
            pull_frame = self._pull_plot(wsp, frame, publication)

            can.cd()
            pull_pad = r.TPad('pull_pad', 'pull pad', 0, 0, 1, 0.3)
            r.SetOwnership(pull_pad, False)
            if publication:
                pull_pad.SetTopMargin(0)
                pull_pad.SetBottomMargin(0.275)
            else:
                pull_pad.SetBottomMargin(0.2)
            pull_pad.Draw()
            pull_pad.SetGridy()
            pull_pad.cd()
            pull_frame.Draw()
            can.add_tobject(pull_pad)

            cans[bin_name] = can

        fit_res = wsp.genobj('fit_res_two_dim')
        n_float_pars = fit_res.floatParsFinal().getSize()
        g_ndf -= n_float_pars
        logging.debug('Floating parameters in fit: {}'.format(n_float_pars))
        logging.info('Global chi2 / ndf = {:.2f} / {}'.format(g_chi2, g_ndf))

        add_info = []

        # loop again over all canvases and put the global chi2 / ndf there
        add_info.append((0.15, 0.875,
                         '(global) #chi^{{2}} / ndf = {:.1f} / {}'.format(g_chi2, g_ndf)))

        if verbose and not publication:
            status = get_var(wsp, '__fit_status__').getVal()
            cov_qual = get_var(wsp, '__cov_qual__').getVal()
            add_info.append((0.15, 0.825,
                             'status = {}, covQual = {}'.format(int(status),
                                                                int(cov_qual))))
            # Check if the simultaneous negative log-likelihood is present in
            # the workspace or create it if it is not yet present.
            # Since the snapshot values of the fit result are already loaded it
            # can be directly evaluated to get the value at the minimum
            sim_nll = get_var(wsp, 'sim_nll')
            if sim_nll is None:
                sim_nll = self._create_nll(wsp, (rf.Extended(True), rf.Offset(True)))
            min_nll = sim_nll.getVal()

            aic = calc_info_crit(fit_res, min_nll)
            bic = calc_info_crit(fit_res, min_nll, data.numEntries())
            add_info.append((0.15, 0.775,
                             'AIC = {:.1f}, BIC = {:.1f}'.format(aic, bic)))


        latex = setup_latex()
        if publication:
            latex.SetTextSize(0.04)
            add_info[0] = (0.185, 0.89, add_info[0][2])


        for bin_name, can in cans.iteritems():
            res_pad = can.attached_tobjects[0]
            res_pad.cd()
            put_on_latex(latex, add_info, ndc=True)
            bin_info = get_bin_info_text(bin_name, self.bins)
            put_on_latex(latex, bin_info, ndc=True)

        # for easier debugging at the moment
        return cans


    def plot_fit_params(self, wsp):
        """
        Plot all free fit parameters onto canvas
        Args:
        wsp (ROOT.RooWorkspace): workspace where the model is stored
        """
        fit_var = get_var(wsp, self.fit_var)

        cans = OrderedDict()

        for bin_name in self.bins:
            frame = fit_var.frame(rf.Title('Fit Results'))

            plotname = '{}_{}'.format(self.full_model, bin_name)
            full_pdf = wsp.pdf(plotname)

            full_pdf.paramOn(frame, rf.Layout(0.1, 0.9, 0.9), rf.Format('NEU', rf.AutoPrecision(2)))

            can = r.TCanvas(create_random_str(32), 'rcan', 600, 600)
            can.cd()
            frame.findObject('{}_paramBox'.format(full_pdf.GetName())).Draw()

            cans[bin_name] = can
        #can.SaveAs(pdfname)

        # returns fit params for each fit

        return cans


    def plot_simvar_graphs(self, wsp, simvars, sym_uncer=False):
        """
        Create TGraphs for all the "simple" proto parameters

        Args:
            wsp (ROOT.RooWorkspace): Workspace holding the fit results from a
                simultaneous fit.
            simvars (list of strings): List of proto parameter names that do not
                have a functional dependency but instead are independent for
                each bin.
            sym_uncer (boolean, optional): Return the graphs with symmetric
                uncertainties instead of asymmetric ones (False)

        Returns:
            list of TGraphAsymmErrors: Graphs of the passed proto_parameters as
                a function of the binning variables. For each bin in either
                direction a graph is created as a function of the other
                direction. The name of the returned graph indicates on which
                variable the graph depends and in which bin of the other
                variable the values are taken from.
                E.g.: 'p_v_var1_bin_0' means that parameter 'p' is plotted as a
                function of 'var1' in bin 0 of 'var2'
        """
        #auxiliary definitions
        bin_var_X = self.bin_vars[0]
        bin_var_Y = self.bin_vars[1]
        #auxiliary dict for getting the parameters in the right order of bins
        p_name = OrderedDict()
        p_name[bin_var_X] = ('j','i')
        p_name[bin_var_Y] = ('i','j')

        #all bin means calculated previously
        bin_means = OrderedDict()
        for bin_var in self.bin_vars:
            for bin_name in self.bins:
                mean_name = '__'.join([bin_var, bin_name])
                bin_means[mean_name] = self.bin_mean(wsp, bin_var, bin_name)

        graphs = []

        for bintovar, bin_var in enumerate(self.bin_vars):
            #get parameter in the right binning order
            p_getname = '{param}_{x_var}_{'+p_name[bin_var][0]+'}_{y_var}_{'+p_name[bin_var][1]+'}'
            b_getname = p_getname[8:]

            for param in simvars:
                vals = []

                #get the values of the parameter for each bin of X,Y
                for i in xrange(len(self.binning[1-bintovar]) - 1):
                    vals.append([get_var(wsp, p_getname.format(param=param, x_var=bin_var_X, y_var=bin_var_Y, i=i, j=j)) for j in xrange(len(self.binning[bintovar]) - 1)])

                #plot graphs as a function of binning var
                graph = self.plot_free_pars(wsp, bin_var, vals, sym_uncer)

                #setting the correct mean (errors adjust accordingly), setting graph name
                for i in xrange(len(self.binning[1-bintovar]) - 1):
                    g_mean = np.array([bin_means['__'.join([bin_var, b_getname.format(x_var=bin_var_X, y_var=bin_var_Y, i=i, j=j)])] for j in xrange(len(self.binning[bintovar]) - 1)])
                    graph[i] = assign_x(graph[i], g_mean)
                    # graph name of the form [y_axis]_v_[x_axis]_bin_[bin of other bin_var]
                    sym_add = '_sym_uncer' if sym_uncer else ''
                    graph[i].SetName('{}_v_{}_bin_{}{}'.format(param, bin_var, i, sym_add))
                    graphs.append(graph[i])

        return graphs


    def plot_comvar_funcs(self, wsp, comvars):
        """
        Plot the parameters that have a functional dependency
        """
        funcs = []

        for el in comvars:
        # func stores the passed TF1, as well as the binning variable on which the function depends
            func, dep_var = self._tf1_helper(wsp, comvars[el][0], comvars[el][1])
            func.SetName('{}_v_{}'.format(el, dep_var))
            funcs.append(func)

        return funcs


    def _tf1_helper(self, wsp, nameold, vrs):
        """
        helper function that takes all vars in a string and replaces them by their
        value in the fit

        NOTE: only works if there is only dependence on the
        mean value of ONE binning variable (otherwise it wouldn't be a tf1)
        """
        v_names = vrs.split(', ')

        listp = [('{}'.format(elm), '['+str(i)+']') for i, elm in enumerate(v_names)]
        name = replace_all(nameold, listp)

        mean_rgx = re.compile(r'<(\w+)>')

        #TODO: error message if more than one element in this set?
        mean_var = list(set(mean_rgx.findall(name)))
        if len(mean_var) != 1:
            # Only warning here. Normally the TF1 constructor should emit
            # another error message because the function doesn't compile below
            logging.warning('Found dependency on more than one variable')

        mean_var = mean_var[0]
        name = name.replace('<{}>'.format(mean_var), 'x')


        f1 = r.TF1(create_random_str(), name,
                   self.binning[self.bin_vars.index(mean_var)][0],
                   self.binning[self.bin_vars.index(mean_var)][-1])
        for i, elm in enumerate(v_names):
            f1.SetParameter(i, get_var(wsp, elm).getVal())
            f1.SetParError(i, get_var(wsp, elm).getError())
            f1.SetParName(i, elm)

        return f1, mean_var


    def bin_mean(self, wsp, var, bin_name):
        """
        Get the mean value of the passed variable in the passed bin
        """
        bin_data = wsp.data('full_data').reduce(
            get_bin_cut(self.bin_cut_vars, self.bins[bin_name])
            )
        # Check if the variable can be taken as it is or if it has to be treated
        # specially
        if var in self.bin_cut_vars:
            meanval = bin_data.mean(get_var(wsp, var))
        else:
            # Check to which bin cut variable it belongs
            form = [v for v in self.bin_cut_vars if var in v]
            if len(form) > 1:
                logging.warning('Cannot uniquely identify to which bin '
                                'definition \'{}\' belongs'.format(var))
            # Create a formula var and add it to the temporary data set to
            # easily obtain the mean value
            binvarf = r.RooFormulaVar('binvar', 'temporary bin variable',
                                      form[0], r.RooArgList(get_var(wsp, var)))
            bin_var = bin_data.addColumn(binvarf)
            meanval = bin_data.mean(bin_var)

        return meanval


    def plot_free_pars(self, wsp, bin_var, vals_var, sym_uncer=False):
        """
        function that does the plotting of a parameter as a function of either
        binning variable
        """
        graph = []

        binb = self.binning[self.bin_vars.index(bin_var)]

        dependent = isinstance(vals_var[0][0], r.RooFormulaVar)
        # Needed for dependent parameter uncertainties
        fit_res = wsp.genobj('fit_res_two_dim')

        for var in vals_var:
            if dependent or sym_uncer:
                vals = np.array(
                    [get_var_err(wsp, v.GetName(), fit_res) for v in var],
                )
                y_val = vals[:, 0].astype(float) # weird conversion error other-
                y_elo = vals[:, 1].astype(float) # wise in constructor below
                y_ehi = y_elo
            else:
                y_val = np.array([v.getVal() for v in var])
                # RooFit reports lower errors with negative sign, but for
                # TGraphAsymmErrors they need to be positive
                y_elo = np.array([-v.getErrorLo() for v in var])
                y_ehi = np.array([v.getErrorHi() for v in var])

            #bin means and errors left to calculate externally because it simplifies the process
            x_val = np.array([0.5 * np.sum(b) for b in zip(binb[:-1], binb[1:])])
            x_err = np.array([0.5 * np.diff(b)[0] for b in zip(binb[:-1], binb[1:])])

            graph.append(r.TGraphAsymmErrors(len(x_val), x_val, y_val, x_err, x_err, y_elo, y_ehi))

        return graph


    def _pull_plot(self, wsp, frame, publication):
        """
        Make the frame containing the pulls
        """
        pulls = frame.pullHist('data_hist', 'full_pdf_curve', True)
        pulls.SetMarkerSize(0.8)
        pull_frame = get_var(wsp, self.fit_var).frame()
        if publication:
            pull_frame.addPlotable(pulls, 'PEX0')
            pulls.SetMarkerSize(0.7)
            pull_frame.GetXaxis().SetTitle(YLABELS.get('chicMass'))
        else:
            pull_frame.addPlotable(pulls, 'P')

        pull_frame.SetTitle("")
        pull_frame.GetYaxis().SetTitle("pull")
        if publication:
            pull_frame.GetXaxis().SetTitleSize(0.12)
            pull_frame.GetYaxis().SetTitleSize(0.12)
            pull_frame.GetXaxis().SetLabelSize(0.12)
            pull_frame.GetYaxis().SetLabelSize(0.12)
            pull_frame.GetYaxis().SetTitleOffset(0.625)
            pull_frame.GetYaxis().SetRangeUser(-3.99, 3.99)
            pull_frame.GetYaxis().SetNdivisions(507)
        else:
            pull_frame.GetXaxis().SetTitleSize(0.08)
            pull_frame.GetYaxis().SetTitleSize(0.08)
            pull_frame.GetXaxis().SetLabelSize(0.08)
            pull_frame.GetYaxis().SetLabelSize(0.08)
            pull_frame.GetYaxis().SetTitleOffset(0.4)
            pull_frame.GetYaxis().SetRangeUser(-5.99, 5.99)

        return pull_frame


    def _plot_bin(self, wsp, full_data, bin_name, bin_borders, n_bins,
                  publication=False):
        """Make the distribution plot for a given bin"""
        data_args = (rf.MarkerSize(0.7), rf.Name('data_hist'))
        fit_var = get_var(wsp, self.fit_var)
        frame = fit_var.frame(rf.Bins(n_bins))

        leg = setup_legend(*self.plot_config['legpos'])
        if publication:
            leg.SetTextSize(0.04)

        cut = get_bin_cut(self.bin_cut_vars, bin_borders)
        full_data.reduce(cut).plotOn(frame, *data_args)
        full_pdf = wsp.pdf(self.full_model + '_' + bin_name)
        full_pdf.plotOn(frame, rf.LineWidth(2), rf.Name('full_pdf_curve'))

        leg.AddEntry(frame.getCurve('full_pdf_curve'), 'Fit result', 'l')

        for name, settings in self.components:
            full_pdf.plotOn(frame, rf.Components(name + '_' + bin_name),
                            rf.LineStyle(settings['line']),
                            rf.LineColor(settings['color']),
                            rf.LineWidth(2),
                            rf.Name(name))
            leg.AddEntry(frame.getCurve(name), settings['label'], 'l')


        # calculate the bin width in MeV
        bw_mev = lambda v, n: (v.getMax() - v.getMin()) / n * 1000
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.GetYaxis().SetTitle('Events / {:.2f} MeV'
                                  .format(bw_mev(fit_var, n_bins)))

        frame.SetTitle("")

        return frame, leg


    def _create_nll(self, wsp, nll_args):
        """
        Create the combined negative log-likelihood (NLL) that is minimized in
        the fit
        """
        logging.info('Creating NLL for fit')
        sim_nll_list = r.RooArgList()

        data = wsp.data('full_data')

        for bin_name, bin_borders in self.bins.iteritems():
            logging.debug('Creating sub NLL for bin {}'.format(bin_name))
            cut = get_bin_cut(self.bin_cut_vars, bin_borders)
            bin_data = data.reduce(cut)

            bin_model = wsp.pdf(self.full_model + '_' + bin_name)
            sim_nll_list.add(bin_model.createNLL(bin_data, *nll_args))

        return r.RooAddition('sum_nll', 'Sum of all NLLs', sim_nll_list)


    def _setup_minimizer(self, nll, verbosity):
        """
        Setup the minimizer
        """
        logging.debug('setting up minimizer')
        minimizer = r.RooMinimizer(nll)
        minimizer.setMinimizerType('Minuit2')
        minimizer.setVerbose(verbosity == -1)
        minimizer.setPrintLevel(verbosity)

        return minimizer


    def _get_func_expr(self, wsp, bin_name, proto_param, param_expr):
        """
        Create the expression that can be passed to the RooWorkspace.factory
        """
        if len(param_expr) != 2:
            logging.error('Expression for proto_param \'{}\' does not have two'
                          ' sub expressions. For functional dependencies the '
                          'functional expression and the parameters for it have'
                          ' to be passed'.format(param_expr[0]))
            return 'INVALID_EXPRESSION'

        func_expr, arg_list = param_expr

        bin_data = wsp.data('full_data').reduce(
            get_bin_cut(self.bin_cut_vars, self.bins[bin_name])
        )

        # regex to check if a string contains a mean-value expression
        mean_rgx = re.compile(r'<(\w+)>')

        mean_vars = set(mean_rgx.findall(func_expr))
        mean_vals = {
            vn: self.bin_mean(wsp, vn, bin_name) for vn in mean_vars
        }

        expr = func_expr
        for var, val in mean_vals.iteritems():
            #Added parenthesis to avoid problems with negative values
            expr = expr.replace('<{}>'.format(var), '('+str(val)+')')

        return r"expr::{}('{}', {})".format(proto_param + '_' + bin_name,
                                            expr,
                                            arg_list)


    def _get_par_expr(self, wsp, bin_name, proto_param, param_expr):
        """
        Create the expression for a variable that is independent in each bin
        """
        #readout method for proto_params of the form
        #'proto_param': '[mid, low, high]'
        return '{}_{}'.format(proto_param, bin_name)+param_expr
