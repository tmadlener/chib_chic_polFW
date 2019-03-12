#!/usr/bin/env python
"""
TODO
"""

import ROOT as r
import ROOT.RooFit as rf
import numpy as np

from collections import OrderedDict

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.roofit_utils import get_var, try_factory
from utils.misc_helpers import parse_binning, get_bin_cut_root, combine_cuts
from utils.plot_helpers import setup_legend, _setup_canvas


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
    return combine_cuts([
        get_bin_cut_root(bv, *bin_borders[i]) for i, bv in enumerate(bin_vars)
    ])


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


class BinnedFitModel(object):
    """
    Fit model for doing a simultaneous fit in multiple bins
    """
    def __init__(self, config):
        """
        TODO
        """
        self.binning = np.array([
            parse_binning(config['binning'][0]),
            parse_binning(config['binning'][1])
        ])
        self.bin_vars = config['bin_vars']
        self.bins = get_bins(self.binning[0], self.binning[1], *self.bin_vars)

        self.fit_var = config['fit_variable']
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

        self.plot_config = config['plot_config']

        self.config = config


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
            for param in self.proto_params:
                wsp_var = get_var(wsp, self.proto_params[param])
                if wsp_var is not None:
                    expr = '{}_{}[{}, {}, {}]'.format(param, bin_name, wsp_var.getVal(),
                                                      wsp_var.getMin(), wsp_var.getMax())
                    success.append(try_factory(wsp, expr))
                else:
                    # TODO: actually handle expressions somehow
                    pass

            for model in self.sub_models:
                name = '{}_{}'.format(model['name'], bin_name)
                expr = model['expression'].format(name, self.fit_var)

                for param in self.proto_params:
                    expr = expr.replace(param, param + '_' + bin_name)

                success.append(try_factory(wsp, expr))

            full_expr = self.proto_fexpr
            for model in self.sub_models:
                full_expr = full_expr.replace(model['name'],
                                              model['name'] + '_' + bin_name)
            for param in self.proto_params:
                full_expr = full_expr.replace(param, param + '_' + bin_name)

            full_expr = full_expr.replace(self.full_model,
                                          self.full_model + '_' + bin_name)

            success.append(try_factory(wsp, full_expr))


        return all(success)


    def fit(self, wsp):
        """
        Import the passed data into the workspace, construct the combined
        negative log-likelihood (NLL) and run the simultaneous fit.

        Args:
            wsp (ROOT.RooWorkspace): The workspace that holds the model and into
                which the data will be imported
            data (ROOT.TTree): The (unbinned) data to which the model should be
                fitted. The binning will be done in this function
        """
        nll_args = ( # TODO: define meaning full args
            rf.NumCPU(2),
            rf.Extended(True)
        )
        sim_nll = self._create_nll(wsp, nll_args)

        minimizer = self._setup_minimizer(sim_nll)
        logging.info('starting migrad')
        minimizer.migrad()
        logging.info('starting minos')
        minimizer.minos()

        # TODO: persisting fit results, etc



    def plot(self, wsp):
        """
        Plot the fit results from the passed workspace

        Args:
            wsp (ROOT.RooWorkspace): The workspace that holds the model and the
                data as well as the fit results
        """
        fit_var = get_var(wsp, self.fit_var)
        data = wsp.data('full_data')

        cans = OrderedDict()

        for bin_name, bin_borders in self.bins.iteritems():
            frame, leg = self._plot_bin(wsp, data, bin_name, bin_borders)
            can = _setup_canvas(None) # returns a TCanvasWrapper
            frame.Draw()
            can.add_tobject(frame)
            leg.Draw()
            can.add_tobject(leg)

            cans[bin_name] = can

        # for easier debugging at the moment
        return cans

        # TODO: pulls


    def _plot_bin(self, wsp, full_data, bin_name, bin_borders):
        """Make the distribution plot for a given bin"""
        data_args = (rf.MarkerSize(0.8),)
        fit_var = get_var(wsp, self.fit_var)
        frame = fit_var.frame(rf.Bins(50))

        leg = setup_legend(*self.plot_config['legpos'])

        cut = get_bin_cut(self.bin_vars, bin_borders)
        full_data.reduce(cut).plotOn(frame, *data_args)
        full_pdf = wsp.pdf(self.full_model + '_' + bin_name)
        full_pdf.plotOn(frame, rf.LineWidth(2), rf.Name('full_pdf_curve'))

        leg.AddEntry(frame.getCurve('full_pdf_curve'), 'sum', 'l')

        for name, settings in self.components:
            full_pdf.plotOn(frame, rf.Components(name + '_' + bin_name),
                            rf.LineStyle(settings['line']),
                            rf.LineColor(settings['color']),
                            rf.LineWidth(2),
                            rf.Name(name))
            leg.AddEntry(frame.getCurve(name), settings['label'], 'l')

        # At least for debugging
        frame.SetTitle(cut)

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
            cut = get_bin_cut(self.bin_vars, bin_borders)
            bin_data = data.reduce(cut)
            bin_model = wsp.pdf(self.full_model + '_' + bin_name)
            sim_nll_list.add(bin_model.createNLL(bin_data, *nll_args))

        return r.RooAddition('sum_nll', 'Sum of all NLLs', sim_nll_list)


    def _setup_minimizer(self, nll):
        """
        Setup the minimizer
        """
        logging.debug('setting up minimizer')
        minimizer = r.RooMinimizer(nll)
        minimizer.setVerbose(False)
        minimizer.setPrintLevel(-1)

        return minimizer

