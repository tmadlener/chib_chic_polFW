#!/usr/bin/env python
"""
Module holding some useful functionality to do selections on pandas.DataFrames

TODO: Make pylint happy (or at least happier) by cleaning up code
"""

import numpy as np

from decorator import decorator

import logging
logger = logging.getLogger()

from utils.misc_helpers import (
    select_bin, get_full_trigger, _get_var, deprecated_soon, make_iterable
)
from utils.data_handling import apply_selections


@decorator
def deprecated_lambda_call(func, *args, **kwargs):
    """
    Decorator that prints out a warning when one of the methods that now return
    a functor is called with the old "style" where they returned a lambda with
    the reight interface directly.
    """
    functor = func(*args, **kwargs)
    if len(args) > 0:
        logger.warn('Using \'{}\' as lambda function will soon be deprecated. '
                     'It is suggested to instead switch to instantiating it '
                     'and use it as a functor.'.format(func.__name__))
        # only pass the dataframe which should be the first arg
        return functor(args[0])

    return functor


def get_cut_funcs(coords):
    """
    Get a cut function from the passed coordinates in the pt-eta plane
    """
    cut_funcs = []

    max_eta = max([c[1] for c in coords])
    # min_pt = min([c[0] for c in coords])

    # contains function, checking if x is somewhere in between l and h
    # without the need to order lower and upper bounds
    cont = lambda x, l, h: ((x > l) & (x < h)) | ((x > h) & (x < l))

    for c1, c2 in zip(coords[:-1], coords[1:]):
        pt1, eta1 = c1
        pt2, eta2 = c2
        # sorting to not have to worry about the order in which the coordinates
        # have been defined
        # 3 cases:
        if pt1 == pt2:
            # to capture the values of pt1, etc by "value" create a new scope in
            # a function. See, e.g.: https://stackoverflow.com/a/13355291
            def pt_cut(pt, eta, pt1=pt1, pt2=pt2, eta1=eta1, eta2=eta2):
                return (pt > pt1) & cont(eta, eta1, eta2)
            cut_funcs.append(pt_cut)

        elif eta1 == eta2:
            def eta_cut(pt, eta, pt1=pt1, pt2=pt2, eta1=eta1, eta2=eta2):
                return (eta > eta1) & cont(pt, pt1, pt2) & (eta < max_eta)
            cut_funcs.append(eta_cut)

        elif pt1 != pt2 and eta1 != eta2:
            def lin_int(pt, eta, pt1=pt1, pt2=pt2, eta1=eta1, eta2=eta2):
                int_eta = eta1 + (pt - pt1) * (eta2 - eta1) / (pt2 - pt1)
                return (eta > int_eta) & cont(eta, eta1, eta2)
            cut_funcs.append(lin_int)

        else:
            print('Something went wrong while defining the selection')

    return cut_funcs


def pt_eta_sel(pt, eta, cuts):
    """Select in the 2D pt-eta plane applying the specified cuts"""
    # since get_cut_funcs returns multiple functions whose product is the or of all functions,
    # start with all events as not selected
    decision = np.zeros(pt.shape, dtype=bool)

    cut_funcs = get_cut_funcs(cuts)
    for cut in cut_funcs:
        decision |= cut(pt, eta)

    return decision


def get_gen_name(name, gen=False):
    """Helper function to get the gen level variable name from the passed variable"""
    if gen:
        return 'gen_' + name
    else:
        return name


def single_muon_sel(cuts, gen=False):
    """Apply a single muon selection (using the AND of both single muons)"""
    class SingleMuonSel(object):
        """Internal helper class to easier handle the requires case"""
        def __init__(self, cuts, gen=False):
            self.cuts = cuts
            self.gen = gen
            self.requires = [get_gen_name(v, gen) for v in ['muNPt', 'muPPt',
                                                            'muNEta', 'muPEta']]

        def __call__(self, dfr):
            mu_name = get_gen_name('mu', self.gen)
            mup_pt, mun_pt = [mu_name + chg + 'Pt' for chg in ['P', 'N']]
            mup_eta, mun_eta = [mu_name + chg + 'Eta' for chg in ['P', 'N']]

            return pt_eta_sel(dfr[mup_pt], dfr[mup_eta].abs(), self.cuts) & \
                pt_eta_sel(dfr[mun_pt], dfr[mun_eta].abs(), self.cuts)

    return SingleMuonSel(cuts, gen)


def photon_sel_(cuts, gen=False):
    """Apply a photon selection
    TODO: rename after deprecation
    """
    class PhotonSel(object):
        """Internal helper class to easier handle the requires case"""
        def __init__(self, cuts, gen=False):
            self.cuts = cuts
            self.gen = gen
            self.requires = [get_gen_name(v, gen) for v in ['photonPt',
                                                            'photonEta']]

        def __call__(self, dfr):
            phot_name = get_gen_name('photon', self.gen)
            return pt_eta_sel(_get_var(dfr, phot_name + 'Pt'),
                              _get_var(dfr, phot_name + 'Eta').abs(),
                              self.cuts)

    return PhotonSel(cuts, gen)


@deprecated_soon('photon_sel_')
def photon_sel(df, cuts, gen=False):
    """TODO: switch to returning a functor after deprecation and deprecate
    'photon_sel_' then
    """
    return photon_sel_(cuts, gen)(df)


TRIGGER_BIT_MAP = {
    # 2012
    'HLT_Dimuon8_Jpsi': 1,
    # 2016
    'HLT_Dimuon10_Jpsi_Barrel': 1 << 0,
    'HLT_Dimuon16_Jpsi': 1 << 1,
    'HLT_Dimuon20_Jpsi': 1 << 2,
    # 2017
    'HLT_Dimuon20_Jpsi_Barrel_Seagulls': 1 << 0,
    'HLT_Dimuon25_Jpsi': 1 << 1,
}

def trigger_sel_(trigger='HLT_Dimuon8_Jpsi'):
    """Apply the trigger selection (reco MC and data only)"""
    class TriggerSel(object):
        """Internal helper class to easier handle the requires case"""
        def __init__(self, trig):
            self.requires = ['trigger']
            self.trigger_bit = TRIGGER_BIT_MAP[get_full_trigger(trig)]

        def __call__(self, dfr):
            return (dfr.trigger & self.trigger_bit) == self.trigger_bit

    return TriggerSel(trigger)


@deprecated_soon('trigger_sel_')
def trigger_sel(df, trigger='HLT_Dimuon8_Jpsi'):
    """Apply the trigger selection (reco and data only)"""
    full_trigger = get_full_trigger(trigger)
    trigger_bit = TRIGGER_BIT_MAP[full_trigger]
    # NOTE: assuming that the trigger decision is stored in 'trigger' column
    return (df.trigger & trigger_bit) == trigger_bit
trigger_sel.requires = ['trigger']

def prompt_sel_(ctau_sigma=2.5):
    """Select the prompt events using a lifetime significance cut"""
    class PromptSel(object):
        """Internal helper class to handle the requires case"""
        def __init__(self, ctau_sig):
            self.sigma = ctau_sig
            self.requires = ['Jpsict', 'JpsictErr']

        def __call__(self, dfr):
            return (dfr.Jpsict.abs() / dfr.JpsictErr) < self.sigma

    return PromptSel(ctau_sigma)


@deprecated_soon('prompt_sel_')
def prompt_sel(dfr, ctau_sigma=2.5):
    """
    Select prompt events using a lifetime significance cut
    """
    return (dfr.Jpsict.abs() / dfr.JpsictErr) < ctau_sigma
prompt_sel.requires = ['Jpsict', 'JpsictErr']


@deprecated_soon('vtx_prob_sel_')
def vtx_prob_sel(df, prob=0.01):
    """Select only the events with vtx probability larger than 0.01"""
    return df.vtxProb > prob
vtx_prob_sel.requires = ['vtxProb']


def vtx_prob_sel_(prob=0.01):
    """
    Create a vertex probability selection functor that selects events that have
    a vtx probability greater than the passed value
    """
    sel_func = lambda d: d.vtxProb > prob
    sel_func.requires = ['vtxProb']
    return sel_func


def deta_sel(df, deta_max=0.015):
    """
    Select only events for which the generated and reconstruced photon eta match
    """
    return (df.photonEta - df.gen_photonEta).abs() < deta_max
deta_sel.requires = ['photonEta', 'gen_photonEta']

def chic_mass_sel(df, min_mass=3.325, max_mass=3.725):
    """Select only those events with a chic mass between min_mass and max_mass"""
    return select_bin('chicMass', min_mass, max_mass)(df)
chic_mass_sel.requires = ['chicMass']

def jpsi_kin_sel_(min_pt=8, max_pt=20, max_rap=1.2, min_rap=0, gen=False):
    """Kinematic selection of jpsi
    TODO: Rename after deprecation
    """
    class JpsiKinSel(object):
        """Internal helper class to help with the """
        def __init__(self, min_pt, max_pt, max_rap, min_rap, gen):
            self.min_pt = min_pt
            self.max_pt = max_pt
            self.max_rap = max_rap
            self.min_rap = min_rap
            self.gen = gen
            self.requires = [get_gen_name(v, gen) for v in ['JpsiPt',
                                                            'JpsiRap']]


        def __call__(self, dfr):
            jpsi_name = get_gen_name('Jpsi', self.gen)

            pt_sel = select_bin(jpsi_name + 'Pt', self.min_pt, self.max_pt)
            rap_sel = select_bin('abs(' + jpsi_name + 'Rap)',
                                 self.min_rap, self.max_rap)

            return pt_sel(dfr) & rap_sel(dfr)

    return JpsiKinSel(min_pt, max_pt, max_rap, min_rap, gen)


@deprecated_lambda_call
def jpsi_kin_sel(*args, **kwargs):
    """TODO: remove after deprecation"""
    return jpsi_kin_sel_(min_pt=kwargs.pop('min_pt', 8),
                         max_pt=kwargs.pop('max_pt', 20),
                         max_rap=kwargs.pop('max_rap', 1.2),
                         gen=kwargs.pop('gen', False))


@deprecated_soon('state_sel')
def chic_state_sel(df, state):
    """(On MC) select either chic1 or chic2 using the pdgId"""
    pdgId = {'chic1': 20443, 'chic2': 445}
    return df.pdgId == pdgId[state]
chic_state_sel.requires = ['pdgId']


def state_sel(state):
    """
    Select a given state by string (using the pdgId in the end)
    """
    class StateSel(object):
        """Internal helper class"""
        def __init__(self, state):
            self.pdgId = {
                'chic1': 20443,
                'chic2': 445
            }[state]
            self.requires = ['pdgId']


        def __call__(self, dfr):
            return dfr.pdgId == self.pdgId


    return StateSel(state)


def flat_pt(pt, eta):
    """Helper function to return a flat pt-eta cut in 2D coords"""
    return (pt, 0), (pt, eta)


def fiducial_cuts():
    """Get the coordinates representing the standard fiducial cuts"""
    coords = (
        (3.0, 1.6), (3.0, 1.4), (3.5, 1.4), (3.5, 1.2), (4.5, 1.2),
        (4.5, 0)
    )
    return coords


def loose_cuts():
    """Loosest cut that does not have any acceptance holes in pt-eta"""
    coords = (
        (2.5, 1.6), (3.5, 1.2), (3.5, 0)
    )

    return coords


def gen_filter_cuts():
    """Gen level filter cut values"""
    return (
        (1.0, 2.5), (1.0, 1.0), (2.5, 1.0), (2.5, 0)
    )


@deprecated_lambda_call
def fiducial_muon_sel(*args, **kwargs):
    """
    The fiducial muon selection
    """
    return single_muon_sel(fiducial_cuts(), gen=kwargs.pop('gen', False))


@deprecated_lambda_call
def loose_muon_sel(*args, **kwargs):
    """The loose muon selection"""
    return single_muon_sel(loose_cuts(), gen=kwargs.pop('gen', False))


def gen_filter_sel(gen=True):
    """The generation level filter muon selection"""
    return single_muon_sel(gen_filter_cuts(), gen)


def get_n_events(data, selections=None, weight=None):
    """
    Get the number of events in the dataframe surviving the passed selection

    If weight is not None the corresponding column will be used as weights,
    unless weight is a function taking the dataframe as only input, than the sum
    of the array returned by that function call will be returned
    """
    sel_data = apply_selections(data, selections)
    if weight is None:
        return sel_data.shape[0]
    weights = _get_var(sel_data, weight)
    return weights.sum()


def toy_reco(df, eff_name='gamma_eff_sm'):
    """
    Select only those events for which the efficiency is larger than 0
    """
    return df.loc[:, eff_name] > 0


def all_sel():
    """
    Selection that selects all events in a way that still satisfies the requires
    helpers
    """
    class AllSel(object):
        """Internal helper class"""
        def __init__(self):
            self.requires = []

        def __call__(self, dfr):
            return np.ones(dfr.shape[0], dtype=bool)

    return AllSel()


def random_sel(n_events=None, fraction=None):
    """
    Randomly select events from the sample.

    NOTE: When applied in combination with other selections, this selection will
    use the unselected sample, so that the results might not be as expected.

    The two arguments are mutually exclusive

    Args:
        n_events (int): Select exactly this number of events
        fraction (float): Select this fraction of events (approximately)
    """
    if n_events is None and fraction is None:
        err = True
        logger.error('Need either n_events or fraction as argument')
    if n_events is not None and fraction is not None:
        logger.error('n_events and fraction are mutually exclusive arguments')

    if n_events is not None:
        def sel_func(data):
            """Function doing the selection"""
            if n_events > data.shape[0]:
                logger.warn('n_events = {} > number of data events {}. Will '
                             'select all events'.format(n_events, data.shape[0]))
                return np.ones(data.shape[0], dtype=bool)

            sel = np.append(
                np.ones(n_events, dtype=bool),
                np.zeros(data.shape[0] - n_events, dtype=bool))
            np.random.shuffle(sel)

            return sel

        func = sel_func

    if fraction is not None:
        if fraction < 0 or fraction > 1:
            logger.warn('Fraction has to be between 0 and 1')
        func = lambda d: np.random.uniform(0, 1, d.shape[0]) < fraction

    func.requires = []
    return func


def collect_requirements(selections):
    """Collect the list of variables that is needed for the selections"""
    if selections is None:
        return []
    variables = set()
    for selection in make_iterable(selections):
        if hasattr(selection, 'requires'):
            for req in selection.requires:
                variables.add(req)
        else:
            sel_name = ''
            if hasattr(selection, '__name__'):
                sel_name = selection.__name__
            elif hasattr(selection, '__class__'):
                sel_name = selection.__class__

            logger.warning('\'{}\' does not have a requires field, possibly '
                            'cannot get the necessary variables'
                            .format(sel_name))

    return list(variables)
