#!/usr/bin/env python
"""
Script to preselect events (e.g. to feed them to the costh binned mass fits
afterwards)
"""

from utils.data_handling import (
    get_dataframe, apply_selections, store_dataframe
)
import utils.selection_functions as sf

VARIABLES = ['trigger', 'Jpsi{Pt,Rap}', '{photon,mu{N,P}}{Pt,Eta}', 'vtxProb',
             'Jpsict{,Err}', 'costh_HX', 'chicMass']

def get_jpsi_sel(jpsi_sel_str):
    """
    Get the jpsi selection function from the passed string
    """
    min_pt, max_pt, max_rap = [float(v) for v in jpsi_sel_str.split(':')]

    return lambda d: sf.jpsi_kin_sel(d, min_pt, max_pt, max_rap)


def main(args):
    """Main"""
    data = get_dataframe(args.infile, columns=VARIABLES, where='trigger > 0')

    selections = (
        sf.loose_muon_sel,
        lambda d: sf.trigger_sel(d, args.trigger),
        sf.vtx_prob_sel,
        get_jpsi_sel(args.jpsi),
        lambda d: sf.photon_sel(d, sf.flat_pt(0.4, 1.5)),
        sf.chic_mass_sel # not strictly necessary?
    )

    sel_data = apply_selections(data, selections)
    store_dataframe(sel_data, args.outfile, 'chic_tuple')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to preselect ')
    parser.add_argument('infile', help='inputfile')
    parser.add_argument('outfile', help='preselected file')
    parser.add_argument('-t', '--trigger', help='trigger',
                        default='Dimuon8_Jpsi')
    parser.add_argument('-j', '--jpsi', default='8:20:1.2',
                        help='jpsi selection. format: minpt:maxpt:maxrap')

    clargs = parser.parse_args()
    main(clargs)
