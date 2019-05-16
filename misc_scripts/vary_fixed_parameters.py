#!/usr/bin/env python
"""
Script to generate fit configuration files for the simultaneous fits with varied
fixed parameters
"""

import json

from itertools import product

from utils.misc_helpers import flatten, cond_mkdir_file


CHI1_VARIATIONS = (
    [{'CBn1_L': 3.2}, {'CBalpha1_L': 0.483}, {'CBalpha1_R': 2.029}],
    [{'CBn1_L': 3.4}, {'CBalpha1_L': 0.459}, {'CBalpha1_R': 1.911}],
    [{'CBn1_L': 3.0}, {'CBalpha1_L': 0.524}, {'CBalpha1_R': 2.124}],
)

CHI2_VARIATIONS = (
    [{'CBn2_L': 3.2}, {'CBalpha2_L': 0.561}, {'CBalpha2_R': 2.213}],
    [{'CBn2_L': 3.4}, {'CBalpha2_L': 0.531}, {'CBalpha1_R': 2.116}],
    [{'CBn2_L': 3.0}, {'CBalpha2_L': 0.587}, {'CBalpha1_R': 2.273}],
)

VARIATIONS = (
    CHI1_VARIATIONS, CHI2_VARIATIONS
)
import pprint


def get_all_variations():
    """
    Get all the possible variations
    """
    return (list(flatten(v)) for v in product(*VARIATIONS))


def get_out_fix_par(variation):
    """Construct a name encoding all the information from the fixed vars"""
    var_list = []
    for param in variation:
        parn, parv = list(param.iteritems())[0]
        var_list.append('{}_{:.3f}'.format(parn, parv).replace('.', 'p').replace('-', 'm'))

    return '_'.join(var_list)



def main(args):
    """Main"""
    with open(args.configfile, 'r') as conff:
        config = json.load(conff)

    all_vars = get_all_variations()

    for var in get_all_variations():
        # Update an existing setting of fixed vars
        var_list = config.get('fix_vars', [])

        for param in var:
            # Since every parameter consists of only one dict entry they can be
            # unpacked easily
            parn, parv = list(param.iteritems())[0]

            # Check if the parameter is already in the list and if so update the
            # value, otherwise put it into the list
            ipar = next((i for i, p in enumerate(var_list) if p.keys()[0] == parn), None)
            if ipar is not None:
                # Update the existing value
                var_list[ipar][parn] = parv
            else:
                var_list.append(param)

        config['fix_vars'] = var_list

        outfile = '{}_fix_{}.json'.format(args.outbase, get_out_fix_par(var))
        cond_mkdir_file(outfile)
        with open(outfile, 'w') as outf:
            json.dump(config, outf, indent=2, sort_keys=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to generate fit config'
                                     ' files with varied fixed parameters for '
                                     'the simultaneous fits')
    parser.add_argument('configfile', help='The template configuration file')
    parser.add_argument('-o', '--outbase', help='The output basename that should'
                        'be used before an auto generated part with the fixed '
                        'parameters is added to the output filename',
                        default='model_variation')

    clargs = parser.parse_args()
    main(clargs)
