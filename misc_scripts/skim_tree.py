#!/usr/bin/env python
"""
Script to skim a TTree from the command line
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

def main(args):
    """Main"""
    infile = r.TFile.Open(args.infile)
    intree = infile.Get(args.tree)

    if args.remove:
        for branch in args.remove:
            if not branch in args.cut:
                if args.verbosity > 1:
                    print('Disabling branch \'{}\''.format(branch))
                intree.SetBranchStatus(branch, 0)

    outfile = r.TFile(args.outfile, 'recreate')
    outtree = intree.CopyTree(args.cut)

    if args.verbosity > 0:
        print('Number of events: Input: {}, Output: {}'
              .format(intree.GetEntries(), outtree.GetEntries()))
        print('Number of branches: Input: {}, Output: {}'
              .format(len(intree.GetListOfBranches()),
                      len(outtree.GetListOfBranches())))

    outtree.Write()
    outfile.Write('', r.TObject.kWriteDelete)
    outfile.Close()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to skim a TTree and '
                                     'write it into a new file')
    parser.add_argument('infile', help='Input file')
    parser.add_argument('outfile', help='Output file')
    parser.add_argument('-t', '--tree', help='Name of the TTree to skim',
                        default='jpsi_tuple')
    parser.add_argument('-c', '--cut', help='Cut string that should be applied'
                        ' while skimming', default='')
    parser.add_argument('-r', '--remove', help='Space separated list of branches'
                        ' that should be removed completely from the output '
                        'TTree', nargs='+', default=[])
    parser.add_argument('-v', '--verbosity', default=1, help='Set the verbosity'
                        ' of the process, 0 to disable all output')

    clargs = parser.parse_args()
    main(clargs)
