#!/usr/bin/env python
"""
Script that generates the root files that are necessary for some tests
"""

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

def _create_trees(tree_names):
    """Create the TTrees"""
    for name in tree_names:
        tree = r.TTree(name, 'test tree')
        tree.Write()


def create_root_file_trees(filename, tree_names):
    """Create a root file and place trees in it"""
    rfile = r.TFile(filename, 'recreate')
    _create_trees(tree_names)
    rfile.Write('', r.TObject.kWriteDelete)
    rfile.Close()


def create_root_file_mixed_contents(filename, tree_names, other_contents):
    """Create a root file with some TTrees and other contents"""
    rfile = r.TFile.Open(filename, 'recreate')
    _create_trees(tree_names)
    for class_name, name in other_contents:
        obj = getattr(r, class_name)()
        obj.SetName(name)
        obj.Write()

    rfile.Write('', r.TObject.kWriteDelete)
    rfile.Close()


def main():
    """Main"""
    create_root_file_trees('multiple_trees.root', ['tree1', 'tree2'])
    create_root_file_trees('one_tree.root', ['tree1'])
    create_root_file_trees('no_tree.root', [])

    create_root_file_mixed_contents('multiple_trees_plus_others.root',
                                    ['tree1', 'tree2'],
                                    [('TH1D', 'hist1'), ('TH1D', 'hist2')])
    create_root_file_mixed_contents('one_tree_plus_others.root',
                                    ['tree1'],
                                    [('TH1D', 'hist1'), ('TH1D', 'hist2')])
    create_root_file_mixed_contents('no_tree_plus_others.root',
                                    [],
                                    [('TH1D', 'hist1'), ('TH1D', 'hist2')])


if __name__ == '__main__':
    main()
