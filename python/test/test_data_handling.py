#!/usr/bin/env python
"""
Tests for the data_handling module
"""

import os
import unittest
from mock import patch

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np
import pandas as pd
import pandas.testing as pdt

import logging
logging.basicConfig(level=logging.FATAL) # disable the error messages from logging

from utils.data_handling import (
    apply_selections, get_treename, list_obj, common_obj
)

class TestApplySelection(unittest.TestCase):
    def setUp(self):
        self.dummy_list =  list('abcaabcaabbccd')
        self.dummy_list2 = list('aaaaabcccbbccb')

        dfr = pd.DataFrame({'A': self.dummy_list, 'B': self.dummy_list2})
        self.dfr = pd.get_dummies(dfr, prefix=['colA', 'colB'])

    @patch('utils.data_handling.logging')
    def test_not_all_callable(self, mock_logger):
        with self.assertRaises(TypeError):
            apply_selections(pd.DataFrame(), [lambda x: x, np.array([0])])


    def test_array_selection(self):
        sel_array = np.ones(self.dfr.shape[0], dtype=bool)
        sel_array[0] = False
        sel_array[-1] = False
        sel_array[2] = False

        sel_dfr = apply_selections(self.dfr, sel_array)
        pdt.assert_frame_equal(sel_dfr, self.dfr[sel_array])


    def test_single_func_selection(self):
        sel_a = lambda df: df.colA_a == 1

        a_positions = [i for i in xrange(len(self.dummy_list))
                       if self.dummy_list[i] == 'a']

        sel_dfr = apply_selections(self.dfr, sel_a)
        pdt.assert_frame_equal(sel_dfr, self.dfr.iloc[a_positions])


    def test_list_func_selection(self):
        sel_a = lambda df: df.colA_a == 1
        sel_d = lambda df: df.colA_d == 1

        sel_dfr = apply_selections(self.dfr, (sel_a, sel_d))
        self.assertTrue(sel_dfr.empty) # no elements where a and d are 1


    def test_negate(self):
        sel_a = lambda df: df.colA_a == 1
        sel_b = lambda df: df.colB_a == 1

        nab_positions = [i for i in xrange(len(self.dummy_list))
                         if not (self.dummy_list[i] == 'a'
                                 and self.dummy_list2[i] == 'a')]

        sel_dfr = apply_selections(self.dfr, (sel_a, sel_b), negate=True)
        pdt.assert_frame_equal(sel_dfr, self.dfr.iloc[nab_positions])


class TestGetTreename(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.join(
            os.environ['CHIB_CHIC_POLFW_DIR'], 'python', 'test', 'test_data'
        )

    @patch('utils.data_handling.logging')
    def test_return_none_mult_trees(self, mock_logger):
        """Test if None is returned on multiple TTrees"""
        mult_tree_files = [
            'multiple_trees.root',
            'multiple_trees_plus_others.root'
        ]
        for filen in mult_tree_files:
            full_path = os.path.join(self.test_data_dir, filen)
            exp_warning = 'Found more than one TTrees in {}: {}'\
                          .format(full_path, ['tree1', 'tree2'])
            self.assertTrue(get_treename(full_path) is None)
            mock_logger.warning.assert_called_with(exp_warning)


    @patch('utils.data_handling.logging')
    def test_return_none_no_trees(self, mock_logger):
        """Test if None is returned on no found TTrees"""
        no_tree_files = [
            'no_tree.root',
            'no_tree_plus_others.root'
        ]
        for filen in no_tree_files:
            full_path = os.path.join(self.test_data_dir, filen)
            exp_warning = 'Found no TTrees in {}'.format(full_path)
            self.assertTrue(get_treename(full_path) is None)
            mock_logger.warning.assert_called_with(exp_warning)


    def test_return_tree_name(self):
        """Test if the returned treename is the expected one for one TTree"""
        one_tree_files = [
            'one_tree.root',
            'one_tree_plus_others.root'
        ]
        for filen in one_tree_files:
            full_path = os.path.join(self.test_data_dir, filen)
            # NOTE: 'tree1' is hardcoded here as it is in the creation script
            self.assertEqual(get_treename(full_path), 'tree1')


class TestListObj(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.join(
            os.environ['CHIB_CHIC_POLFW_DIR'], 'python', 'test', 'test_data'
        )


    def test_list_all(self):
        # NOTE: Only checking one file here
        tfile = r.TFile.Open(os.path.join(self.test_data_dir,
                                          'multiple_trees_plus_others.root'))

        elements = list(list_obj(tfile))
        for elem in ['tree1', 'tree2', 'hist1', 'hist2']:
            self.assertTrue(elem in elements)


    def test_list_type(self):
        # NOTE: Only checking one file here
        tfile = r.TFile.Open(os.path.join(self.test_data_dir,
                                          'one_tree_plus_others.root'))

        elements = list(list_obj(tfile, 'TTree'))
        self.assertEqual(elements, ['tree1'])


    def test_list_filter(self):
        tfile = r.TFile.Open(os.path.join(self.test_data_dir,
                                          'multiple_trees_plus_others.root'))

        elements = list(list_obj(tfile, filter_str='1'))
        for elem in ['tree1', 'hist1']:
            self.assertTrue(elem in elements)
        self.assertEqual(len(elements), 2)


    def test_list_type_filter(self):
        tfile = r.TFile.Open(os.path.join(self.test_data_dir,
                                          'multiple_trees_plus_others.root'))

        elements = list(list_obj(tfile, 'TTree', '1'))
        self.assertEqual(['tree1'], elements)


class TestCommonObj(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.join(
            os.environ['CHIB_CHIC_POLFW_DIR'], 'python', 'test', 'test_data'
        )


    def test_common_nofilters(self):
        files = [r.TFile.Open(os.path.join(self.test_data_dir, f)) for f in [
            'multiple_trees_plus_others.root', 'multiple_trees.root'
        ]]

        # only trees in common in this case
        objs = common_obj(files)
        for elem in ['tree1', 'tree2']:
            self.assertTrue(elem in objs)

        # nothing in common by name returns empty list
        files = [r.TFile.Open(os.path.join(self.test_data_dir, f)) for f in [
            'multiple_trees.root', 'no_tree_plus_others.root'
        ]]
        self.assertEqual([], common_obj(files))


    def test_common_filter_type(self):
        files = [r.TFile.Open(os.path.join(self.test_data_dir, f)) for f in [
            'multiple_trees_plus_others.root', 'multiple_trees.root'
        ]]
        # no type match so empty list here
        objs = common_obj(files, 'TH1')
        self.assertEqual([], objs)

        files = [r.TFile.Open(os.path.join(self.test_data_dir, f)) for f in [
            'multiple_trees_plus_others.root', 'one_tree_plus_others.root'
        ]]
        # there are trees in this file but we only want the TH1
        objs = common_obj(files, 'TH1')
        for elem in ['hist1', 'hist2']:
            self.assertTrue(elem in objs)
        self.assertEqual(len(objs), 2)


    def test_common_string_filter(self):
        files = [r.TFile.Open(os.path.join(self.test_data_dir, f)) for f in [
            'multiple_trees_plus_others.root', 'one_tree_plus_others.root'
        ]]

        objs = common_obj(files, filter_str='1')
        for elem in ['tree1', 'hist1']:
            self.assertTrue(elem in objs)
        self.assertEqual(len(objs), 2)

        objs = common_obj(files, filter_str='this_is_certainly_not_in_any_file')
        self.assertEqual([], objs)


if __name__ == '__main__':
    unittest.main()
