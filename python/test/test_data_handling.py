#!/usr/bin/env python
"""
Tests for the data_handling module
"""

import os
import unittest

import numpy as np
import pandas as pd

import logging
logging.basicConfig(level=logging.FATAL) # disable the error messages from logging

from utils.data_handling import apply_selections, get_treename

class TestApplySelection(unittest.TestCase):
    def setUp(self):
        self.dummy_list =  list('abcaabcaabbccd')
        self.dummy_list2 = list('aaaaabcccbbccb')

        dfr = pd.DataFrame({'A': self.dummy_list, 'B': self.dummy_list2})
        self.dfr = pd.get_dummies(dfr, prefix=['colA', 'colB'])

    def test_not_all_callable(self):
        with self.assertRaises(TypeError):
            apply_selections(pd.DataFrame(), [lambda x: x, np.array([0])])


    def test_array_selection(self):
        sel_array = np.ones(self.dfr.shape[0], dtype=bool)
        sel_array[0] = False
        sel_array[-1] = False
        sel_array[2] = False

        sel_dfr = apply_selections(self.dfr, sel_array)

        # check shape
        self.assertEqual(sel_dfr.shape[0], np.sum(sel_array))
        # ... and contents
        self._compare_dfrs(self.dfr[sel_array], sel_dfr)


    def test_single_func_selection(self):
        sel_a = lambda df: df.colA_a == 1

        a_positions = [i for i in xrange(len(self.dummy_list))
                       if self.dummy_list[i] == 'a']

        sel_dfr = apply_selections(self.dfr, sel_a)
        self.assertEqual(sel_dfr.shape[0], self.dummy_list.count('a'))
        self._compare_dfrs(self.dfr.iloc[a_positions], sel_dfr)


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
        self.assertEqual(sel_dfr.shape[0], len(nab_positions))
        self._compare_dfrs(self.dfr.iloc[nab_positions], sel_dfr)


    def _compare_dfrs(self, dfr1, dfr2):
        """
        Compare two dataframes
        """
        # check basic shape
        self.assertEqual(dfr1.shape, dfr2.shape)
        self.assertTrue(np.all(dfr1 == dfr2))


class TestGetTreename(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.join(
            os.environ['CHIB_CHIC_POLFW_DIR'], 'python', 'test', 'test_data'
        )


    def test_return_none_mult_trees(self):
        """Test if None is returned on multiple TTrees"""
        mult_tree_files = [
            'multiple_trees.root',
            'multiple_trees_plus_others.root'
        ]
        for filen in mult_tree_files:
            full_path = os.path.join(self.test_data_dir, filen)
            self.assertTrue(get_treename(full_path) is None)


    def test_return_none_no_trees(self):
        """Test if None is returned on no found TTrees"""
        no_tree_files = [
            'no_tree.root',
            'no_tree_plus_others.root'
        ]
        for filen in no_tree_files:
            full_path = os.path.join(self.test_data_dir, filen)
            self.assertTrue(get_treename(full_path) is None)


    def test_return_tree_name(self):
        """Test if the returned treename is the expected one for one TTree"""
        one_tree_files = [
            'one_tree.root',
            'one_tree_plus_others.root'
        ]
        for filen in one_tree_files:
            full_path = os.path.join(self.test_data_dir, filen)
            # NOTE: 'tree' is hardcoded here as it is in the creation scrip:t
            self.assertEqual(get_treename(full_path), 'tree')


if __name__ == '__main__':
    unittest.main()
