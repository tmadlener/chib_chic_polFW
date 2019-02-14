#!/usr/bin/env python
"""
Tests for misc_helpers functions
"""

import unittest
from mock import patch

import re

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
import pandas as pd
import numpy as np
import numpy.testing as npt
import pandas.util.testing as pdt

from utils.data_handling import apply_selections

import utils.misc_helpers as mh


def same_binning(bin1, bin2):
    """
    Check if binnings are same taking into account floating point things
    """
    for ibin, borders in enumerate(bin1):
        # NOTE: this does not need to bu super accurate, thus if the borders
        # are within 0.005 we are happy enough (also because otherwise it
        # gets harder to setup the tests)
        npt.assert_allclose(borders, bin2[ibin], rtol=0, atol=5e-3)


def equal_bin_contents(test_class, sel_dfr, binning, var):
    """
    Check if the binning is indeed equally poplated for a given variable
    """
    from utils.data_handling import apply_selections as app_sel

    select_var = lambda df, low, high: (var(df) > low) & (var(df) <= high)
    n_elem = [app_sel(sel_dfr,
                      lambda df: select_var(df, b[0], b[1])).shape[0]
              for b in binning]

    # The criterion for equally binned is half a permille for us (more or less aribtrary)
    unique_occ = np.unique(n_elem)
    test_class.assertLess((np.max(unique_occ) - np.min(unique_occ)) / float(sel_dfr.shape[0]),
                          5e-4)


class TestFlatten(unittest.TestCase):
    def test_already_flat(self):
        in_list = list(xrange(10))
        flat_list = list(mh.flatten(in_list))
        self.assertEqual(flat_list, in_list)


    def test_nested_list(self):
        in_list = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10]]
        flat_list = list(mh.flatten(in_list))
        self.assertEqual(flat_list, range(1, 11))

        in_list = [1, [[2, 3], [4, 5]], 6, 7, [8, 9, [10]]]
        self.assertEqual(flat_list, range(1, 11))


class TestChunks(unittest.TestCase):
    def test_list_smaller_chunk_size(self):
        in_list = range(3)
        chunk_size = 6
        chunks = list(mh.chunks(in_list, chunk_size))

        self.assertEqual(chunks[0], in_list)


    def test_list_larger_chunk_size_exact_possible(self):
        in_list = range(100)
        chunk_size = 20
        chunks = list(mh.chunks(in_list, chunk_size))

        self.assertEqual(len(chunks), 5)
        for chunk in chunks:
            self.assertEqual(len(chunk), 20)


    def test_list_larger_chunk_size_not_exact_possible(self):
        in_list = range(100)
        chunk_size = 27
        chunks = list(mh.chunks(in_list, chunk_size))

        self.assertEqual(len(chunks), 4)
        for chunk in chunks[:-1]:
            self.assertEqual(len(chunk), 27)
        self.assertEqual(len(chunks[-1]), 100 % 27)


    def test_chunks_preserves_order(self):
        in_list = range(100)
        chunk_size = 20

        self.assertEqual(list(mh.flatten(mh.chunks(in_list, chunk_size))),
                         in_list)


class TestGetCosthBinning(unittest.TestCase):
    def setUp(self):
        np.random.seed(12345)
        self.dfr = pd.DataFrame({
            'costh_lin': np.random.uniform(-1, 1, 500000),
            'costh_squ': np.random.uniform(-1, 1, 500000)**2,
        })


    def test_binning_no_selection(self):
        # first test linear case
        n_bins = 5
        w_bin = self.dfr.costh_lin.abs().max() / n_bins # width of bin assuming a uniform distribution
        self.dfr['costh_HX'] = self.dfr.costh_lin
        lin_bins = mh.get_costh_binning(self.dfr, n_bins)
        self.assertEqual(len(lin_bins), n_bins)
        same_binning(lin_bins, [(i*w_bin, (i+1)*w_bin) for i in xrange(n_bins)])
        equal_bin_contents(self, self.dfr, lin_bins, lambda d: d.costh_lin.abs())


        self.dfr['costh_HX'] = self.dfr.costh_squ
        squ_bins = mh.get_costh_binning(self.dfr, n_bins)
        self.assertEqual(len(squ_bins), n_bins)

        same_binning(squ_bins, [((i*w_bin)**2, ((i+1)*w_bin)**2) for i in xrange(n_bins)])
        equal_bin_contents(self, self.dfr, squ_bins, lambda d: d.costh_squ.abs())


    def test_binning_selection(self):
        sel = self.dfr.costh_lin >= 0
        n_bins = 5
        w_bin = self.dfr[sel].costh_lin.abs().max() / n_bins # width of bin assuming a uniform distribution

        self.dfr['costh_HX'] = self.dfr.costh_lin
        lin_bins = mh.get_costh_binning(self.dfr, n_bins, selection=sel)
        self.assertEqual(len(lin_bins), n_bins)
        same_binning(lin_bins, [(i*w_bin, (i+1)*w_bin) for i in xrange(n_bins)])

        equal_bin_contents(self, self.dfr[sel], lin_bins, lambda d: d.costh_lin.abs())

        self.dfr['costh_HX'] = self.dfr.costh_squ
        squ_bins = mh.get_costh_binning(self.dfr, n_bins, selection=sel)
        self.assertEqual(len(squ_bins), n_bins)
        same_binning(squ_bins, [((i*w_bin)**2, ((i+1)*w_bin)**2) for i in xrange(n_bins)])
        equal_bin_contents(self, self.dfr[sel], squ_bins, lambda d: d.costh_squ.abs())


class TestGetEquiPopBins(unittest.TestCase):
    def setUp(self):
        np.random.seed(12345)
        self.dfr = pd.DataFrame({
            'x': np.random.uniform(-1, 1, 500000),
            'y': np.random.uniform(0, 1, 500000)
        })


    def test_get_equi_pop_bins(self):
        n_bins = 5
        w_bin = (self.dfr.y.max() - self.dfr.y.min()) / n_bins

        bins = mh.get_equi_pop_bins(self.dfr, lambda df: df.y, n_bins)
        self.assertEqual(len(bins), n_bins)
        same_binning(bins, [(i*w_bin, (i+1)*w_bin) for i in xrange(n_bins)])
        equal_bin_contents(self, self.dfr, bins, lambda d: d.y)

        w_bin = (self.dfr.x.max() - self.dfr.x.min()) / n_bins
        min_bin = self.dfr.x.min()
        bins = mh.get_equi_pop_bins(self.dfr, lambda df: df.x, n_bins)
        self.assertEqual(len(bins), n_bins)
        same_binning(bins, [(min_bin + i*w_bin,
                                   min_bin + (i+1)*w_bin) for i in xrange(n_bins)])
        equal_bin_contents(self, self.dfr, bins, lambda d: d.x)


    def test_get_equi_pop_bins_abs(self):
        n_bins = 5
        w_bin = (self.dfr.x.abs().max() - self.dfr.x.abs().min()) / n_bins

        bins = mh.get_equi_pop_bins(self.dfr, lambda df: df.x.abs(), n_bins)
        self.assertEqual(len(bins), n_bins)
        same_binning(bins, [(i*w_bin, (i+1)*w_bin) for i in xrange(n_bins)])
        equal_bin_contents(self, self.dfr, bins, lambda d: d.x.abs())


    def test_equi_pop_bins_non_divisor(self):
        n_bins = 13
        w_bin = (self.dfr.y.max() - self.dfr.y.min()) / n_bins

        bins = mh.get_equi_pop_bins(self.dfr, lambda d: d.y, n_bins)
        self.assertEqual(len(bins), n_bins)
        same_binning(bins, [(i*w_bin, (i+1)*w_bin) for i in xrange(n_bins)])
        equal_bin_contents(self, self.dfr, bins, lambda d: d.y)


class TestGetBinMeans(unittest.TestCase):
    def setUp(self):
        np.random.seed(12345)
        self.dfr = pd.DataFrame({
            'x': np.random.uniform(0, 1, 500000)
        })

    def test_get_bin_means(self):
        n_bins = 5
        w_bin = (self.dfr.x.max() - self.dfr.x.min()) / n_bins
        bin_centers = [w_bin / 2 + i * w_bin for i in xrange(n_bins)]

        bins = mh.get_equi_pop_bins(self.dfr, lambda d: d.x, n_bins)
        means = mh.get_bin_means(self.dfr, lambda d: d.x, bins)

        # if we get to the half percent level we are satisfied
        npt.assert_allclose(means, bin_centers, atol=5e-3, rtol=0)


    def test_get_bin_means_selection(self):
        sel = self.dfr.x < 0.8 # should cut off roughly the top 20 %
        n_bins = 4
        w_bin = (self.dfr[sel].x.max() - self.dfr[sel].x.min()) / n_bins
        bin_centers = [w_bin / 2 + i * w_bin for i in xrange(n_bins)]

        bins = mh.get_equi_pop_bins(self.dfr[sel], lambda d: d.x, n_bins)
        means = mh.get_bin_means(self.dfr, lambda d: d.x, bins, sel)
        npt.assert_allclose(means, bin_centers, atol=5e-3, rtol=0)


    def test_get_bin_means_weighted(self):
        # add weights to dataframe for easier handling
        self.dfr['w'] = np.random.uniform(0, 1, 500000)
        binning = [(lo, hi) for lo, hi in zip(np.linspace(0, 1, 5)[:-1],
                                              np.linspace(0, 1, 5)[1:])]

        exp_means = []
        for lo, hi in binning:
            sel_dfr = apply_selections(self.dfr, lambda d: (d.x > lo) & (d.x < hi))
            w_sum = sel_dfr.w.sum()
            wx_sum = np.sum(sel_dfr.x * sel_dfr.w)
            exp_means.append(wx_sum / w_sum)

        means = mh.get_bin_means(self.dfr, 'x', binning, weights=self.dfr.w)
        npt.assert_allclose(means, np.array(exp_means))



class TestUniqueWKey(unittest.TestCase):
    def test_unique_w_key(self):
        tuple_l = [(1, '1'), (2, '2'), (3, '1'), (4, '4'), (5, '5')]
        unique_l = mh.unique_w_key(tuple_l, lambda e: e[1])
        self.assertEqual(unique_l, [(1, '1'), (2, '2'), (4, '4'), (5, '5')])

        unique_l = mh.unique_w_key(tuple_l, lambda e: e[0])
        self.assertEqual(unique_l, tuple_l)

        unique_l = mh.unique_w_key(tuple_l, lambda e: e[0] % 3)
        self.assertEqual(unique_l, [(1, '1'), (2, '2'), (3, '1')])


class TestGetNpFromTMatrix(unittest.TestCase):
    def test_get_symmetric(self):
        n_dim = 3

        # fill a matrix
        matrix = r.TMatrixTSym(float)(n_dim)
        # create a symmetric random matrix
        m_matrix = np.random.uniform(size=(n_dim, n_dim))
        m_matrix = (m_matrix + m_matrix.T) / 2

        # TODO: There should be a way to write directly to the TMatrixT buffer,
        # but I currently can't seem to get the types to match
        for i in xrange(n_dim):
            for j in xrange(n_dim):
                matrix[i][j] = m_matrix[i][j]

        # some (trivial) setup tests
        self.assertEqual(len(m_matrix), matrix.GetNrows())
        self.assertEqual(len(m_matrix[0]), matrix.GetNcols())

        np_matrix = mh.get_np_from_tmatrix(matrix)

        self.assertEqual(m_matrix.shape, np_matrix.shape)
        npt.assert_allclose(np_matrix, m_matrix, atol=1e-6, rtol=0)


    def test_arbitrary_matrix(self):
        n_cols = 3
        n_rows = 4

        matrix = r.TMatrixT(float)(n_rows, n_cols)
        m_matrix = np.random.uniform(size=(n_rows, n_cols))
        for i in xrange(n_rows):
            for j in xrange(n_cols):
                matrix[i][j] = m_matrix[i][j]

        # some (trivial) setup tests
        self.assertEqual(len(m_matrix), matrix.GetNrows())
        self.assertEqual(len(m_matrix[0]), matrix.GetNcols())

        np_matrix = mh.get_np_from_tmatrix(matrix)
        self.assertEqual(m_matrix.shape, np_matrix.shape)
        npt.assert_allclose(np_matrix, m_matrix, atol=1e-6, rtol=0)


class TestGetBinCenters(unittest.TestCase):
    def test_even_binning(self):
        binning = np.linspace(0, 10, 11)
        exp_centers = np.array([i + 0.5 for i in xrange(10)])
        bcenters = mh.get_bin_centers(binning)

        npt.assert_allclose(bcenters, exp_centers)

    # def test_arb_binning(self):
    #     # TODO
    #     pass


class Test_GetVar(unittest.TestCase):
    def setUp(self):
        self.test_df = pd.DataFrame({
            'col1': np.random.uniform(10000, -1, 1),
            'col2': np.random.uniform(10000, 0, 1)
        })


    def test_get_var_func(self):
        """Test if using a function works correctly"""
        var = mh._get_var(self.test_df, lambda d: d.col2)
        pdt.assert_series_equal(var, self.test_df.col2)

        var = mh._get_var(self.test_df, lambda d: d.col1**2)
        pdt.assert_series_equal(var, self.test_df.col1**2)


    def test_get_var_string(self):
        """Test that a string or a list of string gets handled correctly"""
        var = mh._get_var(self.test_df, 'col2')
        pdt.assert_series_equal(var, self.test_df.col2)

        var = mh._get_var(self.test_df, ['col1', 'col2'])
        pdt.assert_frame_equal(var, self.test_df)


    def test_get_var_var(self):
        """Test that when already passed a variable the same is returned"""
        var = mh._get_var(self.test_df, self.test_df.col1)
        pdt.assert_series_equal(var, self.test_df.col1)


    def test_get_var_np_func(self):
        """Test that passing an additional numpy function works"""
        var = mh._get_var(self.test_df, 'col2', np.sqrt)
        pdt.assert_series_equal(var, np.sqrt(self.test_df.col2))

        var = mh._get_var(self.test_df, 'col1', 'abs')
        pdt.assert_series_equal(var, self.test_df.col1.abs())


class TestParseBinning(unittest.TestCase):
    """
    Some tests to illustrate the parse binning function and to check if the
    regexes used to parse are working properly
    """
    def test_parse_list_of_floats(self):
        binning = mh.parse_binning('1,2,3,4,5,6')
        npt.assert_allclose(binning, np.arange(1, 7))

        binning = mh.parse_binning('1.2, 3.4, 5.6,7.8')
        npt.assert_allclose(binning, np.array([1.2, 3.4, 5.6, 7.8]))


    def test_parse_linspace(self):
        binning = mh.parse_binning('1:10,25')
        npt.assert_allclose(binning, np.linspace(1, 10, 25))

        binning = mh.parse_binning('-3.8:4.,30')
        npt.assert_allclose(binning, np.linspace(-3.8, 4.0, 30))


    def test_parse_arange(self):
        binning = mh.parse_binning('1:10:2')
        npt.assert_allclose(binning, np.arange(1, 10, 2))

        binning = mh.parse_binning('-3.1:4.2:0.2')
        npt.assert_allclose(binning, np.arange(-3.1, 4.2, 0.2))


    @patch('utils.misc_helpers.logging')
    def test_parse_fail(self, mock_logger):
        exp_err = 'Cannot handle \'3,4:2\' in trying to parse a binning'
        binning = mh.parse_binning('3,4:2')
        mock_logger.error.assert_called_with(exp_err)
        npt.assert_equal(binning, np.array([]))


        exp_err = 'Could not parse binning from \'3.a, 3.3\' because '\
                  '\'invalid literal for float(): 3.a\''
        binning = mh.parse_binning('3.a, 3.3')
        mock_logger.error.assert_called_with(exp_err)
        npt.assert_equal(binning, np.array([]))


class TestSelectBin(unittest.TestCase):
    """Tests to check the correct working of the GetBinCutDf function"""
    def setUp(self):
        self.dfr = pd.DataFrame({
            'x': np.random.uniform(-1, 1, 500000)
        })


    def test_select_bin(self):
        """
        This mainly tests the correct handling of the used _get_var call in
        combination with parse_func_val
        """
        low, high = -0.5, 0.5
        # test getting variable by function
        sel_rows = mh.select_bin(lambda d: d.x, low, high)(self.dfr)

        pdt.assert_frame_equal(self.dfr[sel_rows],
                               self.dfr[(self.dfr.x > low) & (self.dfr.x < high)])

        # retrieval by string
        low, high = 0, 1
        sel_rows = mh.select_bin('x', low, high)(self.dfr)
        pdt.assert_frame_equal(self.dfr[sel_rows],
                               self.dfr[(self.dfr.x > low) & (self.dfr.x < high)])

        # retrieval with a parseable func val
        low, high = 0.25, 0.45
        sel_rows = mh.select_bin('abs(x)', low, high)(self.dfr)
        pdt.assert_frame_equal(self.dfr[sel_rows],
                               self.dfr[(self.dfr.x.abs() > low) &\
                                        (self.dfr.x.abs() < high)])


class TestFloatingPtRgx(unittest.TestCase):
    def test_float_rgx(self):
        match_cases = {
            s: float(s) for s in ['1.', '-1.', '2.34', '-2.23', '.34', '-.12',
                                  '0.12', '2', '-3', '+1.2', '+.3', '+2.', '+1']
        }
        flt_rgx = mh.float_rgx()
        for string, value in match_cases.iteritems():
            match = re.match(flt_rgx, string)
            self.assertTrue(match)
            self.assertAlmostEqual(value, float(match.group(1)))


    def test_float_rgx_no_match(self):
        nomatch_cases = ['.', '-.', '+.']
        flt_rgx = mh.float_rgx()
        for string in nomatch_cases:
            self.assertFalse(re.match(flt_rgx, string))


    def test_float_rgx_char_separated(self):
        match_cases = ['1p', '2p0', '-3p0', 'p1', '-p2', '+3p1', '-3p4', '0p30']
        flt_rgx = mh.float_rgx(True)
        for string in match_cases:
            match = re.match(flt_rgx, string)
            self.assertTrue(match)
            self.assertAlmostEqual(float(string.replace('p', '.')),
                                   float(match.group(1).replace('p', '.')))


    def test_float_rgx_char_separated_no_match(self):
        nomatch_cases = ['p', '-p', '+p']
        flt_rgx = mh.float_rgx(True)
        for string in nomatch_cases:
            self.assertFalse(re.match(flt_rgx, string))


class TestParseFuncVar(unittest.TestCase):
    def test_parse_func_var_valid(self):
        test_vars = (
            'JpsiPt', 'JpsiRap', 'chic_Mass', 'someFunkyVariable_34'
        )
        test_funcs = ('abs', 'sqrt', 'sin')

        for var in test_vars:
            self.assertEqual((var, None), mh.parse_func_var(var))

            for func in test_funcs:
                self.assertEqual(
                    (var, getattr(np, func)),
                    (mh.parse_func_var(func + '(' + var + ')'))
                )


    @patch('utils.misc_helpers.logging')
    def test_parse_func_var_invalid(self, mock_logger):
        res = mh.parse_func_var('some_funky_func(some_funky_var)')
        mock_logger.error.assert_called_with(
            'Could not find a numpy function named \'some_funky_func\' that '\
            'was parsed from \'some_funky_func(some_funky_var)\''
        )
        self.assertEqual(None, res)


class TestReplaceAll(unittest.TestCase):
    def test_replace_all(self):
        test_string = 'test abc 123'
        repl_pairs = [
            ('test', 'well this worked'),
            ('ab', 'AB'),
            ('abc', 'ABC'),
            ('1', '0')
        ]
        exp_string = 'well this worked ABC 023'
        self.assertEqual(exp_string, mh.replace_all(test_string, repl_pairs))

        test_string = 'What happens if we have two patterns with equal length?'
        repl_pairs = [
            ('app', 'APP'),
            ('ppe', 'I am not used'),
        ]
        exp_string = 'What hAPPens if we have two patterns with equal length?'
        self.assertEqual(exp_string, mh.replace_all(test_string, repl_pairs))

        # a slightly less obvious example where both patterns get replace once
        repl_pairs = [('pe', 'PE'), ('en', 'EN')]
        exp_string = 'What hapPEns if we have two patterns with equal lENgth?'
        self.assertEqual(exp_string, mh.replace_all(test_string, repl_pairs))


    def test_replace_all_reverse(self):
        test_string = 'test ABc foo bar 123'
        repl_pairs = [
            ('ba', 'AB'),
            ('abc', 'ABC'), # this is not used
            ('321', '123'),
            ('bar', 'foo'),
            ('foo', 'bar')
        ]
        exp_string = 'test bac bar foo 321'
        self.assertEqual(exp_string,
                         mh.replace_all(test_string, repl_pairs, reverse=True))

        repl_pairs = [
            ('I will be ignored', 'test'),
            ('But I am used', 'test '),
        ]
        exp_string = 'But I am usedABc foo bar 123'
        self.assertEqual(exp_string,
                         mh.replace_all(test_string, repl_pairs, reverse=True))


    def test_replace_all_reverse_inverts(self):
        test_string = 'test foo bar 123'
        repl_pairs = [
            ('test', 'whatever string we find'),
            ('foo', 'also should not matter'),
            ('bar', 'foo'), # even swapping should work
            ('123', '000')
        ]

        self.assertEqual(test_string,
                         mh.replace_all(mh.replace_all(test_string, repl_pairs),
                                        repl_pairs, reverse=True))


class TestParseSelExpr(unittest.TestCase):
    def setUp(self):
        self.dfr = pd.DataFrame({
            'x': np.random.uniform(-1, 1, 100000)
        })


    def test_parse_sel_expr_double_sel(self):
        """Test that a double sided expression can be parsed"""
        sel_expr = '-0.2 < x < 0.3'
        exp_sel = lambda d: (d.x > -0.2) & (d.x < 0.3)
        sel = mh.parse_sel_expr(sel_expr)
        pdt.assert_frame_equal(apply_selections(self.dfr, sel),
                               apply_selections(self.dfr, exp_sel))

        sel_expr = '0.1 < abs(x) < 0.2'
        exp_sel = lambda d: (d.x.abs() > 0.1) & (d.x.abs() < 0.2)
        sel = mh.parse_sel_expr(sel_expr)
        sel_dfr = apply_selections(self.dfr, sel)
        pdt.assert_frame_equal(apply_selections(self.dfr, sel),
                               apply_selections(self.dfr, exp_sel))


    def test_parse_sel_expr_left_sel(self):
        """Test that selections with only a left condition can be parsed"""
        sel_expr = '0.8 < x'
        exp_sel = lambda d: d.x > 0.8
        sel = mh.parse_sel_expr(sel_expr)
        pdt.assert_frame_equal(apply_selections(self.dfr, sel),
                               apply_selections(self.dfr, exp_sel))


    def test_parse_sel_expr_right_sel(self):
        """Test that selections with only a right condition can be parsed"""
        sel_expr = 'x < 0.2'
        exp_sel = lambda d: d.x < 0.2
        sel = mh.parse_sel_expr(sel_expr)
        pdt.assert_frame_equal(apply_selections(self.dfr, sel),
                               apply_selections(self.dfr, exp_sel))


    @patch('utils.misc_helpers.logging')
    def test_parse_sel_expr_fails(self, mock_logger):
        self.assertFalse(mh.parse_sel_expr('0.2 < x < 0.3 < 4'))
        self.assertFalse(mh.parse_sel_expr('ABC < foo'))
        self.assertFalse(mh.parse_sel_expr('1.2 < noValidFuncExpr(x)'))
        self.assertFalse(mh.parse_sel_expr('0.3 > x'))


class TestIsDivisable(unittest.TestCase):
    def test_is_divisable(self):
        self.assertEqual(mh.is_divisable(10, 5), 2)
        self.assertEqual(mh.is_divisable(60, 12), 5)

    @patch('utils.misc_helpers.logging')
    def test_is_not_divisable(self, mock_logger):
        self.assertEqual(mh.is_divisable(31, 3), None)
        mock_logger.warning.assert_called_with('3 is not a divisor of 31')

if __name__ == '__main__':
    unittest.main()
