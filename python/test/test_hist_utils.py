#!/usr/bin/env python
"""
Test the hist utils package
"""

import unittest
from mock import patch

import numpy as np
import numpy.testing as npt

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from root_numpy import fill_hist

import utils.hist_utils as hu

from utils.misc_helpers import create_random_str

def _get_hist(ndim):
    """Create a histogram of a given dimensionality"""
    if ndim == 3:
        return hu.create_histogram(np.random.uniform(0, 1, (10000, 3)),
                                   (10, 0, 1, 20, 0, 1, 30, 0, 1))
    if ndim == 2:
        return hu.create_histogram(np.random.uniform(0, 1, (10000, 2)),
                                   (10, 0, 1, 20, 0, 1))
    if ndim == 1:
        return hu.create_histogram(np.random.uniform(0, 1, 10000),
                                   (10, 0, 1))
    if ndim >= 4:
        hist = r.THnD('hist4d', '', ndim, np.arange(2, ndim + 2, 1, dtype='i4'),
                      np.zeros(ndim), np.ones(ndim))
        fill_vals = np.random.uniform(0, 1, (10000, ndim))
        for i in xrange(len(fill_vals)):
            hist.Fill(fill_vals[i, :])
        hist.Sumw2()
        return hist



class TestGetArray(unittest.TestCase):
    def _comp_array_hist(self, array, hist, overflow=False, error=False):
        cfunc = 'GetBinError' if error else 'GetBinContent'
        ndim = len(array.shape)

        # indexing offsets for TH objects depends on whether overflow is
        # included or not
        start = 0 if overflow else 1 # Since 0 is underflow bin
        end = 1 if overflow else 0 # Since Nbins + 1 is overflow bin
        off = 0 if overflow else -1 # since "real" bins start at 1 in TH1

        if ndim == 1:
            for i in xrange(start, hist.GetNbinsX() + end + 1):
                self.assertAlmostEqual(array[i + off], getattr(hist, cfunc)(i))

        if ndim == 2:
            for i in xrange(start, hist.GetNbinsX() + end + 1):
                for j in xrange(start, hist.GetNbinsY() + end + 1):
                    self.assertAlmostEqual(array[i + off, j + off],
                                           getattr(hist, cfunc)(i, j))

        if ndim == 3:
            for i in xrange(start, hist.GetNbinsX() + end + 1):
                for j in xrange(start, hist.GetNbinsY() + end + 1):
                    for k in xrange(start, hist.GetNbinsZ() + end + 1):
                        self.assertAlmostEqual(array[i + off, j + off, k + off],
                                               getattr(hist, cfunc)(i, j, k))

        if ndim == 4:
            for i in xrange(start, hist.GetAxis(0).GetNbins() + end + 1):
                for j in xrange(start, hist.GetAxis(1).GetNbins() + end + 1):
                    for k in xrange(start, hist.GetAxis(2).GetNbins() + end + 1):
                        for l in xrange(start, hist.GetAxis(3).GetNbins() + end + 1):
                            self.assertAlmostEqual(array[i+off, j+off, k+off, l+off],
                                                   getattr(hist, cfunc)(np.array([i,j,k,l], dtype='i4')))


    def test_get_array_1dim(self):
        hist = _get_hist(1)
        arr = hu.get_array(hist)
        self.assertEqual(arr.shape[0], hist.GetNbinsX())
        self._comp_array_hist(arr, hist)


    def test_get_array_2dim(self):
        hist = _get_hist(2)
        arr = hu.get_array(hist)
        self.assertEqual(arr.shape[0], hist.GetNbinsX())
        self.assertEqual(arr.shape[1], hist.GetNbinsY())
        self._comp_array_hist(arr, hist)


    def test_get_array_3dim(self):
        hist = _get_hist(3)
        arr = hu.get_array(hist)
        self.assertEqual(arr.shape[0], hist.GetNbinsX())
        self.assertEqual(arr.shape[1], hist.GetNbinsY())
        self.assertEqual(arr.shape[2], hist.GetNbinsZ())
        self._comp_array_hist(arr, hist)


    def test_get_array_4dim(self):
        """Test if THns also work"""
        hist = _get_hist(4)

        arr = hu.get_array(hist)
        self.assertEqual(arr.shape[0], hist.GetAxis(0).GetNbins())
        self.assertEqual(arr.shape[1], hist.GetAxis(1).GetNbins())
        self.assertEqual(arr.shape[2], hist.GetAxis(2).GetNbins())
        self.assertEqual(arr.shape[3], hist.GetAxis(3).GetNbins())
        self._comp_array_hist(arr, hist)


    def test_get_array_w_overflow(self):
        hist = _get_hist(1)
        arr = hu.get_array(hist, overflow=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX() + 2)
        self._comp_array_hist(arr, hist, True)


        hist = _get_hist(2)
        arr = hu.get_array(hist, overflow=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX() + 2)
        self.assertEqual(arr.shape[1], hist.GetNbinsY() + 2)
        self._comp_array_hist(arr, hist, True)


        hist = _get_hist(3)
        arr = hu.get_array(hist, overflow=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX() + 2)
        self.assertEqual(arr.shape[1], hist.GetNbinsY() + 2)
        self.assertEqual(arr.shape[2], hist.GetNbinsZ() + 2)
        self._comp_array_hist(arr, hist, True)

        hist = _get_hist(4)
        arr = hu.get_array(hist, overflow=True)
        self.assertEqual(arr.shape[0], hist.GetAxis(0).GetNbins() + 2)
        self.assertEqual(arr.shape[1], hist.GetAxis(1).GetNbins() + 2)
        self.assertEqual(arr.shape[2], hist.GetAxis(2).GetNbins() + 2)
        self.assertEqual(arr.shape[3], hist.GetAxis(3).GetNbins() + 2)
        self._comp_array_hist(arr, hist, True)


    def test_get_array_w_errors(self):
        hist = _get_hist(1)
        arr = hu.get_array(hist, overflow=True, errors=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX() + 2)
        self._comp_array_hist(arr, hist, True, True)


        hist = _get_hist(2)
        arr = hu.get_array(hist, overflow=True, errors=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX() + 2)
        self.assertEqual(arr.shape[1], hist.GetNbinsY() + 2)
        self._comp_array_hist(arr, hist, True, True)


        hist = _get_hist(3)
        arr = hu.get_array(hist, overflow=True, errors=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX() + 2)
        self.assertEqual(arr.shape[1], hist.GetNbinsY() + 2)
        self.assertEqual(arr.shape[2], hist.GetNbinsZ() + 2)
        self._comp_array_hist(arr, hist, True, True)

        hist = _get_hist(4)
        arr = hu.get_array(hist, overflow=True, errors=True)
        self.assertEqual(arr.shape[0], hist.GetAxis(0).GetNbins() + 2)
        self.assertEqual(arr.shape[1], hist.GetAxis(1).GetNbins() + 2)
        self.assertEqual(arr.shape[2], hist.GetAxis(2).GetNbins() + 2)
        self.assertEqual(arr.shape[3], hist.GetAxis(3).GetNbins() + 2)
        self._comp_array_hist(arr, hist, True, True)



    def test_get_array_w_errors_no_overflow(self):
        hist = _get_hist(1)
        arr = hu.get_array(hist, errors=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX())
        self._comp_array_hist(arr, hist, error=True)

        hist = _get_hist(2)
        arr = hu.get_array(hist, errors=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX())
        self.assertEqual(arr.shape[1], hist.GetNbinsY())
        self._comp_array_hist(arr, hist, error=True)

        hist = _get_hist(3)
        arr = hu.get_array(hist, errors=True)
        self.assertEqual(arr.shape[0], hist.GetNbinsX())
        self.assertEqual(arr.shape[1], hist.GetNbinsY())
        self.assertEqual(arr.shape[2], hist.GetNbinsZ())
        self._comp_array_hist(arr, hist, error=True)

        hist = _get_hist(4)
        arr = hu.get_array(hist, errors=True)
        self.assertEqual(arr.shape[0], hist.GetAxis(0).GetNbins())
        self.assertEqual(arr.shape[1], hist.GetAxis(1).GetNbins())
        self.assertEqual(arr.shape[2], hist.GetAxis(2).GetNbins())
        self.assertEqual(arr.shape[3], hist.GetAxis(3).GetNbins())
        self._comp_array_hist(arr, hist, error=True)



class TestFindBin(unittest.TestCase):
    def test_find_bin_reg_binning(self):
        hist = r.TH1D(create_random_str(8), '', 10, 0, 1)
        binning = hu.get_binning(hist)

        values = np.random.uniform(0, 1, 1000)
        exp_idcs = np.array([hist.FindBin(v) for v in values])
        exp_idcs -= 1 # correct for TH1 indexing starting at 1

        bin_idcs = hu.find_bin(binning, values)
        npt.assert_equal(bin_idcs, exp_idcs)


    def test_find_bin_nonreg_binning(self):
        hist = r.TH1D(create_random_str(8), '', 10, np.linspace(0, 1, 11)**2)
        binning = hu.get_binning(hist)

        values = np.random.uniform(0, 1, 1000)
        exp_idcs = np.array([hist.FindBin(v) for v in values])
        exp_idcs -= 1

        bin_idcs = hu.find_bin(binning, values)
        npt.assert_equal(bin_idcs, exp_idcs)


    @patch('utils.hist_utils.logging')
    def test_find_bin_warning(self, mock_logger):
        exp_warn = 'When trying to find the bin indices at least one value '\
                   'could not be attributed to a bin in the passed binning'
        bins = hu.get_binning(hist = r.TH1D(create_random_str(), '', 10, 0, 1))
        values = np.array([-0.1, 0.2, 0.3, 0.4])
        bin_idcs = hu.find_bin(bins, values)
        mock_logger.warn.assert_called_with(exp_warn)

        values = np.array([0.1, 0.2, 1.3, 0.4, 0.5])
        bin_idcs = hu.find_bin(bins, values)
        mock_logger.warn.assert_called_with(exp_warn)


class TestFromArray(unittest.TestCase):
    def _test_from_array_nd(self, n_dim):
        hist = _get_hist(n_dim)
        arr = hu.get_array(hist)
        axes = 'X'
        if n_dim == 2:
            axes = 'XY'
        if n_dim == 3:
            axes = 'XYZ'
        binning = np.array([hu.get_binning(hist, ax) for ax in axes])

        arr_hist = hu.from_array(arr, binning)

        npt.assert_equal(hu.get_array(arr_hist), arr)
        npt.assert_equal(hu.get_binning(arr_hist, 'X'), hu.get_binning(hist, 'X'))
        if n_dim > 1:
            npt.assert_equal(hu.get_binning(arr_hist, 'Y'), hu.get_binning(hist, 'Y'))
        if n_dim > 2:
            npt.assert_equal(hu.get_binning(arr_hist, 'Z'), hu.get_binning(hist, 'Z'))

        err = hu.get_array(hist, errors=True)
        arr_err_hist = hu.from_array(arr, binning, errors=err)
        npt.assert_equal(hu.get_array(arr_err_hist), arr)
        npt.assert_equal(hu.get_array(arr_err_hist, errors=True), err)


    def _test_from_array_nd_w_overflow(self, n_dim):
        hist = _get_hist(n_dim)
        arr = hu.get_array(hist, overflow=True)
        axes = 'X'
        if n_dim == 2:
            axes = 'XY'
        if n_dim == 3:
            axes = 'XYZ'
        binning = np.array([hu.get_binning(hist, ax) for ax in axes])

        arr_hist = hu.from_array(arr, binning)

        npt.assert_equal(hu.get_array(arr_hist, overflow=True), arr)
        npt.assert_equal(hu.get_binning(arr_hist, 'X'), hu.get_binning(hist, 'X'))
        if n_dim > 1:
            npt.assert_equal(hu.get_binning(arr_hist, 'Y'), hu.get_binning(hist, 'Y'))
        if n_dim > 2:
            npt.assert_equal(hu.get_binning(arr_hist, 'Z'), hu.get_binning(hist, 'Z'))

        err = hu.get_array(hist, errors=True, overflow=True)
        arr_err_hist = hu.from_array(arr, binning, errors=err)
        npt.assert_equal(hu.get_array(arr_err_hist, overflow=True), arr)
        npt.assert_equal(hu.get_array(arr_err_hist, overflow=True, errors=True), err)


    def test_from_array(self):
        self._test_from_array_nd(1)
        self._test_from_array_nd(2)
        self._test_from_array_nd(3)


    def test_from_array_w_overflow(self):
        self._test_from_array_nd_w_overflow(1)
        self._test_from_array_nd_w_overflow(2)
        self._test_from_array_nd_w_overflow(3)


    def test_raises_errors(self):
        # binning, array mismatch
        with self.assertRaises(ValueError):
            hu.from_array(np.random.uniform(0, 1, 10), np.linspace(0, 1, 30))

        with self.assertRaises(ValueError):
            hu.from_array(np.random.uniform(0, 1, (10, 10)),
                          np.array([np.linspace(0, 1, 11), np.linspace(0, 1, 12)]))

        with self.assertRaises(ValueError):
            hu.from_array(np.random.uniform(0, 1, (10, 10, 10)),
                          np.array([np.linspace(0, 1, 11), np.linspace(0, 1, 11),
                                    np.linspace(0, 1, 12)]))

        with self.assertRaises(ValueError):
            hu.from_array(np.random.uniform(0, 1, (10, 2)), np.linspace(0, 1, 30))

        # wrong error dimension
        with self.assertRaises(ValueError):
            hu.from_array(np.random.uniform(0, 1, 10), np.linspace(0, 1, 11),
                          errors=np.random.uniform(0, 1, 11))

        with self.assertRaises(ValueError):
            hu.from_array(np.random.uniform(0, 1, (10, 10)),
                          np.array([np.linspace(0, 1, 11), np.linspace(0, 1, 11)]),
                          errors=np.random.uniform(0, 1, (11, 2)))

        # unhandleable array dimension
        with self.assertRaises(ValueError):
            hu.from_array(np.random.uniform(0, 1, (10, 4, 4, 5)), np.random.uniform(0, 1, 10))


    def test_handles_binning(self):
        # Test if binning arrays are handled correctly (type conversion for
        # ROOT to understand)
        hist = hu.from_array(np.random.uniform(0, 1, 10), np.arange(0, 11, 1))
        npt.assert_equal(hu.get_binning(hist, 'X'), np.arange(0, 11, 1))


class TestProject(unittest.TestCase):
    def test_project_3d_to_2d(self):
        hist_3d = _get_hist(3)
        # populate overflow bins to make sure that they are treated as expected
        fill_hist(hist_3d, np.random.uniform(-1, 0, (100, 3)))
        fill_hist(hist_3d, np.random.uniform(1, 2, (100, 3)))
        val3d, err3d = hu.get_array(hist_3d), hu.get_array(hist_3d, errors=True)

        hist_xy = hu.project(hist_3d, 'xy')
        val, err = hu.get_array(hist_xy), hu.get_array(hist_xy, errors=True)
        npt.assert_equal(val, np.sum(val3d, axis=2))
        npt.assert_equal(err, np.sqrt(np.sum(err3d**2, axis=2)))

        hist_yz = hu.project(hist_3d, 'yz')
        val, err = hu.get_array(hist_yz), hu.get_array(hist_yz, errors=True)
        npt.assert_equal(val, np.sum(val3d, axis=0))
        npt.assert_equal(err, np.sqrt(np.sum(err3d**2, axis=0)))

        hist_zx = hu.project(hist_3d, 'zx')
        val, err = hu.get_array(hist_zx), hu.get_array(hist_zx, errors=True)
        npt.assert_equal(val, np.sum(val3d, axis=1).T)
        npt.assert_equal(err, np.sqrt(np.sum(err3d**2, axis=1)).T)

        hist_yx = hu.project(hist_3d, 'yx')
        val, err = hu.get_array(hist_yx), hu.get_array(hist_yx, errors=True)
        npt.assert_equal(val, np.sum(val3d, axis=2).T)
        npt.assert_equal(err, np.sqrt(np.sum(err3d**2, axis=2)).T)


    def test_project_3d_to_1d(self):
        hist3d = _get_hist(3)
        # populate overflow bins to make sure that they are treated as expected
        fill_hist(hist3d, np.random.uniform(-1, 0, (100, 3)))
        fill_hist(hist3d, np.random.uniform(1, 2, (100, 3)))
        val3d, err3d = hu.get_array(hist3d), hu.get_array(hist3d, errors=True)

        hist_x = hu.project(hist3d, 'x')
        val, err = hu.get_array(hist_x), hu.get_array(hist_x, errors=True)
        npt.assert_equal(val, np.sum(val3d, axis=(1,2)))
        npt.assert_equal(err, np.sqrt(np.sum(err3d**2, axis=(1,2))))

        hist_y = hu.project(hist3d, 'y')
        val, err = hu.get_array(hist_y), hu.get_array(hist_y, errors=True)
        npt.assert_equal(val, np.sum(val3d, axis=(0,2)))
        npt.assert_equal(err, np.sqrt(np.sum(err3d**2, axis=(0,2))))

        hist_z = hu.project(hist3d, 'z')
        val, err = hu.get_array(hist_z), hu.get_array(hist_z, errors=True)
        npt.assert_equal(val, np.sum(val3d, axis=(0,1)))
        npt.assert_equal(err, np.sqrt(np.sum(err3d**2, axis=(0,1))))


    def test_project_2d_to_1d(self):
        hist2d = _get_hist(2)
        # populate overflow bins to make sure that they are treated as expected
        fill_hist(hist2d, np.random.uniform(-1, 0, (100, 2)))
        fill_hist(hist2d, np.random.uniform(1, 2, (100, 2)))
        val2d, err2d = hu.get_array(hist2d), hu.get_array(hist2d, errors=True)

        hist_x = hu.project(hist2d, 'x')
        val, err = hu.get_array(hist_x), hu.get_array(hist_x, errors=True)
        npt.assert_equal(val, np.sum(val2d, axis=1))
        npt.assert_equal(err, np.sqrt(np.sum(err2d**2, axis=1)))


        hist_y = hu.project(hist2d, 'y')
        val, err = hu.get_array(hist_y), hu.get_array(hist_y, errors=True)
        npt.assert_equal(val, np.sum(val2d, axis=0))
        npt.assert_equal(err, np.sqrt(np.sum(err2d**2, axis=0)))


    def test_project_4d_to_3d(self):
        # Fill a test histogram
        hist4d = _get_hist(4)
        val4d, err4d = hu.get_array(hist4d), hu.get_array(hist4d, errors=True)

        hist_123 = hu.project(hist4d, [1, 2, 3])
        self.assertTrue(isinstance(hist_123, r.TH3))
        val, err = hu.get_array(hist_123), hu.get_array(hist_123, errors=True)
        npt.assert_allclose(val, np.sum(val4d, axis=0))
        npt.assert_allclose(err, np.sqrt(np.sum(err4d**2, axis=0)))

        hist_021 = hu.project(hist4d, [0, 2, 1])
        self.assertTrue(isinstance(hist_021, r.TH3))
        val, err = hu.get_array(hist_021), hu.get_array(hist_021, errors=True)
        npt.assert_allclose(val, np.swapaxes(np.sum(val4d, axis=3), 1, 2))
        npt.assert_allclose(err, np.swapaxes(np.sqrt(np.sum(err4d**2, axis=3)), 1, 2))


    def test_project_4d_to_2d(self):
        hist4d = _get_hist(4)
        val4d, err4d = hu.get_array(hist4d), hu.get_array(hist4d, errors=True)

        hist_01 = hu.project(hist4d, [0, 1])
        self.assertTrue(isinstance(hist_01, r.TH2))
        val, err = hu.get_array(hist_01), hu.get_array(hist_01, errors=True)
        # For some reason the sum does not give the expected shape here,
        # However transposing helps
        # But the projected histogram has the expected shape and contens
        # TODO: investigate why
        npt.assert_allclose(val, np.sum(val4d, axis=(2, 3)).T)
        npt.assert_allclose(err, np.sqrt(np.sum(err4d**2, axis=(2, 3))).T)


    def test_project_4d_to_1d(self):
        hist4d = _get_hist(4)
        val4d, err4d = hu.get_array(hist4d), hu.get_array(hist4d, errors=True)

        hist_0 = hu.project(hist4d, 0)
        self.assertTrue(isinstance(hist_0, r.TH1))
        val, err = hu.get_array(hist_0), hu.get_array(hist_0, errors=True)
        npt.assert_allclose(val, np.sum(val4d, axis=(1, 2, 3)))
        npt.assert_allclose(err, np.sqrt(np.sum(err4d**2, axis=(1,2,3))))


    def test_project_5d_to_4d(self):
        """Test that also returning into a THn works"""
        hist5d = _get_hist(5)
        val5d, err5d = hu.get_array(hist5d), hu.get_array(hist5d, errors=True)

        hist_0234 = hu.project(hist5d, [0, 2, 3, 4])
        val, err = hu.get_array(hist_0234), hu.get_array(hist_0234, errors=True)
        npt.assert_allclose(val, np.sum(val5d, axis=1))
        npt.assert_allclose(err, np.sqrt(np.sum(err5d**2, axis=1)))

        hist_0314 = hu.project(hist5d, [0, 3, 1, 4])
        val, err = hu.get_array(hist_0314), hu.get_array(hist_0314, errors=True)
        npt.assert_allclose(val, np.swapaxes(np.sum(val5d, axis=2), 1, 2))
        npt.assert_allclose(err, np.swapaxes(np.sqrt(np.sum(err5d**2, axis=2)), 1, 2))


class TestUncerHist(unittest.TestCase):
    def test_abs_uncer(self):
        for i in xrange(1, 4):
            hist = _get_hist(i)
            uncer_hist = hu.uncer_hist(hist, abs_uncer=True)
            npt.assert_allclose(np.sqrt(hu.get_array(hist)),
                                hu.get_array(uncer_hist))


    def test_rel_uncer(self):
        for i in xrange(1, 4):
            hist = _get_hist(i)
            uncer_hist = hu.uncer_hist(hist)
            vals = hu.get_array(hist)
            exp_uncers = np.zeros_like(vals)
            np.divide(1, np.sqrt(vals), out=exp_uncers, where=vals!=0)
            npt.assert_allclose(exp_uncers, hu.get_array(uncer_hist))


class TestRebin(unittest.TestCase):
    """
    Tests for checking rebinning of histograms. Since all the heavy lifting
    in these cases is done by ROOT (and hopefully tested) there, only some basic
    checks are done here:

    - Input histogram remains unchanged
    - The number of bins corresponds to the ones that are expected
    - Errors are caught and reported correctly
    """
    def test_rebin_4d(self):
        hist = _get_hist(4)
        # NOTE: only very limited possibilities to do an actual rebinning
        # here
        re_hist = hu.rebin(hist, [(0, 1), (2, 2)])
        # check original hist is unchanged
        self.assertEqual(hu._get_nbins(hist), (2, 3, 4, 5))
        self.assertEqual(hu._get_nbins(re_hist), (1, 3, 2, 5))


    def test_rebin_3d(self):
        hist = _get_hist(3)
        # some more possibilities here
        re_hist = hu.rebin(hist, [(0, 2), (1, 5), (2, 6)])

        self.assertEqual(hu._get_nbins(hist), (10, 20, 30))
        self.assertEqual(hu._get_nbins(re_hist), (2, 5, 6))

        # check that char indices work
        re_hist = hu.rebin(hist, [('x', 5), ('Y', 2), ('z', 3)])
        self.assertEqual(hu._get_nbins(hist), (10, 20, 30))
        self.assertEqual(hu._get_nbins(re_hist), (5, 2, 3))

        # check that in principle also mixed indices work
        re_hist = hu.rebin(hist, [('x', 5), ('Y', 2), (2, 3)])
        self.assertEqual(hu._get_nbins(hist), (10, 20, 30))
        self.assertEqual(hu._get_nbins(re_hist), (5, 2, 3))


    def test_rebin_2d(self):
        hist = _get_hist(2)
        # only checking the basic case here, since basically only the return
        # statement is different from the 3d case
        re_hist = hu.rebin(hist, [(0, 2), (1, 4)])
        self.assertEqual(hu._get_nbins(hist), (10, 20))
        self.assertEqual(hu._get_nbins(re_hist), (2, 4))


    def test_rebin_1d(self):
        hist = _get_hist(1)
        # only checking the basic case here, since basically only the return
        # statement is different from the 3d case
        re_hist = hu.rebin(hist, [('x', 2)])
        self.assertEqual(hu._get_nbins(hist), (10,))
        self.assertEqual(hu._get_nbins(re_hist), (2,))


class TestRebin1DBinning(unittest.TestCase):
    """
    Tests checking that binning into a non-uniform binning works
    """
    @patch('utils.hist_utils.logging')
    def test_non_compatible_binning(self, mock_logger):
        hist = _get_hist(1)
        non_comp_binning = np.linspace(0, 1, 7)

        exp_err = 'Cannot rebin histogram with binning {} to target binning {}'
        self.assertTrue(hu.rebin_1d_binning(hist, non_comp_binning) is None)
        mock_logger.error.assert_called_with(exp_err.format(hu.get_binning(hist),
                                                            non_comp_binning))

    def test_uniform_rebinning(self):
        hist = _get_hist(1)
        targ_bin = np.linspace(0, 1, 6)

        # since hist has 10 bin, rebin(hist, 5) should return the same histogram
        # as the one when we rebin it according to a custom binning
        npt.assert_allclose(hu.get_array(hu.rebin_1d_binning(hist, targ_bin)),
                            hu.get_array(hu.rebin(hist, [(0, 5)])))


    def test_non_uniform_binning(self):
        hist = _get_hist(1)
        targ_bin = np.array([0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])

        # manually rebin this histogram here, since we now what we want
        vals = hu.get_array(hist)
        targ_vals = np.array([vals[0], vals[1], vals[2] + vals[3],
                              vals[4] + vals[5], vals[6] + vals[7],
                              vals[8] + vals[9]])

        npt.assert_allclose(hu.get_array(hu.rebin_1d_binning(hist, targ_bin)),
                            targ_vals)

    def test_non_full_coverage(self):
        """
        Check that if the target binning does not cover the full range of the
        original histogram things still work
        """
        hist = _get_hist(1)
        targ_bin = np.linspace(0.1, 0.6, 6)
        targ_vals = hu.get_array(hist)[1:6]

        npt.assert_allclose(hu.get_array(hu.rebin_1d_binning(hist, targ_bin)),
                            targ_vals)


if __name__ == '__main__':
    unittest.main()
