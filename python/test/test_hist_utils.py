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


if __name__ == '__main__':
    unittest.main()
