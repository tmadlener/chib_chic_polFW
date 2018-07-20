#!/usr/bin/env python
"""
Tests for symbolic functions
"""
import unittest
import sympy as sp

from utils.symbolic import lth_2, lth_1

class TestLambdaThetas(unittest.TestCase):
    def test_lth_2(self):
        """Test if lth_2 returns the correct values for some exemplary inputs"""
        # unpolarized case
        self.assertEqual(lth_2(R1=sp.S(2)/5, R2=sp.S(2)/5), 0)
        # maximally negative
        self.assertEqual(lth_2(R1=0, R2=0), -sp.S(3)/5)
        # maximally positive
        self.assertEqual(lth_2(R1=0, R2=1), 1)
        # |Jz|=1 pure state has lambd_theta = -1/3
        self.assertEqual(lth_2(R1=1, R2=0), -sp.S(1)/3)


    def test_lth_1(self):
        """Test if lth_1 returns the correct values for some exemplary inputs"""
        # unpolarized case
        self.assertEqual(lth_1(R=sp.S(2)/3), 0)
        # maximally negative
        self.assertEqual(lth_1(R=1), -sp.S(1)/3)
        # maximally positive
        self.assertEqual(lth_1(R=0), 1)


if __name__ == '__main__':
    unittest.main()
