import unittest
import numpy as np


class TestRoots(unittest.TestCase):
    rootsystem = None

    """
    Test that the sums of two roots is calculated correctly
    """

    # @unittest.skip("test")
    def test_roots_sum(self):
        rr = self.rootsystem
        for i in range(rr.nr_roots):
            for j in range(rr.nr_roots):
                tmp = rr.index(rr.roots[i] + rr.roots[j])
                self.assertEqual(tmp, rr.rsum[i][j])

    """
    Test that the short roots are determined correctly
    """

    # @unittest.skip("test")
    def test_short_roots(self):
        rr = self.rootsystem
        sr = rr.determine_short_roots()
        for i in range(rr.nr_roots):
            self.assertEqual(rr.short[i], sr[i])

    """
    Test that -roots[i] is at position i+nr_pos_roots
    for all positive roots root[i]
    """

    def test_pos_neg_roots(self):
        rr = self.rootsystem
        for i in range(rr.nr_pos_roots):
            self.assertTrue(np.array_equal(rr.roots[i], -rr.roots[i + rr.nr_pos_roots]))
