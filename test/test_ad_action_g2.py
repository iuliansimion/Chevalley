import unittest
import os
import numpy as np
import test_ad_action
import chevalley as chv


class TestAdG2(test_ad_action.TestAd):

    def setUp(self):
        os.chdir('..')
        t = chv.LieTypes.get_lie_type("g", 2)
        self.ad = t.ad_action
        n = t.dim
        self.unu = np.ones((n, n), dtype=np.int)
        self.zero = np.zeros((n, n), dtype=np.int)

    def tearDown(self):
        os.chdir('test')


if __name__ == '__main__':
    unittest.main()
