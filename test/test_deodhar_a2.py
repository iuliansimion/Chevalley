import unittest
import numpy as np
import os
import test_deodhar
import chevalley as chv


class TestDeodharA2(test_deodhar.TestDeodhar):

    def setUp(self):
        os.chdir('..')
        t = chv.LieTypes.get_lie_type("a", 2)
        self.deodhar = t.deodhar
        self.ad = t.ad_action
        n = t.dim
        self.unu = np.ones((n, n), dtype=np.int)
        self.zero = np.zeros((n, n), dtype=np.int)

    def tearDown(self):
        os.chdir('test')


if __name__ == '__main__':
    unittest.main()
