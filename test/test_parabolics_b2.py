import unittest
import os
import test_parabolics
import chevalley as chv


class TestParabolicsB2(test_parabolics.TestParabolics):

    def setUp(self):
        os.chdir('..')
        t = chv.LieTypes.get_lie_type("b", 2)
        self.parabolics = t.parabolics

    def tearDown(self):
        os.chdir('test')


if __name__ == '__main__':
    unittest.main()
