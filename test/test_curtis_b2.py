import unittest
import os
import test_curtis
import chevalley as chv


class TestCurtisB2(test_curtis.TestCurtis):

    def setUp(self):
        os.chdir('..')
        t = chv.LieTypes.get_lie_type("b", 2)
        self.curtis = t.curtis

    def tearDown(self):
        os.chdir('test')


if __name__ == '__main__':
    unittest.main()
