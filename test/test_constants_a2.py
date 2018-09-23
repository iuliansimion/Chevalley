import unittest
import os
import test_constants
import chevalley as chv


class TestConstantsA2(test_constants.TestConstants):

    def setUp(self):
        os.chdir('..')
        t = chv.LieTypes.get_lie_type("a", 2)
        self.roots = t.rootsystem
        self.constants = t.constants

    def tearDown(self):
        os.chdir('test')


if __name__ == '__main__':
    unittest.main()
