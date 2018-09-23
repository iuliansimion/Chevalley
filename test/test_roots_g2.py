import unittest
import os
import test_roots
import chevalley as chv


class TestRootsG2(test_roots.TestRoots):

    def setUp(self):
        os.chdir('..')
        self.rootsystem = chv.LieTypes.get_lie_type("g", 2).rootsystem

    def tearDown(self):
        os.chdir('test')


if __name__ == '__main__':
    unittest.main()
