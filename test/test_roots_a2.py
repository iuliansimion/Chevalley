import unittest
import os
import test_roots
import chevalley as chv


class TestRootsA2(test_roots.TestRoots):

    def setUp(self):
        os.chdir('..')
        self.rootsystem = chv.LieTypes.get_lie_type("a", 2).rootsystem

    def tearDown(self):
        os.chdir('test')


if __name__ == '__main__':
    unittest.main()
