import unittest
import os
import test_group
import chevalley as chv


class TestGroupA2(test_group.TestGroup):

    def setUp(self):
        os.chdir('..')
        t = chv.LieTypes.get_lie_type("a", 2)
        self.group = t.group
        self.variables = t.var

    def tearDown(self):
        os.chdir('test')
        del self.group
        del self.variables


if __name__ == '__main__':
    unittest.main()
