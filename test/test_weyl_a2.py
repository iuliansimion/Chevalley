import unittest
import os
import test_weyl
import chevalley as chv


class TestWeylA2(test_weyl.TestWeyl):

    def setUp(self):
        os.chdir('..')
        t = chv.LieTypes.get_lie_type("a", 2)
        self.rang = 2
        t.loadDistinguishedExpressions()
        self.weyl = t.weyl

    def tearDown(self):
        os.chdir('test')
        del self.weyl


if __name__ == '__main__':
    unittest.main()
