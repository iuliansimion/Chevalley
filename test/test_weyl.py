import unittest
import numpy as np


class TestWeyl(unittest.TestCase):
    weyl = None
    rang = 0

    """
    Test that the sums of two roots is calculated correctly
    """

    def test_perm_ord(self):
        perm = self.weyl.perm
        ords = self.weyl.perm_ord
        lung = self.weyl.nr_perm
        unu = perm[0]
        for i in range(lung):
            tmp = perm[i]
            o = 1
            while not np.array_equal(tmp, unu):
                tmp = tmp[perm[i]]
                o = o + 1
            self.assertEqual(o, ords[i])

    """
    Test that the card of plen_set equals the length of the word
    """

    def test_plen_set(self):
        lung = self.weyl.nr_perm
        pls = self.weyl.plen_set
        wlen = self.weyl.word_len
        for i in range(lung):
            self.assertEqual(len(pls(i)), wlen[i])

    """
    Test that the word expressions are reduced,
    --- that the length of the words is correct
    """

    def test_word_len(self):
        lung = self.weyl.nr_perm
        w = self.weyl.word
        wlen = self.weyl.word_len
        for i in range(lung):
            self.assertEqual(len(w(i)), wlen[i])

    """
    Test that the identity element is the firs one in list
    """

    def test_id_pos0(self):
        p0 = self.weyl.perm[0]
        nr = self.weyl.roots.nr_roots
        pid = [i for i in range(nr)]
        self.assertTrue(np.array_equal(p0, pid))

    """
    Test that reflections correspond to roots correctly
    --- we check that s_a(a)=-a
    --- and that the index of wprod(s_a) equals that of s_a
    """

    def test_refs_roots(self):
        refs = self.weyl.refs
        nr = self.weyl.roots.nr_roots
        npr = self.weyl.roots.nr_pos_roots
        p = self.weyl.perm
        wp = self.weyl.wprod
        for i in range(nr):
            if i < npr:
                self.assertEqual(p[refs[i]][i], i + npr)
            else:
                self.assertEqual(p[refs[i]][i], i - npr)
            self.assertEqual(wp([i]), refs[i])

    """
    Test that the words are in left descent order:
    i.e. that removing any number of terms on left,
    we get a word in self.weyl.words
    """

    def test_left_descent(self):
        nr = self.weyl.nr_perm
        word = self.weyl.word
        windex = self.weyl.windex
        for i in range(nr):
            w = word(i)
            for j in range(len(w)):
                ww = w[j:len(w)]
                self.assertTrue(windex(ww) > -1)

    """
    Test that the permutations are the right product of reflections in words
    """

    def test_perm_word(self):
        nrp = self.weyl.nr_perm
        word = self.weyl.word
        perm = self.weyl.perm
        refs = self.weyl.refs
        nr = self.weyl.roots.nr_roots
        self.assertTrue(np.array_equal(perm[0], [i for i in range(nr)]))
        for i in range(1, nrp):
            w = word(i)
            tmp = perm[refs[w[0]]]
            for s in w[1:len(w)]:
                tmp = tmp[perm[refs[s]]]
            self.assertTrue(np.array_equal(perm[i], tmp))

    """
    Test that distinguished expressions are
    --- 1. expressions
    --- 2. distinguished
    """

    # @unittest.skip("ttt")
    def test_dist_expr(self):
        self.weyl.load_dist_expr()
        de = self.weyl.dist_expr
        pprod = self.weyl.pprod
        lung = self.weyl.word_len
        word = self.weyl.word
        wprod = self.weyl.wprod
        refs = self.weyl.refs
        for c in de:
            x = c[0][0]
            xw = word(x)
            y = c[0][1]
            z = c[0][2]
            if x == 0:
                if y == z:
                    self.assertEqual(c[1], [[[], [[], [], []]]])
                else:
                    self.assertEqual(c[1], [])
            else:
                sz = pprod(refs[xw[0]], z)
                for e in c[1]:
                    # A
                    if 0 in e[1][0]:
                        self.assertTrue(lung[z] < lung[sz])
                    # B or C
                    if 0 in e[1][1] or 0 in e[1][2]:
                        self.assertTrue(lung[z] > lung[sz])
                    # check that expr*y=z
                    xx = [i for i in e[0] if i != -1]
                    # permutation corresponding to subexpression
                    xxperm = wprod(xx)
                    self.assertEqual(pprod(xxperm, y), z)
                    #
                    # distinguished condition
                    #
                    for j in range(len(e[0])):
                        if e[0][j] == -1:
                            xx = [i for i in e[0][j + 1:len(e[0])] if i != -1]
                            xxperm = wprod(xx)
                            xxpermy = pprod(xxperm, y)
                            sxxpermy = pprod(refs[xw[j]], xxpermy)
                            self.assertTrue(lung[sxxpermy] < lung[xxpermy])

    ###############################################
    #
    #
    # Testing for properties
    #
    #
    ###############################################
    """
    Check that for
    - fixed w
    - simple root a such that wa<0
    - root r<0 such that wr>0
    - if a+r is root then w(a+r)>0

    Doesn't hold: in B2

    """

    def check_simple_r_plus_phix(self, x):
        phix = self.weyl.plen_set(x)
        rsum = self.weyl.roots.sum
        perm = self.weyl.perm
        npr = self.weyl.roots.nr_pos_roots
        for s in range(self.rang):
            for r in phix:
                minus = perm[x][r]
                ss = rsum[s][minus]
                if ss != -1 and s != r:
                    tmp = perm[x][ss]
                    self.assertTrue(tmp < npr)

    @unittest.skip("Property doesn't hold: in B2")
    def test_simple_r_plus_phix(self):
        for i in range(self.weyl.nr_perm):
            self.check_simple_r_plus_phix(i)
