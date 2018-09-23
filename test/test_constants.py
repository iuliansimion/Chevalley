import unittest


class TestConstants(unittest.TestCase):
    roots = None
    constants = None

    """
    Test that on extra special pairs the sign of Nrs is +1
    """

    def test_extra_special(self):
        esp = self.constants.roots.extra_special_pairs()
        N = self.constants.N
        for e in esp:
            self.assertTrue(N[e[0]][e[1]] > 0)

    def test_absNrs(self):
        tmp = self.constants.absNrs()
        N = self.constants.N
        nr = self.constants.roots.nr_roots
        for r in range(nr):
            for s in range(nr):
                if s != r:
                    self.assertEqual(abs(N[r][s]), tmp[r][s])

    """
    Check Carter p.55 Theorem 4.1.2 (i)
    """

    def check_fill_i(self, S):
        nr = self.constants.roots.nr_roots
        for r in range(nr):
            for s in range(nr):
                if r != s and S[r][s] != -S[s][r] != 0:
                    return False
        return True

    """
    Check Carter p.55 Theorem 4.1.2 (ii)
    """

    def check_fill_ii(self, S):
        rsum = self.constants.roots.rsum
        nr = self.constants.roots.nr_roots
        pnr = self.constants.roots.nr_pos_roots
        for r1 in range(nr):
            for r2 in range(nr):
                if rsum[r1][r2] != -1:
                    # determine -(r1+r2)
                    if rsum[r1][r2] >= pnr:
                        r3 = rsum[r1][r2] - pnr
                    else:
                        r3 = rsum[r1][r2] + pnr
                    # now r1+r2+r3=0
                    self.assertEqual(S[r1][r2], S[r2][r3])
                    self.assertEqual(S[r2][r3], S[r3][r1])
        return True

    """
    Check Carter p.55 Theorem 4.1.2 (iii)
    """

    def check_fill_iii(self, S):
        pnr = self.constants.roots.nr_pos_roots
        for r in range(pnr):
            for s in range(pnr):
                self.assertEqual(S[r][s], -S[r + pnr][s + pnr])
                self.assertEqual(S[r][s + pnr], -S[r + pnr][s])
                self.assertEqual(S[r + pnr][s], -S[r][s + pnr])
        return True

    """
    Check Carter p.55 Theorem 4.1.2 (iv)
    """

    def check_fill_iv(self, S, absN):
        rsum = self.constants.roots.rsum
        nr = self.constants.roots.nr_roots
        index = self.constants.roots.index
        rr = self.constants.roots.roots
        lss = self.constants.roots.ls_square
        short = self.constants.roots.short
        npo = self.constants.roots.no_pair_opp
        pdi = self.constants.roots.pair_distinct
        for r1 in range(nr):
            for r2 in range(nr):
                for r3 in range(nr):
                    r4 = index(-(rr[r1] + rr[r2] + rr[r3]))
                    if r4 != -1 and npo(r1, r2, r3, r3) and pdi(r1, r2, r3, r4):
                        # now r1+r2+r3+r4=0
                        if not short[rsum[r1][r2]]:
                            num1 = lss
                        else:
                            num1 = 1
                        if not short[rsum[r2][r3]]:
                            num2 = lss
                        else:
                            num2 = 1
                        if not short[rsum[r3][r1]]:
                            num3 = lss
                        else:
                            num3 = 1

                        print(r1, r2, r3, r4)
                        s = 0
                        s = s + S[r1][r2] * absN[r1][r2] * S[r3][r4] * absN[r3][r4] / num1
                        s = s + S[r2][r3] * absN[r2][r3] * S[r1][r4] * absN[r1][r4] / num2
                        s = s + S[r3][r1] * absN[r3][r1] * S[r2][r4] * absN[r2][r4] / num3
                        self.assertEqual(s, 0)

        return True

    """
    Check that the signs corresponding to non-zero Nrs are non-zero
    """

    def check_sgn_complete(self, S, absN):
        nr = self.constants.roots.nr_roots
        for r in range(nr):
            for s in range(nr):
                if r != s:
                    if absN[r][s] != 0:
                        self.assertTrue(S[r][s] != 0)

    """
    Test for calculating with Carter p.55 Theorem 4.1.2
    """

    def test_Nrs(self):
        absN = self.constants.absNrs()
        S = self.constants.sgnNrs()
        self.check_sgn_complete(S, absN)
        self.assertTrue(self.check_fill_i(S))
        self.assertTrue(self.check_fill_ii(S))
        self.assertTrue(self.check_fill_iii(S))
        self.assertTrue(self.check_fill_iv(S, absN))

    """
    Test for calculation of Mrsi
    - Carter p.61
    - in particular test that the calculated constants are integers
    """

    def test_Mrsi(self):
        M1 = self.constants.M
        M2 = self.constants.Mrsi()
        nr = self.constants.roots.nr_roots
        index = self.constants.roots.index
        rr = self.constants.roots.roots
        for r in range(nr):
            for s in range(nr):
                i = 1
                tmp1 = 1
                tmp2 = 1
                while index((i - 1) * rr[r] + rr[s]) != -1:
                    v = index((i - 1) * rr[r] + rr[s])
                    if self.constants.N[r][v] == 0:
                        break
                    tmp1 = tmp1 * self.constants.N[r][v]
                    tmp2 = tmp2 * i
                    self.assertTrue(int(tmp1 / tmp2) * 1.0 == tmp1 / tmp2)
                    self.assertEqual(int(tmp1 / tmp2), M1[r][s][i - 1])
                    self.assertEqual(int(tmp1 / tmp2), M2[r][s][i - 1])
                    i = i + 1

    """
    Test that the eta s satisfy the properties in
    Theorem 6.4.3 (Carter p.95)
    note: test both eta and etars()
    """

    @unittest.skip("TODO")
    def test_etars(self):
        self.assertTrue(False)
