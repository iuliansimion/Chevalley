import unittest


class TestAd(unittest.TestCase):
    ad = None

    unu = 1
    zero = 0

    """
    Check that LEFT conjugation by A1 cancels on A2
    """

    def check_conjugate_left(self, A1, A2):
        tmp1 = self.ad.group.conjugate_left([A1], [A2])
        tmp1 = self.ad.ad_mat(tmp1)
        tmp2 = [A1, A2] + self.ad.group.invert([A1])
        tmp2 = self.ad.ad_mat(tmp2)
        tmp1 = tmp1 * self.unu
        tmp2 = tmp2 * self.unu
        #
        # tmp1=tmp1==tmp2
        #
        tmp1 = (tmp1 - tmp2) * self.unu
        tmp1 = tmp1 == self.zero
        self.assertTrue(False not in tmp1)

    """
    Check that RIGHT conjugation by A1 cancels on A2
    """

    def check_conjugate_right(self, A1, A2):
        tmp1 = self.ad.group.conjugate_right([A1], [A2])
        tmp1 = self.ad.ad_mat(tmp1)
        tmp2 = self.ad.group.invert([A2]) + [A1, A2]
        tmp2 = self.ad.ad_mat(tmp2)
        tmp1 = tmp1 * self.unu
        tmp2 = tmp2 * self.unu
        #
        # tmp1=tmp1==tmp2
        #
        tmp1 = (tmp1 - tmp2) * self.unu
        tmp1 = tmp1 == self.zero
        self.assertTrue(False not in tmp1)

    """
    TODO: Remove the condition r!=-s
    """

    # @unittest.skip("test")
    def test_conjugate_uu_left(self):
        nr = self.ad.group.roots.nr_roots
        npr = self.ad.group.roots.nr_pos_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                if r not in [s, s + npr, s - npr]:
                    # print(">>>>>>",r,s)
                    self.check_conjugate_left(["u", r, xx[0]], ["u", s, xx[1]])

    """
    TODO: Remove the condition r!=-s
    """

    # @unittest.skip("test")
    def test_conjugate_uu_right(self):
        nr = self.ad.group.roots.nr_roots
        npr = self.ad.group.roots.nr_pos_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                if r not in [s, s + npr, s - npr]:
                    # print(">>>>>>",r,s)
                    self.check_conjugate_right(["u", r, xx[0]], ["u", s, xx[1]])

    # @unittest.skip("test")
    def test_conjugate_hu_left(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_left(["h", r, xx[0]], ["u", s, xx[1]])

    # @unittest.skip("test")
    def test_conjugate_uh_right(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_right(["u", r, xx[0]], ["h", s, xx[1]])

    # @unittest.skip("test!!")
    def test_conjugate_nu_left(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_left(["n", r, xx[0]], ["u", s, xx[1]])

    # @unittest.skip("test!!")
    def test_conjugate_un_right(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_right(["u", r, xx[0]], ["n", s, xx[1]])

    # @unittest.skip("test!!")
    def test_conjugate_nh_left(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_left(["n", r, xx[0]], ["h", s, xx[1]])

    # @unittest.skip("test!!")
    def test_conjugate_nh_right(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_right(["h", r, xx[0]], ["n", s, xx[1]])

    # @unittest.skip("test!!")
    def test_conjugate_nn_left(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_left(["n", r, xx[0]], ["n", s, xx[1]])

    def test_conjugate_nn_right(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            for s in range(nr):
                # print(">>>>>>",r,s)
                self.check_conjugate_right(["n", r, xx[0]], ["n", s, xx[1]])

    """
    Check decomposition:
    [["u",r,term]]=[["u",mr,term**-1],["n",mr,-term**-1],["u",mr,term**-1]]
    """

    def check_u_rmr(self, u):
        tmp1 = self.ad.ad_mat([u])
        mr = self.ad.roots.minus_r(u[1])
        term = u[2]
        uu = [["u", mr, term ** -1], ["n", mr, -term ** -1], ["u", mr, term ** -1]]
        tmp2 = self.ad.ad_mat(uu)
        tmp1 = tmp1 * self.unu
        tmp2 = tmp2 * self.unu
        tmp1 = tmp1 == tmp2
        self.assertTrue(False not in tmp1)

    """
    TODO: Remove the condition r!=-s
    """

    # @unittest.skip("test")
    def test_u_rmr(self):
        nr = self.ad.group.roots.nr_roots
        xx = self.ad.group.var.x
        for r in range(nr):
            self.check_u_rmr(["u", r, xx[0]])
