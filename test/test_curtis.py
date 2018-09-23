import unittest


class TestCurtis(unittest.TestCase):
    curtis = None

    """
    Check that zUy is in UNT-form
    """

    # @unittest.skip("old")
    def test_zUy_Bruhat_form(self):
        de = self.curtis.dist_expr_p
        for i in range(len(de)):
            for j in range(len(de[i][1])):
                print(">>>>>>", i, j, ">>", de[i][0][0], de[i][0][1], de[i][0][2])
                [uy, uxu] = self.curtis.prepare_zUy_UxU(de[i][2], j)
                self.assertTrue(self.curtis.bruhat.in_Bruhat_form(uy, n_coef=-1, check_Uw=False))
                self.assertTrue(self.curtis.bruhat.in_Bruhat_form(uxu, n_coef=1, check_Uw=False))

    """
    Check that UxU is in strict Bruhat-form
    """

    @unittest.skip("old")
    def test_UxU_Bruhat_form(self):
        de = self.curtis.dist_expr_p
        for i in range(len(de)):
            for j in range(len(de[i][1])):
                print(">>>>>>", i, j, ">>", de[i][0][0], de[i][0][1], de[i][0][2])
                tmp = self.curtis.UxU[i][j]
                tmp = tmp[0] + tmp[1] + tmp[2]
                self.assertTrue(self.curtis.group.in_Bruhat_form(tmp))

    """
    Necessary condition for Hecke algebra of GG-rep to be abelian
    """

    @unittest.skip("old")
    def test_dist_sym(self):
        de = self.curtis.dist_expr_p
        tmp = [i[0] for i in de]
        for e in tmp:
            ee = [e[1], e[0], e[2]]
            self.assertTrue(ee in tmp)

    ############################################
    #
    # Observations
    #
    ############################################
    """
    Test that the toral part in zUy comes only from inverting y
    !!! Fails (for A2)
    ---> so we don't get exactly y^-1
    ---> the order also matters: yty or yyt
    ---> but both yty and yyt are standard toral
    """

    # @unittest.skip("test")
    @unittest.skip("old")
    def test_yi(self):
        de = self.curtis.dist_expr_p
        for i in range(len(de)):
            for j in range(len(de[i][1])):
                print(">>>>>>", i, j, ">>", de[i][0][0], de[i][0][1], de[i][0][2])
                tmp = self.curtis.zUyi[i][j]
                yt = tmp[1] + tmp[2]
                yt = self.curtis.group.canonic_nt(yt)
                y = de[i][0][1]
                y = self.curtis.weyl.word(y)
                y = self.curtis.group.w_to_n(y)
                yyt = self.curtis.group.canonic_nt(y + yt)
                self.assertTrue(self.curtis.group.all_t(yyt))
                yty = self.curtis.group.canonic_nt(yt + y)
                self.assertTrue(self.curtis.group.all_t(yty))
