import unittest
import sympy


class TestDeodhar(unittest.TestCase):
    deodhar = None
    ad = None

    unu = 1
    zero = 0

    """
    Check that the D cells and the DI cells are the same in their ad-action
    """

    # @unittest.skip("test")
    def check_D_DI(self, i, j):
        d = self.deodhar.cell_D(i, j)
        # print(d)
        di = self.deodhar.cell_DI(i, j)
        # print(di)
        tmp1 = self.ad.ad_mat(d)
        tmp2 = self.ad.ad_mat(di)
        tmp1 = tmp1 * self.unu
        tmp2 = tmp2 * self.unu
        # print(tmp1-tmp2)
        # print(self.unu)
        tmp1 = (tmp1 - tmp2) * self.unu
        tmp1 = tmp1 == self.zero
        self.assertTrue(False not in tmp1)

    def test_D_DI(self):
        de = self.deodhar.weyl.dist_expr
        # for i in range(len(de)):
        for i in range(len(de) - 3):
            for j in range(len(de[i][1])):
                print(">>>>>>", i, j, ">>", de[i][0][0], de[i][0][1], de[i][0][2])
                self.check_D_DI(i, j)
        #
        # ATENTIE LA FOR-ul asta
        #
        # print("SKIPPED last three!")

    """
    Check that:
    - the D cells and the UyiU cells are the same in their ad-action
    - UyiU is in Bruhat form
    - yi is correct
    """

    def check_D_UyiU(self, i, j):
        d = self.deodhar.cell_D(i, j)
        uyiu = self.deodhar.cell_UyiU(i, j)
        tmp1 = self.ad.ad_mat(d)
        de = self.deodhar.weyl.dist_expr[i]
        # z
        z = de[0][2]
        z = self.deodhar.group.w_to_n(self.deodhar.group.weyl.word(z))
        tmp2 = self.ad.ad_mat(z + uyiu)
        tmp1 = tmp1 * self.unu
        tmp2 = tmp2 * self.unu
        # print(tmp1-tmp2)
        tmp1 = (tmp1 - tmp2) * self.unu
        for a in range(len(tmp1)):
            for b in range(len(tmp1)):
                tmp1[a][b] = sympy.simplify(tmp1[a][b])
        tmp1 = tmp1 == self.zero
        self.assertTrue(False not in tmp1)
        self.assertTrue(self.deodhar.bruhat.in_Bruhat_form(uyiu, n_coef=-1))
        w = self.deodhar.group.e_to_w(uyiu)
        w = self.deodhar.weyl.wprod(w)
        w = self.deodhar.weyl.perm_inv[w]
        # print("w:",w)
        # print("dist_expr:",self.deodhar.weyl.dist_expr[i])
        self.assertEqual(w, self.deodhar.weyl.dist_expr[i][0][1])

    # @unittest.skip("test")
    def test_D_UyiU(self):
        de = self.deodhar.weyl.dist_expr
        # for i in range(len(de)):
        for i in range(len(de) - 3):
            for j in range(len(de[i][1])):
                print(">>>>>>", i, j, ">>", de[i][0][0], de[i][0][1], de[i][0][2])
                self.check_D_UyiU(i, j)
        #
        # ATENTIE LA FOR-ul asta
        #
        # print("SKIPPED last three!")

    """
    Check that:
    - the D cells and the Ux cells are the same in their ad-action
    - Ux is in Bruhat form
    - x is correct
    """

    def check_D_Ux(self, i, j):
        d = self.deodhar.cell_DI(i, j)
        ux = self.deodhar.cell_Ux(i, j)
        # print(ux)
        tmp1 = self.ad.ad_mat(d)
        tmp2 = self.ad.ad_mat(ux)
        tmp1 = tmp1 * self.unu
        tmp2 = tmp2 * self.unu
        # print(tmp1-tmp2)
        tmp1 = (tmp1 - tmp2) * self.unu
        for a in range(len(tmp1)):
            for b in range(len(tmp1)):
                tmp1[a][b] = sympy.simplify(tmp1[a][b])
        tmp1 = tmp1 == self.zero
        self.assertTrue(False not in tmp1)
        self.assertTrue(self.deodhar.bruhat.in_Bruhat_form(ux))
        w = self.deodhar.group.e_to_w(ux)
        w = self.deodhar.weyl.wprod(w)
        # w=self.deodhar.weyl.perm_inv[w]
        # self.assertEqual(w,self.deodhar.weyl.dist_expr[i][0][1])
        self.assertEqual(w, self.deodhar.weyl.dist_expr[i][0][0])

    # @unittest.skip("test")
    def test_D_Ux(self):
        de = self.deodhar.weyl.dist_expr
        # for i in range(len(de)):
        for i in range(len(de) - 3):
            for j in range(len(de[i][1])):
                print(">>>>>>", i, j, ">>", de[i][0][0], de[i][0][1], de[i][0][2])
                self.check_D_Ux(i, j)
        #
        # ATENTIE LA FOR-ul asta
        #
        # print("SKIPPED last three!")

    """
    Check that zUy is in UNT-form
    """

    @unittest.skip("test")
    def check_zUy_Bruhat_form(self, i, j):
        de = self.deodhar.weyl.dist_expr
        tmp = self.deodhar.cell_zUy(de[i][0][0], de[i][0][1], de[i][0][2], j)
        tmp = tmp[0] + tmp[1] + tmp[2]
        self.assertTrue(self.deodhar.bruhat.in_Bruhat_form(tmp, n_coef=-1, check_Uw=False))

    @unittest.skip("old")
    def test_zUy_Bruhat_form(self):
        de = self.deodhar.weyl.dist_expr
        for i in range(len(de)):
            for j in range(len(de[i][1])):
                # print(">>>>>>",i,j,">>",de[i][0][0],de[i][0][1],de[i][0][2])
                self.check_zUy_Bruhat_form(i, j)

    """
    Check that UxU is in strict Bruhat-form
    """

    def check_UxU_Bruhat_form(self, i, j):
        de = self.deodhar.weyl.dist_expr
        tmp = self.deodhar.cell_UxU(de[i][0][0], de[i][0][1], de[i][0][2], j)
        tmp = tmp[0] + tmp[1] + tmp[2]
        self.assertTrue(self.deodhar.bruhat.in_Bruhat_form(tmp))

    @unittest.skip("old")
    def test_UxU_Bruhat_form(self):
        de = self.deodhar.weyl.dist_expr
        for i in range(len(de)):
            for j in range(len(de[i][1])):
                # print(">>>>>>",i,j,">>",de[i][0][0],de[i][0][1],de[i][0][2])
                self.check_UxU_Bruhat_form(i, j)
