import numpy as np


class LieAlgebra:
    type = 0

    # constants
    constants = 0
    # underlying roots
    roots = 0
    # Chevalley Basis
    basis = 0
    # we also need variables
    var = 0

    def __init__(self, t):
        self.type = t

        self.roots = self.type.rootsystem
        self.constants = self.type.constants
        self.basis = []
        for r in range(self.roots.nr_roots):
            self.basis.append(["e", r, 1])
        for r in range(self.type.rank):
            self.basis.append(["h", r, 1])

        self.var = self.type.var

    """
    Remove the elements in the vector v with coefficients 0
    """

    @staticmethod
    def clean(v):
        result = []
        for a in v:
            if a[2] != 0:
                result += [a]
        return result

    """
    For an h=2r/(r,r) atom, decompose it as a sum of h_s for s simple
    
    !!! requires that the simple roots are in order [1,0,.. [0,1,.. [0,0,1..
    """

    def canonic_h_atom(self, h):
        npr = self.roots.nr_pos_roots
        r = h[1]
        if r >= npr:
            return self.canonic_h_atom(["h", r - npr, -h[2]])
        # now r is positive
        result = []
        rr = self.roots.roots[r]
        short = self.roots.short
        ls = self.roots.ls_square
        if short[r]:
            for s in range(self.type.rank):
                if rr[s] != 0:
                    if short[s]:
                        result.append(["h", s, rr[s] * h[2]])
                    else:
                        result.append(["h", s, rr[s] * ls * h[2]])
        else:
            for s in range(self.type.rank):
                if rr[s] != 0:
                    if short[s]:
                        term = ls ** (-1)
                        # print(term)
                        if type(term) == float:
                            if int(term) * 1.0 == term:
                                term = int(term)
                        term = term * h[2]
                        # print(type(term))
                        if type(term) == float:
                            if int(term) * 1.0 == term:
                                term = int(term)
                        # print(term)
                        term = rr[s] * term
                        # print(type(term))
                        #
                        # implement this in other places as well!!!
                        #
                        if type(term) == float or type(term) == np.float64:
                            if int(term) * 1.0 == term:
                                term = int(term)
                        # print(term)
                        result.append(["h", s, term])
                    else:
                        result.append(["h", s, rr[s] * h[2]])
        return result

    """
    Replaces all the h in v by their canonic form
    """

    def canonic_h(self, v):
        result = []
        for a in v:
            if a[0] == "h":
                result += self.canonic_h_atom(a)
            else:
                result += [list(a)]
        return result

    """
    Convert list of coefficients to vector
    """

    def l_to_v(self, lis):
        nr = self.roots.nr_roots
        result = []
        for i in range(nr):
            if lis[i] != 0:
                result += [["e", i, lis[i]]]
        for i in range(self.type.rank):
            if lis[nr + i] != 0:
                result += [["h", i, lis[nr + i]]]
        return result

    """
    Add similar elements of the basis together
    --- return list of coefficients
    """

    def canonic_list(self, v):
        v = list(v)
        v = self.clean(v)
        v = self.canonic_h(v)
        nr = self.roots.nr_roots
        result = [0] * (nr + self.type.rank)
        for a in v:
            if a[0] == "e":
                result[a[1]] += a[2]
            elif a[0] == "h":
                result[a[1] + nr] += a[2]
            else:
                print("lie_algebra.canonic: this should not be!")
                return False

        return result

    """
    Add similar elements of the basis together
    --- return vector
    """

    def canonic(self, v):
        result = self.canonic_list(v)
        return self.l_to_v(result)

    """
    Lie bracket for atoms
    """

    def ad_atom(self, a, b):
        A = self.constants.A
        N = self.constants.N
        rsum = self.roots.sum
        if a[0] == "h":
            if b[0] == "h":
                return []
            if b[0] == "e":
                return [[b[0], b[1], A[a[1]][b[1]] * b[2] * a[2]]]
            else:
                print("chv_lie_algebra.ad_atom: this should not be!")
        if a[0] == "e":
            if b[0] == "h":
                return [[a[0], a[1], -A[b[1]][a[1]] * b[2] * a[2]]]
            if b[0] == "e":
                mr = self.roots.minus_r(a[1])
                if b[1] == mr:
                    return [["h", a[1], b[2] * a[2]]]
                else:
                    rs = rsum[a[1]][b[1]]
                    if rs != -1:
                        return [["e", rs, N[a[1]][b[1]] * b[2] * a[2]]]
                    else:
                        return []
            else:
                print("chv_lie_algebra.ad_atom: this should not be!")

    """
    Lie bracket atom on vector
    """

    def ad_atom_v(self, a, v):
        result = []
        for i in v:
            result = result + self.ad_atom(a, i)
        return result

    """
    Lie bracket element on vector
    """

    def ad(self, e, v):
        result = []
        for i in e:
            for j in v:
                result = result + self.ad_atom(i, j)
        return result

    """
    Adjoint matrix for a lie algebra element e
    """

    def ad_mat(self, e):
        b = self.basis
        result = []
        for x in b:
            tmp = self.ad(e, [x])
            # print(g,x,tmp)
            tmp = self.canonic_list(tmp)
            result.append(tmp)
        result = np.array(result)
        result = result.transpose()
        return result

    """
    e_{\alpha_1}+e_{\alpha_2}+...+e_{\alpha_rang}
    """

    def regular_rep(self):
        result = []
        for i in range(self.type.rank):
            result.append(["e", i, 1])
        return result

    ############################################
    #
    # Generic elements
    #
    ############################################
    def generic_e_pos(self, x=0):
        if x == 0:
            x = self.var.x
        npr = self.roots.nr_pos_roots
        result = []
        for i in range(npr):
            result.append(["e", i, x[i]])
        return result

    def generic_e_neg(self, x=0):
        if x == 0:
            x = self.var.x
        npr = self.roots.nr_pos_roots
        result = []
        for i in range(npr):
            result.append(["e", i + npr, x[i + npr]])
        return result

    def generic_h(self, x=0):
        if x == 0:
            x = self.var.x
        r = self.type.rank
        npr = self.roots.nr_pos_roots
        result = []
        for i in range(r):
            result.append(["h", i, x[i + 2 * npr]])
        return result

    def generic_element(self, x=0):
        if x == 0:
            x = self.var.x
        return self.generic_e_pos(x) + self.generic_e_neg(x) + self.generic_h(x)
