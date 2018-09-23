import numpy as np
import sympy as sp


class AdAction:
    type = 0

    constants = 0
    # underlying Lie algebra
    lie = 0
    # underlying group
    group = 0
    # underlying roots
    roots = 0

    def __init__(self, t):  # with_dist_expr=True,with_weyl_group=True):
        self.type = t

        self.lie = self.type.lie
        self.roots = self.type.rootsystem
        self.group = self.type.group
        self.constants = self.type.constants

    """
    Adjoint action of u_r on e_s
    --- Carter p.61
    --- !!! atentie, in Mrsi the index is shifted
    """

    def ad_ue(self, u, e):
        r = u[1]
        mr = self.roots.minus_r(r)
        if r == e[1]:
            return [[e[0], e[1], e[2]]]
        if mr == e[1]:
            result = [[e[0], e[1], e[2]]]
            result += [["h", r, u[2] * e[2]]]
            result += [["e", r, -u[2] ** 2 * e[2]]]
            return result

        s = e[1]
        Mrs = self.constants.M[r][s]
        rr = self.roots.roots[r]
        ss = self.roots.roots[s]
        result = [["e", e[1], e[2]]]
        for i in range(len(Mrs)):
            irs = self.roots.index((i + 1) * rr + ss)
            result.append(["e", irs, Mrs[i] * (u[2] ** (i + 1)) * e[2]])
        return result

    """
    Adjoint action of u_r on h_s
    --- Carter p.61
    """

    def ad_uh(self, u, h):
        r = u[1]
        s = h[1]
        Asr = self.constants.A[s][r]
        result = [["h", s, h[2]]]
        if Asr != 0:
            result += [["e", r, -Asr * u[2] * h[2]]]
        return result

    """
    Adjoint action u_r on atom
    """

    def ad_u_atom(self, u, a):
        if a[0] == "e":
            return self.ad_ue(u, a)
        if a[0] == "h":
            return self.ad_uh(u, a)

    """
    Adjoint action u_r on vector
    """

    def ad_u(self, u, v):
        result = []
        for a in v:
            result += self.ad_u_atom(u, a)
        return result

    """
    Adjoint action of h_r on e_s
    """

    def ad_he(self, t, e):
        r = t[1]
        s = e[1]
        Ars = self.constants.A[r][s]
        return [["e", s, (t[2] ** Ars) * e[2]]]

    """
    Adjoint action h_r on atom
    """

    def ad_h_atom(self, t, a):
        if a[0] == "e":
            return self.ad_he(t, a)
        if a[0] == "h":
            return [list(a)]

    """
    Adjoint action h_r on vector
    """

    def ad_h(self, t, v):
        result = []
        for a in v:
            result += self.ad_h_atom(t, a)
        return result

    """
    Adjoint action of t_r on e_s
    """

    def ad_te(self, t, e):
        r = t[1]
        s = e[1]
        if r == self.roots.minus_r(s):
            #
            # IS THIS OK FOR NEGATIVE ROOTS?
            #
            return [["e", s, t[2] * e[2]]]
        else:
            return [["e", s, e[2]]]

    """
    Adjoint action t_r on atom
    """

    def ad_t_atom(self, t, a):
        if a[0] == "e":
            return self.ad_te(t, a)
        if a[0] == "h":
            return [list(a)]

    """
    Adjoint action t_r on vector
    """
    # IS THIS OK FOR NEGATIVE ROOTS?
    #    def ad_t(self,t,v):
    #        result=[]
    #        for a in v:
    #            result+=self.ad_t_atom(t,a)
    #        return result

    """
    Adjoint action of n_r on e_s
    """

    def ad_ne(self, n, e):
        r = n[1]
        s = e[1]
        etars = self.constants.eta[r][s]
        wr = self.group.weyl.refs[r]
        pwr = self.group.weyl.perm[wr]
        result = [["e", pwr[s], etars * e[2]]]
        if n[2] != 0:
            result = self.ad_he(["h", r, n[2]], result[0])
        return result

    """
    Adjoint action of n_r on h_s
    """

    def ad_nh(self, n, h):
        r = n[1]
        s = h[1]
        wr = self.group.weyl.refs[r]
        pwr = self.group.weyl.perm[wr]
        return [["h", pwr[s], h[2]]]

    """
    Adjoint action n_r on atom
    """

    def ad_n_atom(self, n, a):
        if a[0] == "e":
            return self.ad_ne(n, a)
        if a[0] == "h":
            return self.ad_nh(n, a)

    """
    Adjoint action n_r on vector
    """

    def ad_n(self, n, v):
        result = []
        for a in v:
            result += self.ad_n_atom(n, a)
        return result

    """
    Adjoint action of an atoms u,t,n on vector v
    """

    def ad_atom(self, a, v):
        if a[0] == "u":
            return self.ad_u(a, v)
        if a[0] == "h":
            return self.ad_h(a, v)
        #        if a[0]=="t":
        #            return self.ad_t(a,v)
        if a[0] == "n":
            return self.ad_n(a, v)
        print("atom=", a)
        print("ad_action.ad_atom: this should not be!")
        return False

    """
    Adjoint action of the group element g on a vector v
    --- we act with the rightmost first so the list needs to be reverted
    """

    def ad(self, g, v):
        result = list(v)
        gg = list(g)
        gg.reverse()
        for x in gg:
            result = self.ad_atom(x, result)
        return result

    """
    simplify the expressions in a matrix
    --- for sympy symbols
    """

    @staticmethod
    def simplify_mat(mat):
        result = []
        for i in range(len(mat)):
            result.append([])
            for ee in mat[i]:
                e = sp.simplify(ee)
                if type(e) == float or \
                        type(e) == np.int64 or type(e) == sp.numbers.Float:
                    if int(e) * 1.0 == e:
                        # print(i,e)
                        e = int(e)
                result[i].append(e)
        return np.array(result)

    """
    Adjoint matrix for a group element g
    """

    def ad_mat(self, g):
        b = self.lie.basis
        result = []
        for x in b:
            tmp = self.ad(g, [x])
            # print(g,x,tmp)
            tmp = self.lie.canonic_list(tmp)
            result.append(tmp)
        result = np.array(result)
        result = result.transpose()
        return result
