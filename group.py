import sympy


class Group:
    type = 0

    #    tip="none"
    #    rang=0
    #    label="none"

    constants = 0
    # underlying roots
    roots = 0
    # underlying weyl group
    weyl = 0
    # we also need variables
    var = 0

    def __init__(self, t, with_dist_expr=True, with_weyl_group=True):

        self.type = t

        self.roots = t.rootsystem
        if with_weyl_group:
            self.weyl = t.weyl
        if with_dist_expr:
            self.type.loadDistinguishedExpressions()
        self.constants = t.constants
        self.var = t.var

    """
    Invert an element of type "u" or "h" or "n"
    --- inverse of u: Carter p.87
    --- inverse of t: Carter p.92
    --- inverse of n: Carter p.96 + Carter p.102 Theorem 7.2.2 (see also p.103)
    TEST this: especially type n
    """

    @staticmethod
    def invert_atom(E):
        if E[0] == "u":
            return [[E[0], E[1], -E[2]]]

        if E[0] in ["t", "h"]:
            arg = E[2] ** -1
            if type(arg) == float:
                if int(arg) * 1.0 == arg:
                    arg = int(arg)
            return [[E[0], E[1], arg]]
        # if E[0]=="n" and E[2]==1:
        #    return [[E[0],E[1],-1]]
        if E[0] == "n" and E[2] == -1:
            return [[E[0], E[1], 1]]
        if E[0] == "n":
            arg = E[2] ** -1
            if type(arg) == float:
                if int(arg) * 1.0 == arg:
                    arg = int(arg)
            return [[E[0], E[1], 1], ["h", E[1], -arg]]
        print("invert_element: not implemented")
        return False

    """
    Invert a product of atoms
    """

    def invert(self, LE):
        # we don't want the initial list to be reversed
        LE = list(LE)
        LE.reverse()
        result = []
        for i in LE:
            for j in self.invert_atom(i):
                result.append(j)
        return result

    ############################################
    #
    # U U
    #
    ############################################
    """
    For unipotent atoms, calculate a*b*a^-1*b^-1
    eg. commutator_uu(["u",1,xx[1]],["u",2,xx[2]])
    returns list of atoms
    """

    def commutator_uu(self, S, R):
        if S[1] == self.roots.minus_r(R[1]):
            print("group.commutator_uu: case not implemented!")
            return

        C = self.constants.C
        index = self.roots.index
        roots = self.roots.roots
        comm = []
        for c in C[R[1]][S[1]]:
            irjs = c[0] * roots[R[1]] + c[1] * roots[S[1]]
            irjs = index(irjs)
            if irjs == -1:
                print("group.commutator_uu: This should not be")
            comm.append(["u", irjs, c[2] * ((R[2]) ** c[0]) * ((-S[2]) ** c[1])])
        return comm

    """
    Calculate [us][ur]=[ur][us][[us],[ur]]
    --- see Carter p. 76 for details
    --- with commutator on the RIGHT

    eg. commute_uu(["u",1,xx[1]],["u",2,xx[2]])
    TODO: test comm should equal
    invert_elem(commutator_uu(invert_elem(S),invert_elem(R)))
    """

    def commute_uu_right(self, S, R):
        if S[1] == self.roots.minus_r(R[1]):
            print("group.commute_uu_right: case not implemented!")
            return

        C = self.constants.C
        index = self.roots.index
        roots = self.roots.roots
        comm = []
        for c in C[R[1]][S[1]]:
            irjs = c[0] * roots[R[1]] + c[1] * roots[S[1]]
            irjs = index(irjs)
            comm.append(["u", irjs, c[2] * ((-R[2]) ** c[0]) * (S[2] ** c[1])])
        # comm could equal []
        return [R, S] + comm

    """
    Calculate [us][ur]=[ur][us][[us],[ur]]
    --- see Carter p. 76 for details
    --- with commutator on the LEFT
    """

    def commute_uu_left(self, S, R):
        if S[1] == self.roots.minus_r(R[1]):
            print("group.commute_uu_left: case not implemented!")
            return

        comm = self.commutator_uu(S, R)
        # comm could equal []
        return comm + [R, S]

    """
    Calculate [u1][u2][u1]^-1=...
    --- using element on the left

    eg. conjugate_uu_left(["u",1,xx[1]],["u",2,xx[2]])
    """

    def conjugate_uu_left(self, U1, U2):
        if U1[1] == self.roots.minus_r(U2[1]):
            print("group.conjugate_uu_left: case not implemented!")
            return

        mu2 = self.roots.minus_r(U2[1])
        if U1[1] == mu2:
            # print("grou.conjugate_uu_left: not implemented +/- r!")
            term = U2[2]
            term = term ** -1
            if type(term) == float:
                if int(term) * 1.0 == term:
                    term = int(term)
            return [["u", U1[1], term + U1[2]], ["n", U1[1], -term], ["u", U1[1], term - U1[2]]]
        comm = self.commutator_uu([U2[0], U2[1], -U2[2]], U1)
        return [U2] + comm

    """
    Calculate [u2]^-1[u1][u2]=...
    --- using element on the right

    eg. commute.conjugate_uu_right(["u",1,xx[1]],["u",2,xx[2]])
    """

    def conjugate_uu_right(self, U1, U2):
        if U1[1] == self.roots.minus_r(U2[1]):
            print("group.conjugate_uu_right: case not implemented!")
            return

        mu2 = self.roots.minus_r(U2[1])
        if U1[1] == mu2:
            # print("grou.conjugate_uu_left: not implemented +/- r!")
            term = U1[2]
            term = term ** -1
            if type(term) == float:
                if int(term) * 1.0 == term:
                    term = int(term)
            return [["u", U2[1], term - U2[2]], ["n", U2[1], -term], ["u", U2[1], term + U2[2]]]
        comm = self.commutator_uu([U1[0], U1[1], -U1[2]], [U2[0], U2[1], -U2[2]])
        return [U1] + comm

    """
    Remove
    --- Us with parameter 0
    --- []
    """

    @staticmethod
    def clean_u(LU):
        result = []
        for u in LU:
            if u:
                if u[0] == "u":
                    if u[2] != 0:
                        result.append(u)
                else:
                    result.append(u)
        return result

    """
    Reduce list of u s by multiplying elements corresponding
    to the same root if they are consecutive

    eg. simplify_u([["u",1,xx[1]],["u",1,xx[2]]])
    TEST this
    """

    def simplify_u(self, LU):
        # we don't want the initial list to be reversed
        LU = list(LU)
        LU = self.clean_u(LU)
        modif = True
        result = []
        while modif:
            modif = False
            result = []
            i = 0
            while i < len(LU) - 1:
                if LU[i][0] != "u":
                    result.append(LU[i])
                else:
                    tmp = LU[i][2]
                    j = i + 1
                    while j < len(LU) and \
                            LU[i][0] == "u" and \
                            LU[j][0] == "u" and \
                            LU[j][1] == LU[i][1]:
                        tmp = tmp + LU[j][2]
                        j = j + 1
                    result.append(["u", LU[i][1], tmp])
                    if j > i + 1:
                        modif = True
                        i = j - 1
                i = i + 1
            # don't forget the last element if it was not merged
            if i < len(LU):
                result.append(LU[len(LU) - 1])
            result = self.clean_u(result)
            LU = result
        return result

    """
    Test if all atoms are unipotent
    """

    @staticmethod
    def all_u(L):
        for e in L:
            if e[0] != "u":
                return False
        return True

    """
    Test if all atoms are unipotent for positive roots
    """

    def all_up(self, L):
        npr = self.roots.nr_pos_roots
        for e in L:
            if e[0] != "u":
                return False
            if e[1] >= npr:
                return False
        return True

    """
    Bring product of unipotent atoms in canonic form,
    i.e. ordered increasingly on the order of roots

    eg. canonic_u([["u",2,xx[1]],["u",1,xx[0]]])
    """

    def canonic_u(self, LU):
        # we don't want the initial list to be reversed
        LU = list(LU)
        if not LU:
            return []

        LU = self.simplify_u(LU)

        if not self.all_u(LU):
            print("group.canonic_u: not all elements are unipotent!")
            return LU

        if not self.roots.in_half_plane([u[1] for u in LU]):
            print("group.canonic_u: elements are not all in one half-plane!")
            return LU

        result = []
        modif = True
        while modif:
            modif = False
            result = []
            i = 0
            # go up to len(LU)-2
            while i < len(LU) - 1:
                if LU[i][1] > LU[i + 1][1]:
                    result = result + self.commute_uu_right(LU[i], LU[i + 1])
                    i = i + 2
                    modif = True
                else:
                    result.append(LU[i])
                    i = i + 1
            # don't forget the last element
            if i == len(LU) - 1:
                result.append(LU[len(LU) - 1])
            result = self.simplify_u(result)
            LU = result
        return result

    """
    Bring positive-u-subexpressions in canonic form
    """

    def canonic_pos_u_part(self, L):
        L = list(L)
        result = []
        i = 0
        npr = self.roots.nr_pos_roots
        while i < len(L):
            # skip non-u and negative u
            while i < len(L) and (L[i][0] != "u" or L[i][1] >= npr):
                result.append(L[i])
                i = i + 1
            j = i + 1
            while j < len(L) and L[j][0] == "u" and L[j][1] < npr:
                j = j + 1
            # print(">>>>>>",i,j,L[i:j])
            result = result + self.canonic_u(L[i:j])
            i = j
        return result

    ############################################
    #
    # T T
    #
    ############################################
    """
    Remove
    --- Ts with parameter 1
    --- []
    """

    @staticmethod
    def clean_h(L):
        result = []
        for t in L:
            if t:
                if t[0] == "h":
                    if t[2] != 1:
                        result.append(t)
                else:
                    result.append(t)
        return result

    """
    Remove
    --- Ts with parameter 1
    --- []
    """

    @staticmethod
    def clean_t(L):
        result = []
        for t in L:
            if t:
                if t[0] == "t":
                    if t[2] != 1:
                        result.append(t)
                else:
                    result.append(t)
        return result

    """
    Test if all atoms are standard toral
    """

    @staticmethod
    def all_h(L):
        for e in L:
            if e[0] != "h":
                return False
        return True

    """
    Test if all atoms are standard toral
    """

    @staticmethod
    def all_t(L):
        for e in L:
            if e[0] != "t":
                return False
        return True

    """
    Test if all atoms are standard toral
    """

    @staticmethod
    def all_th(L):
        for e in L:
            if e[0] != "t" and e[0] != "h":
                return False
        return True

    """
    Reduce list of t s by multiplying elements corresponding
    to the same root if they are consecutive
    --- only if all elements are standard toral
    --- it works differently from group.simplify_u

    eg. simplify_h([["h",1,xx[1]],["h",1,xx[2]]])
    TEST this
    """

    def simplify_h(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        L = self.clean_h(L)

        if not self.all_h(L):
            print("group.simplify_h: not all elements are standard h-toral")
            return L

        # select the roots which appear:
        rr = list(set([i[1] for i in L]))
        result = [["h", i, 1] for i in rr]

        for t in L:
            i = rr.index(t[1])
            term = t[2]
            if type(term) is float:
                if int(term) * 1.00 == term:
                    term = int(term)
            result[i][2] = result[i][2] * term

        return self.clean_h(result)

    def simplify_t(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        L = self.clean_t(L)

        if not self.all_t(L):
            print("group.simplify_t: not all elements are standard t-toral")
            return L

        # select the roots which appear:
        rr = list(set([i[1] for i in L]))
        result = [["t", i, 1] for i in rr]

        for t in L:
            i = rr.index(t[1])
            term = t[2]
            if type(term) is float:
                if int(term) * 1.00 == term:
                    term = int(term)
            result[i][2] = result[i][2] * term

        return self.clean_t(result)

    """
    Decompose h_r(t) as product of h_s(tt) with s a simple root

    TEST: needs test that the simple roots are the first elements in
    roots, ordered as [1,0,0,.. [0,1,0.. [0,0,1..
    """

    def decompose_h_atom(self, t):
        if t[0] != "h":
            print("group.decompose_h_atom: only for t atoms!")
            return t

        ls = self.roots.ls_square
        short = self.roots.short
        rang = self.type.rank

        r = t[1]
        R = self.roots.roots[r]
        t = t[2]

        result = []
        for i in range(rang):
            if R[i] != 0:
                if short[r]:
                    if short[i]:
                        # print("a")
                        e = R[i]
                        if int(e) * 1.00 == e:
                            e = int(e)
                        term = t ** e
                    else:
                        # print("b")
                        e = R[i] * ls
                        if int(e) * 1.00 == e:
                            e = int(e)
                        term = t ** e
                else:
                    if short[i]:
                        # print("c")
                        e = R[i] * (ls ** -1)
                        if int(e) * 1.00 == e:
                            e = int(e)
                        term = t ** e
                    else:
                        # print("d")
                        e = R[i]
                        if int(e) * 1.00 == e:
                            e = int(e)
                        term = t ** e

                if type(term) == float:
                    if int(term) * 1.00 == term:
                        term = int(term)
                result.append(["h", i, term])
        return result

    """
    convert h_r with r simple to product of t_s for s simple
    !!! Assumes that simple roots are of the form [1,0,.. [0,1,..
    """

    def simple_h_to_t_atom(self, h):
        A = self.constants.A
        Ar = A[h[1]]
        result = []
        for s in range(self.type.rank):
            if Ar[s] != 0:
                result += [["t", s, h[2] ** Ar[s]]]
        return result

    """
    convert h_r to product of t_s for s simple
    """

    def h_to_t_atom(self, h):
        hh = self.decompose_h_atom(h)
        result = []
        for x in hh:
            result += self.simple_h_to_t_atom(x)
        return result

    """
    Product of toral elements in L calculated in canonic form
    i.e. as product of toral elements corresponding to simple roots
    --- uses group.decompose_h_atom()

    eg. canonic_h([["h",1,xx[1]],["h",2,xx[2]]])
    """

    def canonic_h(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        L = self.clean_h(L)

        if not self.all_h(L):
            print("group.canonic_h: not all elements are standard toral")
            return L

        result = []
        for t in L:
            result = result + self.decompose_h_atom(t)

        result = self.simplify_h(result)
        return result

    """
    converts product of t and h to a canonical product of t
    """

    def canonic_th(self, L):
        if not self.all_th(L):
            print("group.canonic_th: not all atoms are t or h!")
            return

        result = []
        LL = list(L)
        for a in LL:
            if a[0] == "h":
                result += self.h_to_t_atom(a)
            else:
                result += [a]

        return self.simplify_t(result)

    ############################################
    #
    # T U
    #
    ############################################
    """
    aba^-1b^-1
    eg. commutator_hu(["h",1,3],["u",0,1])
    """

    def commutator_hu(self, T, U):
        A = self.constants.A
        return [["u", U[1], U[2] * (T[2] ** A[T[1]][U[1]] - U[2])]]

    """
    [ta][ub][ta^-1][ub^-1]=[ub] 
    --- commutator on the left

    eg. commute_hu_left(["h",1,xx[1]],["u",2,xx[2]])
    """

    def commute_hu_left(self, T, U):
        comm = self.commutator_hu(T, U)
        return comm + [U, T]

    """
    [ha][ub][ha^-1]=[ub]
    ---  conjugate with LEFT atom

    eg. conjugate_hu_left(["h",1,xx[1]],["u",2,xx[2]])
    """

    def conjugate_hu_left(self, T, U):
        A = self.constants.A
        term = (T[2] ** A[T[1]][U[1]])
        if type(term) is float:
            if int(term) * 1.0 == term:
                term = int(term)
        return [["u", U[1], U[2] * term]]

    """
    [ta][ub][ta^-1]=[ub]
    ---  conjugate with LEFT atom

    eg. conjugate_tu_left(["h",1,xx[1]],["u",2,xx[2]])
    !!! Assumes that the simple roots are [1,0,.. , [0,1,..
    """

    def conjugate_tu_left(self, T, U):
        r = self.roots.roots[U[1]]
        return [["u", U[1], (T[2] ** int(r[T[1]])) * U[2]]]

    #        if T[1]==self.roots.minus_r(U[1]):
    #            return [["u",U[1],(T[2]**(-1))*U[2]]]
    #            #print("group.conjugate_tu_left: case not implemented!")
    #            #return
    #
    #        if T[1]==U[1]:
    #            return [["u",U[1],T[2]*U[2]]]#
    #
    #        return [["u",U[1],U[2]]]

    """
    [ha^-1][ub][ha]=[ub]
    ---  conjugate with RIGHT atom

    eg. conjugate_uh_right(["n",1,xx[1]],["h",2,xx[2]])
    """

    def conjugate_uh_right(self, U, T):
        TT = self.invert_atom(T)
        return self.conjugate_hu_left(TT[0], U)

    """
    [ta^-1][ub][ta]=[ub]
    ---  conjugate with RIGHT atom

    eg. conjugate_ut_right(["u",1,xx[1]],["t",2,xx[2]])
    """

    def conjugate_ut_right(self, U, T):
        TT = self.invert_atom(T)
        return self.conjugate_tu_left(TT[0], U)

    ############################################
    #
    # N N
    #
    ############################################
    """
    Split n(t) as h(t)n
    Bring n_b(t) to form n_ah(s) where a is a positive root, so b=+/-a
    SEE Observation tests, there is something to prove
    """

    def canonic_n_atom(self, N):
        npr = self.roots.nr_pos_roots
        if N[0] != "n":
            print("group.canonic_n_atom: this is for atomic n!")
            return N
        # if N[2]==1:
        #    return [N]
        # return [["h",N[1],N[2]],["n",N[1],1]]
        if N[2] == 1 and N[1] < npr:
            return [N]
        if N[1] < npr:
            return [["h", N[1], N[2]], ["n", N[1], 1]]
        # here N[1]>npr
        npr = self.roots.nr_pos_roots
        rn = N[1] - npr
        term = -N[2] ** -1
        if type(term) == float:
            if int(term) * 1.0 == term:
                term = int(term)
        Nn = ["n", rn, term]
        return self.canonic_n_atom(Nn)

    """
    Bring a product of n s in the form n1n2..t1t2..
    --- is the same as group.simplify_nh()
    """

    def simplify_n(self, L):
        return self.simplify_nh(L)

    """
    Given w in the Weyl group give product of n(1) in the group
    """

    @staticmethod
    def w_to_n(w):
        result = []
        for i in w:
            result.append(["n", i, 1])
        return result

    """
    Extract the word w in the Weyl group from
    the subword E in G given by n s
    """

    @staticmethod
    def e_to_w(E):
        result = []
        for i in E:
            if i[0] == "n":
                result.append(i[1])
        return result

    """
    Given w in the Weyl group give product of n(-1) in the group
    """

    @staticmethod
    def w_to_ni(w):
        result = []
        for i in w:
            result.append(["n", i, -1])
        return result

    """
    Cancel of the form na(1)na(-1) or na(-1)na(1)
    --- there can also be other elements
    """

    @staticmethod
    def cancel_n(L):
        result = []
        i = 0
        while i < len(L) - 1:
            if L[i][0] == "n" and \
                    L[i + 1][0] == "n" and \
                    L[i][1] == L[i + 1][1] and \
                    L[i][2] == 1 and \
                    L[i + 1][2] == -1:
                i = i + 2
            elif L[i][0] == "n" and \
                    L[i + 1][0] == "n" and \
                    L[i][1] == L[i + 1][1] and \
                    L[i][2] == -1 and \
                    L[i + 1][2] == 1:
                i = i + 2
            else:
                result.append(L[i])
                i = i + 1
        # don't forget the last element in case it didn't cancel
        if i < len(L):  # and L[i][0]=="n": it can be different from n
            result.append(L[i])

        return result

    """
    cancel_n and
    replace n_{-r}(t) with n_{r}(-t^-1) for positive roots r
    """

    def clean_n(self, L):
        L = list(L)
        L = self.cancel_n(L)
        npr = self.roots.nr_pos_roots
        i = 0
        while i < len(L):
            if L[i][0] == "n" and L[i][1] >= npr:
                term = -L[i][2] ** -1
                if type(term) == float:
                    if int(term) * 1.0 == term:
                        term = int(term)
                L[i] = [L[i][0], L[i][1] - npr, term]
            i = i + 1
        return L

    """
    [na(q)][nb(t)][na(q)^-1]=[hc(eta)][nc(s)] where c=n(b) 
    --- using atom on the LEFT
    --- na(q)=ha(q)na(1) (Carter p.96)
    --- !! I cannot use canonic_nh() here because that method uses this one

    eg. conjugate_nn_left(["n",1,xx[1]],["n",2,xx[2]])
    TEST: that conjugate_nn_left(["n",1,1],["n",1,1])==["n",1,1]
    --- it can be in ANOTHER form !!! like [['t', 4, -1], ['n', 4, 1]]
    """

    def conjugate_nn_left(self, Nr, Ns):
        perm = self.weyl.perm
        refs = self.weyl.refs
        eta = self.constants.eta

        if Nr[2] == 1 and Ns[2] == 1:
            ns = perm[refs[Nr[1]]][Ns[1]]
            return [["h", ns, eta[Nr[1]][Ns[1]]], ["n", ns, 1]]

        # split Nr and Ns:
        Hr = ["h", Nr[1], Nr[2]]
        Nr = ["n", Nr[1], 1]

        Hs = ["h", Ns[1], Ns[2]]
        Ns = ["n", Ns[1], 1]

        Hrs = self.conjugate_nh_left(Nr, Hs)

        Nrs = self.conjugate_nn_left(Nr, Ns)

        Nrs[1] = self.conjugate_hn_left(Hr, Nrs[1])[0]

        result = self.simplify_h([Hrs[0], Nrs[0]])
        # the toral part can end up being 1
        result = result + [Nrs[1]]

        #
        # se remarks at beginning of method
        #
        # result=self.simplify_nh(result) doesn't help much

        return result

    """
    n^-1wn
    --- use atom on the RIGHT
    --- na(q)=ha(q)na(1) (Carter p.96)
    --- TEST: requires that self.invert_atom(Ns) has at most two elements

    xr=["n",1,xx[1]]
    xs=["n",2,xx[2]]
    a=invert_elem(xs) a
    b=conjugate_hu_left(a[1],xr) b

    eg. conjugate_nn_right(["n",1,xx[1]],["n",2,xx[2]])
    """

    def conjugate_nn_right(self, Nr, Ns):
        if Ns[0] != "n":
            print("group.conjugate_nn_right: this methods is for atoms")
            return Ns

        NNs = self.invert_atom(Ns)

        if len(NNs) == 2:
            Nr = self.conjugate_hn_left(NNs[1], Nr)[0]
        Nr = self.conjugate_nn_left(NNs[0], Nr)

        return Nr

    """
    Bring product of n s and h s in the form wh
    with w corresponding to an element in lWord
    --- the same as group.canonic_nh()
    --- parameter is a list of n

    w^-1*n=h  =>  n=wh .. so I determine h and return!!!: wh
    """

    def canonic_n(self, L):
        return self.canonic_nh(L)

    """
    Same as canonic_n but write n_r(-1)=n_r^-1 instead of n_r
    """

    def canonic_ni(self, L):
        return self.canonic_nih(L)

    ############################################
    #
    # N T H
    #
    ############################################
    """
    Test if all atoms are h s or n s
    """

    @staticmethod
    def all_hn(L):
        for e in L:
            if e[0] != "h" and e[0] != "n":
                return False
        return True

    """
    Test if all atoms are standard toral or n s
    """

    @staticmethod
    def all_thn(L):
        for e in L:
            if e[0] != "h" and e[0] != "n" and e[0] != "t":
                return False
        return True

    """
    aba^-1b^-1
    """

    def commutator_nh(self, N, T):
        perm = self.weyl.perm
        refs = self.weyl.refs
        c = perm[refs[N[1]]][T[1]]
        term = T[2] ** -1
        if type(term) == float:
            if int(term) * 1.0 == term:
                term = int(term)
        return self.canonic_h([["h", c, T[2]], ["h", T[1], term]])

    """
    [na][tb][na^-1][tb^-1]=[tc][tb^-1] where c=na(b)
    --- commutator on the left

    eg. commute_nt(["n",1,xx[1]],["h",2,xx[2]])
    """

    def commute_nh_left(self, N, T):
        comm = self.commutator_nh(N, T)
        return comm + [T, N]

    """
    [na][tb][na^-1]=[tc] where c=na(b)
    --- us the atom on the LEFT

    eg. conjugate_nh_left(["n",1,xx[1]],["h",2,xx[2]])
    """

    def conjugate_nh_left(self, N, T):
        perm = self.weyl.perm
        refs = self.weyl.refs
        c = perm[refs[N[1]]][T[1]]
        return [["h", c, T[2]]]

    """
    [na][tb][na^-1]
    !!! b is simple
      =[tb] if b!=a
      =[tb][hb] if a=b NO!!! this is only when a is simple
    --- us the atom on the LEFT

    eg. conjugate_nt_left(["n",1,xx[1]],["h",2,xx[2]])
    """

    def conjugate_nt_left(self, N, T):
        r = self.roots.roots[N[1]]
        p = int(-r[T[1]])
        if p != 0:
            return [list(T), ["h", N[1], T[2] ** p]]

        return [list(T)]

    """
    [na^-1][tb][na]=[tc] where c=na(b)
    --- us the atom on the RIGHT
    --- TEST: needs that group.inverse of atom n puts n on left most position

    eg. conjugate_hn_right(["h",1,xx[1]],["n",2,xx[2]])
    """

    def conjugate_hn_right(self, T, N):
        NN = self.invert_atom(N)
        return self.conjugate_nh_left(NN[0], T)

    """
    [na^-1][tb][na]=[tc] where c=na(b)
    --- us the atom on the RIGHT
    --- TEST: needs that group.inverse of atom n puts n on left most position

    eg. conjugate_tn_right(["h",1,xx[1]],["n",2,xx[2]])
    """

    def conjugate_tn_right(self, T, N):
        NN = self.invert_atom(N)
        return self.conjugate_nt_left(NN[0], T)

    """
    Conjugate n by h, hnh^-1
    --- use the LEFT atom
    --- it conjugates like if n was u

    eg. conjugate_hn_left(["h",1,-1],["n",0,-1])
    TEST: needs testing 
    OBSERVATION: conjugating ["n",0,-1] by ["h",1,-1] inverts the n
    """

    def conjugate_hn_left(self, T, N):
        U = ["u", N[1], N[2]]
        U = self.conjugate_hu_left(T, U)
        return [["n", N[1], U[0][2]]]

    """
    Conjugate n by h, h^-1nh
    --- using the atom on the RIGHT
    --- it conjugates like if n was u

    eg. conjugate_hn_left(["n",0,-1],["h",1,-1])
    OBSERVATION: conjugating ["n",0,-1] by ["h",1,-1] inverts the n
    """

    def conjugate_nh_right(self, N, T):
        TT = self.invert_atom(T)
        return self.conjugate_hn_left(TT[0], N)

    """
    Conjugate n by t, tnt^-1
    --- use the LEFT atom
    --- it conjugates like if n was u

    eg. conjugate_tn_left(["t",1,-1],["n",0,-1])
    TEST: needs testing 
    """

    def conjugate_tn_left(self, T, N):
        U = ["u", N[1], N[2]]
        U = self.conjugate_tu_left(T, U)
        return [["n", N[1], U[0][2]]]

    """
    Conjugate n by t, t^-1nt
    --- using the atom on the RIGHT
    --- it conjugates like if n was u

    eg. conjugate_tn_left(["n",0,-1],["t",1,-1])
    """

    def conjugate_nt_right(self, N, T):
        TT = self.invert_atom(T)
        return self.conjugate_tn_left(TT[0], N)

    """
    Bring a product of n s and h s in the form n1n2..simplify_h(h1h2..)
    """

    def simplify_nh(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        L = self.clean_h(L)

        # remove na(1)na(-1) or na(-1)na(1)
        # L=self.cancel_n(L)
        L = self.clean_n(L)

        if not self.all_hn(L):
            print("group.simplify_nh: elements are not all t or n!")
            return L

        # split the n's
        result = []
        for e in L:
            if e[0] == "n":
                result += self.canonic_n_atom(e)
            else:
                result.append(e)

        modif = True
        while modif:
            modif = False
            for i in range(len(result) - 1):
                if result[i][0] == "h" and result[i + 1][0] == "n":
                    a = self.conjugate_hn_right(result[i], result[i + 1])[0]
                    result[i] = result[i + 1]
                    result[i + 1] = a
                    modif = True
        # now all h ar on the right from n s

        # find the last n in list
        i = 0
        lung = len(result)
        while i < lung and result[i][0] == "n":
            i = i + 1

        # bring the h-part in canonic form
        # trying to cancel has no effect because elements we have n(1)s
        result = result[0:i] + self.canonic_h(result[i:len(result)])
        return result

    """
    Bring product of n s and h s in the form wh
    with w corresponding to an element in lWord

    w^-1*n=h  =>  n=wh .. so I determine h and return!!!: wh
    NO if n contains reflections in non-simple roots, this doesn't work
    SO:
    n^-1*w=h  =>  n=wh^-1
    ALSO:
    w is a product of simple refs:
    ---> so I have to convert n_{-r} in something positive
    ---> I replace such elements by w0 n_{-r} w0^-1

    (YES see closure paper)CHECK IF n_{-r}(t)=n_{r}(-t^-1)
    !!!! DOESN't WORK IF L contains Nr with negative r
    uses simplify which converts Nr with negative r to Nr with positive r
    """

    def canonic_nh(self, L):

        if not L:
            return []

        # we don't want the initial list to be reversed
        L = list(L)

        if not self.all_hn(L):
            print(L)
            print("canonic_nh: couldn't find canonical from!")
            return L

        L = self.simplify_nh(L)
        # separate the n s (which are at the beginning now
        i = 0
        lung = len(L)
        while i < lung and L[i][0] == "n":
            i = i + 1
        LN = L[0:i]
        LH = L[i:len(L)]

        # refs=self.weyl.refs
        # find the reduced word corresponding to the n's
        w = [i[1] for i in LN]
        w = self.weyl.wprod(w)
        w = self.weyl.word(w)
        w = self.w_to_n(w)
        # hold w for later use,
        # else the invert(..) will invert ww as well
        ww = [i for i in w]
        #
        # se remark at beginning of method
        # w=self.invert(w)
        #
        LN = self.invert(LN)
        h = LN + w
        h = self.simplify_nh(h)

        modif = True
        while modif:
            modif = False
            # do I need this:?
            i = 0
            while i < len(h) - 1 and h[i][0] == "h":
                i = i + 1
            # here h[i] is the first n
            pair = False
            if i < len(h) - 1 and h[i + 1][0] == "n" and h[i + 1][1] == h[i][1]:
                pair = True
            while i < len(h) - 1 and not pair:
                if h[i + 1][0] == "n":
                    if h[i + 1][1] != h[i][1]:
                        # conjugate by the one on the RIGHT !!!
                        # in order to reach something that cancels
                        a = self.conjugate_nn_right(h[i], h[i + 1])
                        h = h[0:i] + [h[i + 1]] + a + h[i + 2:len(h)]
                        # I don't know if there is a t or not in a
                        i = i + 1
                        while h[i][0] == "h":
                            i = i + 1
                        # now i is again pointing to my n
                    else:
                        pair = True
                # h[i+1] is t
                else:
                    a = self.conjugate_nh_left(h[i], h[i + 1])
                    h = h[0:i] + a + [h[i]] + h[i + 2:len(h)]
                    i = i + 1

            if i < len(h) - 1 and h[i][0] != "n":
                print("canonic_n: this should not be!", i, len(h), h[i], h)
            if i < len(h) - 1 and h[i + 1][0] == "n" and h[i + 1][1] == h[i][1]:
                if h[i][2] == 1 and h[i + 1][2] == 1:
                    h = h[0:i] + [["h", h[i][1], -1]] + h[i + 2:len(h)]
                else:  # and repeat the loop
                    if h[i][2] != 1:
                        h = h[0:i] + [["h", h[i][1], h[i][2]], ["n", h[i][1], 1]] + h[i + 1:len(h)]
                        i = i + 1  # asta trebuie pt al doilea
                    i = i + 1
                    if h[i][2] != 1:
                        h = h[0:i] + [["h", h[i][1], h[i][2]], ["n", h[i][1], 1]] + h[i + 1:len(h)]
                modif = True
        #
        # se remark at beginning of method
        #
        # print(h)
        h = self.invert(h)
        # don't forget about LH
        h = self.canonic_h(h + LH)
        return ww + h

    """
    Bring the consecutive products of n and t in canonic forms
    and leave the u as they are
    """

    def canonic_nh_consecutive(self, L):
        result = []
        i = 0
        while i < len(L):
            if L[i][0] == "u":
                result.append(L[i])
                i = i + 1
            else:
                i0 = i
                while i < len(L) and L[i][0] != "u":
                    i = i + 1
                result = result + self.canonic_nh(L[i0:i])
        return result

    """
    Same as canonic_n but write n_r(-1)=n_r^-1 instead of n_r
    """

    def canonic_nih(self, N):
        # we don't want the initial list to be reversed
        N = list(N)
        N = self.canonic_n(N)
        i = 0
        NN = []
        # invert n s
        # suntem in forma canonica: NH
        while i < len(N) and N[i][0] == "n":
            NN = NN + [["h", N[i][1], -1], ["n", N[i][1], -1]]
            i = i + 1
        # append h s
        while i < len(N):
            NN.append(N[i])
            i = i + 1

        # conjugate n's to the left
        modif = True
        while modif:
            modif = False
            i = 0
            while i < len(NN) - 1:
                if NN[i][0] == "h" and NN[i + 1][0] == "n":
                    a = self.conjugate_hn_right(NN[i], NN[i + 1])[0]
                    NN[i] = NN[i + 1]
                    NN[i + 1] = a
                    modif = True
                i = i + 1
        i = 0
        while i < len(NN) and NN[i][0] == "n":
            i = i + 1
        H = NN[i:len(NN)]
        NN = NN[0:i]
        return NN + self.canonic_h(H)

    ############################################
    #
    # N U
    #
    ############################################

    """
    [na(q)][ub(t)][na(q)^-1]=[uc(s)] where c=n(b) 
    --- commutator on the LEFT

    na(q)=ha(q)na(1) (Carter p.96)
    
    eg. conjugate_nu_left(["n",1,xx[1]],["u",2,xx[2]])
    TEST: !!
    """

    def commute_nu_left(self, N, U):
        perm = self.weyl.perm
        refs = self.weyl.refs
        eta = self.constants.eta

        ns = perm[refs[N[1]]][U[1]]
        conj = ["u", ns, eta[N[1]][U[1]] * U[2]]
        conj = self.conjugate_hu_left(["h", N[1], N[2]], conj)
        # if ns==U[1]:
        #    comm=["u",conj[1],conj[2]-U[2]]
        #    if comm[2]==0:
        #        return [U,N]
        return conj + [N]

    """
    nun^-1 ... [na(q)][ub(t)][na(q)^-1]=[uc(s)] where c=n(b) 
    --- use atom on the LEFT
    --- na(q)=ha(q)na(1) (Carter p.96)

    eg. commute.conjugate_nu_left(["n",1,xx[1]],["u",2,xx[2]])
    """

    def conjugate_nu_left(self, N, U):
        perm = self.weyl.perm
        refs = self.weyl.refs
        eta = self.constants.eta

        ns = perm[refs[N[1]]][U[1]]
        result = ["u", ns, eta[N[1]][U[1]] * U[2]]
        result = self.conjugate_hu_left(["h", N[1], N[2]], result)
        return result

    """
    n^-1un .. na(q)^-1ub(t)na(q)=uc(s) where c=n(b) 
    --- use atom on the RIGHT
    --- na(q)=ha(q)na(1) (Carter p.96)

    eg. commute.conjugate_nu_right(["n",1,xx[1]],["u",2,xx[2]])
    """

    def conjugate_un_right(self, U, N):
        # we don't want the initial list to be reversed
        U = list(U)

        NN = self.invert_atom(N)
        result = [U]
        NN.reverse()
        for i in NN:
            if i[0] == "n":
                result = self.conjugate_nu_left(i, result[0])
            elif i[0] == "h":
                result = self.conjugate_hu_left(i, result[0])
            else:
                print("conjugate_nu_right: this should not be!")
                return False
        return result

    ############################################
    #
    # U T N
    #
    ############################################
    """
    Conjugate ABA^-1 for atoms A and B
    """

    def conjugate_left_pair(self, A, B):
        if A[0] == "u":
            if B[0] == "u":
                return self.conjugate_uu_left(A, B)
            # if B[0]=="h":
            # if B[0]=="n":
        if A[0] == "h":
            if B[0] == "u":
                return self.conjugate_hu_left(A, B)
            if B[0] in ["h", "t"]:
                return [B]
            # if B[0]=="n":
        if A[0] == "t":
            if B[0] == "u":
                return self.conjugate_tu_left(A, B)
            if B[0] in ["h", "t"]:
                return [B]
            # if B[0]=="n":
        if A[0] == "n":
            if B[0] == "u":
                return self.conjugate_nu_left(A, B)
            if B[0] == "h":
                return self.conjugate_nh_left(A, B)
            if B[0] == "t":
                return self.conjugate_nt_left(A, B)
            if B[0] == "n":
                return self.conjugate_nn_left(A, B)

        print("group.conjugate_left_pair: not implemented! for ", A, B, "\n")
        return False

    """
    Conjugate the elements in L by A: ALA^-1
    --- A is on LEFT
    --- A is an atom
    """

    def conjugate_left_atom(self, A, L):
        result = []
        for i in L:
            result += self.conjugate_left_pair(A, i)
        return result

    """
    Conjugate the elements in L by the elements in A: ALA^-1
    --- A is on LEFT
    --- A is not an atom
    --- needs all conjugation configurations
    TEST: not tested
    """

    def conjugate_left(self, A, L):
        # we don't want the initial list to be reversed
        A = list(A)
        result = list(L)
        A.reverse()
        for a in A:
            result = self.conjugate_left_atom(a, result)
        return result

    """
    Conjugate A^-1BA for atoms A and B
    """

    def conjugate_right_pair(self, B, A):
        if A[0] == "u":
            if B[0] == "u":
                return self.conjugate_uu_right(B, A)
            #
            # if B[0]=="h":
            # return self.conjugate_hu_right(B,A)
            # if B[0]=="n":
        if A[0] == "h":
            if B[0] == "u":
                return self.conjugate_uh_right(B, A)
            if B[0] == "h":
                return [B]
            # if B[0]=="n":
        if A[0] == "t":
            if B[0] == "u":
                return self.conjugate_ut_right(B, A)
            if B[0] in ["h", "t"]:
                return [B]
            # if B[0]=="n":
        if A[0] == "n":
            if B[0] == "u":
                return self.conjugate_un_right(B, A)
            if B[0] == "h":
                return self.conjugate_hn_right(B, A)
            if B[0] == "t":
                return self.conjugate_tn_right(B, A)
            if B[0] == "n":
                return self.conjugate_nn_right(B, A)

        print("group.conjugate_right_pair: not implemented!\n")
        return False

    """
    Conjugate the elements in L by A: A^-1LA
    --- A is on RIGHT
    --- A is an atom
    """

    def conjugate_right_atom(self, L, A):
        result = []
        for i in L:
            result += self.conjugate_right_pair(i, A)
        return result

    """
    Conjugate the elements in L by the elements in A: A^-1LA
    --- A is on RIGHT
    --- A is not an atom
    --- needs all conjugation configurations
    TEST: not tested
    """

    def conjugate_right(self, L, A):
        # we don't want the initial list to be reversed
        A = list(A)
        result = list(L)
        # !!! the list doesn't have to be reverted
        # A.reverse()
        for a in A:
            result = self.conjugate_right_atom(result, a)
        return result

    """
    Move all the u s on the LEFT
    """

    def move_u_left(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        modif = True
        while modif:
            modif = False
            for i in range(len(L) - 1):
                if L[i][0] != "u" and L[i + 1][0] == "u":
                    a = self.conjugate_left([L[i]], [L[i + 1]])[0]
                    L[i + 1] = L[i]
                    L[i] = a
                    modif = True
        return L

    """
    Move all the u s on the RIGHT
    """

    def move_u_right(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        modif = True
        while modif:
            modif = False
            for i in range(len(L) - 1):
                if L[i][0] == "u" and L[i + 1][0] != "u":
                    a = self.conjugate_right([L[i]], [L[i + 1]])[0]
                    L[i + 1] = L[i]
                    L[i] = a
                    modif = True
        return L

    """
    Move all the h s on the LEFT
    """

    def move_h_left(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        modif = True
        while modif:
            modif = False
            for i in range(len(L) - 1):
                if L[i][0] == "u" and L[i + 1][0] == "h":
                    a = self.conjugate_right([L[i]], [L[i + 1]])[0]
                    L[i] = L[i + 1]
                    L[i + 1] = a
                    modif = True
                if L[i][0] == "n" and L[i + 1][0] == "h":
                    a = self.conjugate_left([L[i]], [L[i + 1]])[0]
                    L[i + 1] = L[i]
                    L[i] = a
                    modif = True
        return L

    """
    Move all the h s on the RIGHT
    """

    def move_h_right(self, L):
        # we don't want the initial list to be reversed
        L = list(L)
        modif = True
        while modif:
            modif = False
            for i in range(len(L) - 1):
                if L[i][0] == "h" and L[i + 1][0] == "u":
                    a = self.conjugate_left([L[i]], [L[i + 1]])[0]
                    L[i + 1] = L[i]
                    L[i] = a
                    modif = True
                if L[i][0] == "h" and L[i + 1][0] == "n":
                    a = self.conjugate_right([L[i]], [L[i + 1]])[0]
                    L[i] = L[i + 1]
                    L[i + 1] = a
                    modif = True
        return L

    """
    Successively conjugate atom A by elements in E
    moving it from RIGHT to LEFT 
    while leaving commutator on the spot
    the result is [A']+L'
    --- mostly for unipotent elements
    --- TEST that after left_conjugation the main part stays on LEFT
    """

    def move_atom_left(self, E, A):
        # we don't want the initial list to be reversed
        E = list(E)
        A = list(A)
        E.reverse()
        result = []
        for e in E:
            tmp = self.conjugate_left([e], [A])
            result = tmp[1:len(tmp)] + [e] + result
            A = tmp[0]
        return [A] + result

    """
    Successively conjugate atom A by elements in E
    moving it from LEFT to RIGHT
    while leaving commutator on the spot
    the result is L'+[A']
    --- mostly for unipotent elements
    --- TEST that after left_conjugation the main part stays on LEFT
    --- --- doesn't happen for unipotents
    """

    def move_atom_right(self, A, E):
        result = []
        for e in E:
            tmp = self.conjugate_right([A], [e])
            tmp2 = self.move_atom_right(tmp[0], tmp[1:len(tmp)])
            result = result + [e] + tmp2[0:len(tmp2) - 1]
            A = tmp[0]
        return result + [A]

    """
    Truncate the expression of a product of u to those elements
    ["u",s,*] with s a simple root
    """

    def truncate_u_sr(self, U):
        if not self.all_u(U):
            print("group.truncate_u_sr: not all atoms are unipotent")
            return U
        rang = self.type.rank
        result = []
        for u in U:
            if u[1] < rang:
                result.append(u)
        return result

    ############################################
    #
    # Generic elements
    #
    ############################################
    def generic_u_pos(self, x=0):
        if x == 0:
            x = self.var.x
        npr = self.roots.nr_pos_roots
        result = []
        for i in range(npr):
            result.append(["u", i, x[i]])
        return result

    def generic_h(self, x=0):
        if x == 0:
            x = self.var.z
        result = []
        for i in range(self.type.rank):
            result.append(["h", i, x[i]])
        return result

    def generic_t(self, x=0):
        if x == 0:
            x = self.var.z
        result = []
        for i in range(self.type.rank):
            result.append(["t", i, x[i]])
        return result

    """
    UwU=U1s1.U2s2...U
    Return a parametrization of Ux=U1s1.U2s2...
    """

    def generic_Uw(self, w, xx=0):
        w = self.weyl.word(w)
        if xx == 0:
            xx = self.var.x
        result = []
        for r in range(len(w)):
            result.append(["u", w[r], xx[r]])
            result.append(["n", w[r], 1])

        return result

    """
    U[w]U=U_w[w]U
    Return a parametrization of U_w
    !!! U_w is on the left, so I have to use plen_set(w^-1)
    """

    def generic_Uw2(self, w, xx=0):
        wi = self.weyl.perm_inv[w]
        rr = self.weyl.plen_set(wi)
        if xx == 0:
            xx = self.var.x
        result = []
        for i in range(len(rr)):
            result.append(["u", rr[i], xx[rr[i]]])
        return result

    """
    U[w]U=U_w[w]U
    Return a parametrization of U_w[w] (also w)
    in the form U1s1U2s2...
    !!! U_w is on the left, so I have to use plen_set(w^-1)
    """

    def generic_Uww(self, w, xx=0):
        wi = self.weyl.perm_inv[w]
        wiw = self.weyl.word(wi)
        # rr=self.weyl.plen_set(wi)
        if xx == 0:
            xx = self.var.x
        result = []
        for i in range(len(wiw)):
            result.append(["u", wiw[i], xx[i]])
            result.append(["n", wiw[i], -1])
        return result

    @staticmethod
    def simplify_params(e):
        for i in range(len(e)):
            e[i] = [e[i][0], e[i][1], sympy.simplify(e[i][2])]
        return e

    ############################################
    #
    # Latex
    #
    ############################################
    def latex_atom(self, a):
        if not a:
            return "1"
        # result=a[0]+"_{"+str(a[1]+1)+"}("+str(a[2])+")"
        result = a[0] + "_{" + str(a[1] + 1) + "}(" + self.var.latex(a[2]) + ")"
        return result

    def latex(self, L):
        if not L:
            return "1"
        result = ""
        for a in L:
            result = result + self.latex_atom(a)
        return result
