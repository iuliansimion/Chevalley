import unittest


class TestGroup(unittest.TestCase):
    group = None
    variables = None

    """
    Test that inverting an atom twice is the identity map
    """

    @unittest.skip("test")
    def test_invert_atom(self):
        self.assertTrue(False)

    """
    Test that inverting a products of pairs of atoms twice is the identity map
    """

    @unittest.skip("test")
    def test_invert(self):
        self.assertTrue(False)

    """
    Check that conjugation by N commutes with taking the canonic form of T
    """

    def check_conjugate_nh_canonic(self, N, T):
        conj_nh_left = self.group.conjugate_nh_left
        canonic = self.group.canonic_h
        decompose_atom = self.group.decompose_h_atom
        # first conjugate with n
        tmp1 = conj_nh_left(N, T)
        tmp1 = canonic(tmp1)
        # first decompose
        tmp2 = [conj_nh_left(N, t)[0] for t in decompose_atom(T)]
        tmp2 = [t for i in tmp2 for t in decompose_atom(i)]
        tmp2 = canonic(tmp2)
        self.assertEqual(tmp1, tmp2)

    def test_conjugate_nh_canonic(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for n in range(nr):
            for t in range(nr):
                self.check_conjugate_nh_canonic(["n", n, xx[0]], ["h", t, xx[1]])

    """
    Check the decomposition [na]=[ua][u-a][u]
    """

    def check_conjugate_uu_nu(self, N, U):
        npr = self.group.roots.nr_pos_roots
        V = ["u", N[1], N[2]]
        if N[1] < npr:
            Vm = ["u", N[1] + npr, -N[2] ** -1]
        else:
            Vm = ["u", N[1] - npr, -N[2] ** -1]
        # conjugate with u
        tmp1 = self.group.conjugate_left([V], [U])
        tmp1 = self.group.conjugate_left([Vm], tmp1)
        tmp1 = self.group.conjugate_left([V], tmp1)
        tmp1 = self.group.canonic_u(tmp1)
        # conjugate with n
        tmp2 = self.group.conjugate_left([N], [U])
        self.assertEqual(tmp1, tmp2)

    def test_conjugate_uu_nu(self):
        nr = self.group.roots.nr_roots
        npr = self.group.roots.nr_pos_roots
        xx = self.variables.x
        for n in range(nr):
            for u in range(nr):
                if n not in [u, u + npr, u - npr]:
                    self.check_conjugate_uu_nu(["n", n, xx[0]], ["u", u, xx[1]])

    """
    Check that commutation on left and on right cancel for u
    """

    def check_commute_uu_left_right(self, U1, U2):
        tmp1 = self.group.commute_uu_left(U1, U2)
        tmp2 = self.group.commute_uu_right(U1, U2)
        self.assertEqual([], self.group.canonic_u(tmp1 + self.group.invert(tmp2)))

    def test_commute_uu_left_right(self):
        nr = self.group.roots.nr_roots
        npr = self.group.roots.nr_pos_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                if r not in [s, s + npr, s - npr]:
                    self.check_commute_uu_left_right(["u", r, xx[0]], ["u", s, xx[1]])

    """
    Check the decomposition of h to t
    """

    def check_decomposition_h_to_t(self, h, u):
        dec = self.group.simple_h_to_t_atom(h)
        tmp1 = self.group.conjugate_left([h], [u])
        tmp2 = self.group.conjugate_left(dec, [u])
        self.assertEqual(tmp1, tmp2)
        tmp1 = self.group.conjugate_right([u], [h])
        tmp2 = self.group.conjugate_right([u], dec)
        self.assertEqual(tmp1, tmp2)

    def test_decomposition_h_to_t(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        #
        # the t are only for simple roots
        # but when decomposing we get product of ts with s simple
        #
        for r in range(nr):
            for s in range(nr):
                # if not r in [s,s+npr,s-npr]:
                self.check_decomposition_h_to_t(["h", r, xx[0]], ["u", s, xx[1]])

    """
    Check that left and right conjugation by A1 cancels on A2
    """

    def check_conjugate_left_right(self, A1, A2):
        tmp = self.group.conjugate_left([A1], [A2])
        tmp = self.group.conjugate_right(tmp, [A1])
        if self.group.all_u(tmp):
            tmp = self.group.canonic_u(tmp)
            A2 = [A2]
        elif self.group.all_hn(tmp):
            tmp = self.group.canonic_nh(tmp)
            # the same is true for A2
            A2 = self.group.canonic_nh([A2])
        else:
            A2 = [A2]
        self.assertEqual(A2, tmp)

    def test_conjugate_uu_left_right(self):
        nr = self.group.roots.nr_roots
        npr = self.group.roots.nr_pos_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                if r not in [s, s + npr, s - npr]:
                    self.check_conjugate_left_right(["u", r, xx[0]], ["u", s, xx[1]])

    def test_conjugate_hu_left_right(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                self.check_conjugate_left_right(["h", r, xx[0]], ["u", s, xx[1]])

    def test_conjugate_tu_left_right(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        #
        # the t are only for simple roots
        #
        for r in range(self.group.type.rank):
            for s in range(nr):
                self.check_conjugate_left_right(["t", r, xx[0]], ["u", s, xx[1]])

    def test_conjugate_nu_left_right(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                self.check_conjugate_left_right(["n", r, xx[0]], ["u", s, xx[1]])

    def test_conjugate_nh_left_right(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                self.check_conjugate_left_right(["n", r, xx[0]], ["h", s, xx[1]])

    def test_conjugate_nn_left_right(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                self.check_conjugate_left_right(["n", r, xx[0]], ["n", s, xx[1]])

    """
    Check that conjugation is `associative', i.e. that it is a group action
    """

    def check_conjugate_left_left(self, A, B, C):
        conjugate_left = self.group.conjugate_left
        tmp1 = conjugate_left([B], [C])
        tmp1 = conjugate_left([A], tmp1)
        tmp2 = conjugate_left(conjugate_left([A], [B]), conjugate_left([A], [C]))
        self.assertEqual(tmp1, tmp2)

    def test_conjugate_nhu_left_left(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                for v in range(nr):
                    self.check_conjugate_left_left(["n", r, xx[0]], ["h", v, xx[1]], ["u", s, xx[2]])

    """
    Check the action of n(t) and h(t)n(1) on u
    """

    def check_nh_decomposition(self, N, U):
        conjugate_left = self.group.conjugate_left
        canonic_n = self.group.canonic_n
        tmp1 = conjugate_left([N], [U])
        T = ["h", N[1], N[2]]
        N1 = ["n", N[1], 1]
        tmp2 = conjugate_left([N1], [U])
        tmp2 = conjugate_left([T], tmp2)
        tmp3 = conjugate_left(canonic_n([N]), [U])
        self.assertEqual(tmp1, tmp2)
        self.assertEqual(tmp2, tmp3)

    # @unittest.skip("test")
    def test_nh_decomposition(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                self.check_nh_decomposition(["n", r, xx[0]], ["u", s, xx[2]])

    """
    Check that move_atom_left/right for unipotent elements
    - here we move u_{r} left and then right of a 
    generic U(pos_roots) times a generic T
    """

    def check_move_atom_u_left_right(self, r):
        U = self.group.generic_u_pos()
        T = self.group.generic_h()
        u = ["u", 0, self.group.var.y[0]]
        tmp = self.group.move_atom_left(U + T, u)
        tmp = self.group.move_atom_right(tmp[0], tmp[1:len(tmp)])
        uu = tmp[len(tmp) - 1]
        self.assertEqual(u, uu)

    def test_move_atom_u_left_right(self):
        npr = self.group.roots.nr_pos_roots
        for r in range(npr):
            self.check_move_atom_u_left_right(r)

    ############################################
    #
    # OBSERVATIONS
    #
    ############################################
    """
    Check action of n_{r}(t) and n_{-r}(-t^-1), are they equal?
    Apparently so, works for A2, B2 and G2:
    TO PROVE: this
    """

    # @unittest.skip("test")
    def check_neg_root_n(self, N, U):
        npr = self.group.roots.nr_pos_roots
        tmp1 = self.group.conjugate_left([N], [U])
        if N[1] < npr:
            nr = N[1] + npr
        else:
            nr = N[1] - npr
        term = -N[2] ** -1
        if type(term) == float:
            if int(term) * 1.0 == term:
                term = int(term)
        Nn = ["n", nr, term]
        tmp2 = self.group.conjugate_left([Nn], [U])
        self.assertTrue(tmp1, tmp2)

    def test_neg_root_n(self):
        nr = self.group.roots.nr_roots
        xx = self.variables.x
        for r in range(nr):
            for s in range(nr):
                self.check_neg_root_n(["n", r, xx[0]], ["u", s, xx[2]])
