import sympy


class Curtis:
    type = 0

    # module for computing zUy and UxU
    deodhar = 0
    # Bruhat form
    bruhat = 0
    # the Chevalley group
    group = 0
    # the Weyl group
    weyl = 0

    # standard parabolics
    para = 0
    # distinguished expressions for standard parabolics
    dist_expr_p = 0
    # Deodhar cells
    D = 0
    # Deodhar cells DI-form
    DI = 0
    # Deodhar cells in zUyi form
    zUyi = 0
    # Deodhar cells in UxU form
    UxU = 0

    # the toral elements for the basis of the Hecke algebra of a GG-rep
    # given explicitly in derived classes
    tori = []
    # a second list of the same tori with "primed" variables
    tori2 = []
    # a third list of the same tori with "double primed" variables
    tori3 = []

    def __init__(self, t):
        self.type = t

        self.deodhar = self.type.deodhar
        self.bruhat = self.type.bruhat
        self.group = self.type.group
        self.weyl = self.type.weyl

        self.para = self.type.parabolics
        self.dist_expr_p = self.extract_para_dist_expr()
        # needs dist_expr_p:
        # self.load_cells()

    """
    Selecting those distinguished expressions corresponding to
    standard parabolic subgroups
    """

    def extract_para_dist_expr(self):
        de = self.weyl.dist_expr
        w0w = self.para.w0w
        result = []
        for i in range(len(de)):
            e = de[i]
            if e[0][0] in w0w and \
                    e[0][1] in w0w and \
                    e[0][2] in w0w:
                result.append(e + [i])
        return result

    """
    Select cells corresponding to dist_expr_p
    --- needs dist_expr_p
    """

    def load_cells(self):
        dep = self.dist_expr_p
        self.D = []
        self.DI = []
        self.zUyi = []
        self.UxU = []
        for e in dep:
            pos = e[len(e) - 1]
            tmpD = []
            tmpDI = []
            tmpzUyi = []
            tmpUxU = []
            for j in range(len(e[1])):
                # D and zUyi
                uyiu = self.deodhar.cell_UyiU(pos, j)
                tmpzUyi.append(uyiu)
                # DI and UxU
                uxu = self.deodhar.cell_Ux(pos, j)
                tmpUxU.append(uxu)
            self.D.append(tmpD)
            self.DI.append(tmpDI)
            self.zUyi.append(tmpzUyi)
            self.UxU.append(tmpUxU)

    """
    prepare the two forms of the cell
    """

    def prepare_zUy_UxU(self, ii, j):
        de = self.weyl.dist_expr
        x = de[ii][0][0]
        y = de[ii][0][1]
        z = de[ii][0][2]

        nx = self.group.w_to_n(self.weyl.word(x))
        ny = self.group.w_to_n(self.weyl.word(y))
        nz = self.group.w_to_n(self.weyl.word(z))

        ty = self.para.w0w.index(y)
        ty = self.tori2[ty]
        tyi = self.group.invert(ty)
        ytyi = self.group.conjugate_left(ny, tyi)

        tz = self.para.w0w.index(z)
        tz = self.tori3[tz]
        ztz = self.group.conjugate_left(nz, tz)

        uyiu = self.deodhar.cell_UyiU(ii, j)
        uxu = self.deodhar.cell_Ux(ii, j)

        uyiu = self.bruhat.split_strict_Bruhat(uyiu, n_coef=-1)

        ytyi0 = ytyi + self.group.invert(uyiu[2])

        uxu = self.bruhat.split_strict_Bruhat(uxu)

        uxu[0] = self.group.conjugate_left(ztz, uxu[0])
        ztzx = self.group.conjugate_right(ztz, nx)
        if nx != uxu[1]:
            print("curtis.prepare_zUy_UxU: this should not be!")
        uxu[3] = uxu[3] + self.group.invert(uyiu[3])
        uxu[3] = self.group.conjugate_right(uxu[3], ytyi0)
        uxu[2] = uxu[2] + ztzx + ytyi0

        uy = uyiu[0] + uyiu[1]
        uxu = uxu[0] + uxu[1] + self.group.canonic_th(uxu[2]) + self.group.canonic_u(uxu[3])

        for i in range(len(uy)):
            uy[i] = [uy[i][0], uy[i][1], sympy.simplify(uy[i][2])]

        for i in range(len(uxu)):
            uxu[i] = [uxu[i][0], uxu[i][1], sympy.simplify(uxu[i][2])]

        return [uy, uxu]

    """
    Get condition for toral elements to represent the same cell
    --- we need t0 in zUyi*t0
    --- we need t00 in Uxt00U

    [z*tz][U][(y*ty)^-1]t
    = [tz^(z^-1)][z][U][y^-1][(ty^-1)^(y^-1)]
    = [tz^(z^-1)][zUyi][t0^-1][(ty^-1)^(y^-1)]
    = [tz^(z^-1)][UxU][t0^-1][(ty^-1)^(y^-1)]
    = [tz^(z^-1)][U][x][t00][U][t0^-1][(ty^-1)^(y^-1)]
    """

    def structure_equation(self, i, j):
        x = self.dist_expr_p[i][0][0]
        y = self.dist_expr_p[i][0][1]
        z = self.dist_expr_p[i][0][2]
        # copiem ca sa nu modificam
        zUyi = [list(e) for e in self.zUyi[i][j]]
        UxU = [list(e) for e in self.UxU[i][j]]

        xx = self.weyl.word(x)
        xx = self.group.w_to_n(xx)
        yy = self.weyl.word(y)
        yy = self.group.w_to_n(yy)
        zz = self.weyl.word(z)
        zz = self.group.w_to_n(zz)

        #
        # toral part for y
        #
        # the order is important
        # this is the correct order to get t0 on the right
        t0 = yy + zUyi[1] + zUyi[2]
        t0 = self.group.canonic_nt(t0)
        if not self.group.all_t(t0):
            print("curtis.structure_equation: This should not be! (t0)")

        #
        # toral part for x
        #
        xxi = self.group.invert(xx)
        # the order is important
        # this is the correct order to get t0 on the right
        t00 = xxi + UxU[1] + UxU[2]
        t00 = self.group.canonic_nt(t00)
        if not self.group.all_t(t00):
            print("curtis.structure_equation: This should not be! (t00)")

        #
        # tz and ty
        #
        tz = self.para.w0w.index(z)
        # use the second set of variables for z
        tz = self.tori2[tz]
        ty = self.para.w0w.index(y)
        ty = self.tori[ty]
        # bring to other form
        # left U
        zztz = self.group.conjugate_left(zz, tz)
        UxU[0] = self.group.conjugate_left(zztz, UxU[0])
        xxizztz = self.group.conjugate_right(zztz, xxi)
        # right U
        t0i = self.group.invert(t0)
        UxU[3] = self.group.conjugate_right(UxU[3], t0i)
        tyi = self.group.invert(ty)
        yytyi = self.group.conjugate_left(yy, tyi)
        UxU[3] = self.group.conjugate_right(UxU[3], yytyi)

        tt = xxizztz + t00 + t0i + yytyi
        tt = self.group.canonic_t(tt)
        return [tt, zUyi, UxU]

    """
    Truncate the unipotent part
    and bring the two forms of the cells in the right form for
    the structure constants of the Hecke algebra of a GG-rep
    """

    def Hecke_GG_form(self, i, j):
        [tt, zUyi, UxU] = self.structure_equation(i, j)
        Uyz = self.group.truncate_u_sr(zUyi[0])
        #
        # just added !!! non-standard
        #
        # no Uyz=self.group.invert(Uyz)
        # no Uyz=self.group.canonic_u(Uyz)
        # no Uyz=self.group.truncate_u_sr(Uyz)

        Ux_left = self.group.truncate_u_sr(UxU[0])
        Ux_right = self.group.truncate_u_sr(UxU[3])

        Ux = Ux_left + Ux_right
        Ux = self.group.invert(Ux)
        Ux = self.group.canonic_u(Ux)
        Ux = self.group.truncate_u_sr(Ux)

        U = Ux + Uyz
        U = self.group.canonic_u(U)
        U = self.group.truncate_u_sr(U)

        return [tt, zUyi, UxU, U]

    """
    Produce a report for the j-th cell in the i-th case
    """

    def report(self, i, j):
        [uy, uxu] = self.prepare_zUy_UxU(i, j)
        uy = self.bruhat.split_strict_Bruhat(uy, n_coef=-1)
        uxu = self.bruhat.split_strict_Bruhat(uxu)

        de = self.weyl.dist_expr[i]
        word = self.weyl.word
        latex = self.group.latex
        truncate = self.group.truncate_u_sr
        print("############################")
        print("CASE: ", i, j)
        print("CONFIGURATION: ", de[0])
        print("DIST EXPR: ", de[1][j])
        print("------------------")
        print("Z: ", word(de[0][2]))
        print("Y: ", word(de[0][1]))
        print("X: ", word(de[0][0]))
        print("------------------")
        print("U in zUyi:")
        print("U1: ", latex(truncate(uy[0])))
        print("U in UxU:")
        print(uxu)
        print("U2: ", latex(truncate(uxu[0])))
        print("U3: ", latex(truncate(uxu[3])))
        print("------------------")
        print("Condition on toral element:")
        print("A) ", latex(uxu[2]))
        print("------------------")
        print("U to evaluate psi on:")
        Ux_left = truncate(uxu[0])
        Ux_right = truncate(uxu[3])
        Ux = Ux_left + Ux_right
        Ux = self.group.invert(Ux)
        Ux = self.group.canonic_u(Ux)
        Ux = truncate(Ux)

        U = Ux + uy[0]
        U = self.group.canonic_u(U)
        U = truncate(U)
        U = self.group.simplify_params(U)
        print(U)
        print(latex(U))
        print("############################")

    """
    Produce a report for the j-th cell in the i-th case
    """

    def report_file(self, i, j):
        f_name = "data/" + self.type.label + "/reports/" + str(i) + str(j) + ".rep"
        f_name = f_name.lower()
        f = open(f_name, "w")
        # [tt,zUyi,UxU,U]=self.Hecke_GG_form(i,j)
        [uy, uxu] = self.prepare_zUy_UxU(i, j)
        uy = self.bruhat.split_strict_Bruhat(uy, n_coef=-1)
        uxu = self.bruhat.split_strict_Bruhat(uxu)

        de = self.weyl.dist_expr[i]
        word = self.weyl.word
        latex = self.group.latex
        truncate = self.group.truncate_u_sr
        f.write("############################\n")
        f.write("CASE: " + str(i) + str(j) + "\n")
        f.write("CONFIGURATION: " + str(de[0]) + "\n")
        f.write("DIST EXPR: " + str(de[1][j]) + "\n")
        f.write("------------------")
        f.write("Z: " + str(word(de[0][2])) + "\n")
        # f.write("Y^-1t0: ",zUyi[1]+zUyi[2])
        f.write("Y: " + str(word(de[0][1])) + "\n")
        # f.write("Xt00: ",UxU[1]+UxU[2])
        f.write("X: " + str(word(de[0][0])) + "\n")
        f.write("------------------\n")
        f.write("U in zUyi:")
        f.write("U1: " + latex(truncate(uy[0])) + "\n")
        f.write("U2: " + latex(truncate(uxu[0])) + "\n")
        f.write("U in UxU:")
        f.write("U3: " + latex(truncate(uxu[3])) + "\n")
        f.write("------------------\n")
        f.write("Condition on toral element:\n")
        f.write("A) " + latex(uxu[2]) + "\n")
        f.write("------------------\n")
        f.write("U to evaluate psi on:\n")
        Ux_left = truncate(uxu[0])
        Ux_right = truncate(uxu[3])
        Ux = Ux_left + Ux_right
        Ux = self.group.invert(Ux)
        Ux = self.group.canonic_u(Ux)
        Ux = truncate(Ux)

        U = Ux + uy[0]
        U = self.group.canonic_u(U)
        U = truncate(U)
        U = self.group.simplify_params(U)
        f.write(latex(U) + "\n")
        f.write("############################\n")
        f.close()

    """
    Returns the index in the list dist_expr_p of the case c
    """

    def index(self, c):
        de = self.dist_expr_p
        tmp = [i[0] for i in de]
        return tmp.index(c)

    def latex_dist_expr(self, i, j):
        de = self.weyl.dist_expr[i][1][j]
        result = "$" + str([i + 1 for i in de[0]]) + "$"
        result += " (of type "
        t = ""
        vari = ""
        for k in range(len(de[0])):
            if k in de[1][0]:
                t += "A"
                vari += "$x_{" + str(k + 1) + "}\in k$, "
            elif k in de[1][1]:
                t += "B"
                vari += "$x_{" + str(k + 1) + "}\in k^{\\ast}$, "
            elif k in de[1][2]:
                t += "C"
                vari += "$x_{" + str(k + 1) + "}=1$, "
            else:
                print("curtis.latex_dist_expr: this should not be!")
                return
        result += t + ") " + vari
        return result

    """
    Produce a report for the j-th cell in the i-th case
    """

    def report_latex(self, i):
        ii = self.dist_expr_p[i][2]
        w0w = list(self.para.w0w)
        #
        # atentie inversez ultimul cu primul element aici
        #
        tmp = w0w[3]
        w0w[3] = w0w[2]
        w0w[2] = tmp

        case = [w0w.index(k) for k in self.dist_expr_p[i][0]]
        case_str = "".join([str(k) for k in case])
        fname = "latex/" + self.type.label + "/" + case_str + ".tex"
        f = open(fname, "w+")
        f.write("\subsection{" + case_str + "}\n")
        f.write("\label{" + case_str + "}\n")

        for j in range(len(self.dist_expr_p[i][1])):
            f.write(self.latex_dist_expr(ii, j) + ":\n")
            self.report_latex_sub(ii, j, f, [self.para.w0w.index(k) for k in self.dist_expr_p[i][0]])  # case)
        other_case = case_str[1] + case_str[0] + case_str[2]
        f.write("Should equal \eqref{" + other_case + "}\n")
        f.close()

    def report_latex_sub(self, i, j, f, case):
        # [tt,zUyi,UxU,U]=self.Hecke_GG_form(i,j)
        [uy, uxu] = self.prepare_zUy_UxU(i, j)
        uy = self.bruhat.split_strict_Bruhat(uy, n_coef=-1)
        uxu = self.bruhat.split_strict_Bruhat(uxu)

        latex = self.group.latex
        truncate = self.group.truncate_u_sr

        f.write("$$" + latex(self.tori[case[0]]) + "=" + latex(uxu[2]) + "$$\n")
        Ux_left = truncate(uxu[0])
        Ux_right = truncate(uxu[3])
        Ux = Ux_left + Ux_right
        Ux = self.group.invert(Ux)
        Ux = self.group.canonic_u(Ux)
        Ux = truncate(Ux)

        U = Ux + uy[0]
        U = self.group.canonic_u(U)
        U = truncate(U)
        U = self.group.simplify_params(U)
        f.write("$$\sum\psi(" + latex(U) + ")$$\n")

    def report_latex_files(self):
        w0w = list(self.para.w0w)

        # atentie inversez ultimul cu primul element aici
        #
        tmp = w0w[3]
        w0w[3] = w0w[2]
        w0w[2] = tmp
        result = []
        for i in range(len(self.dist_expr_p)):
            case = [w0w.index(k) for k in self.dist_expr_p[i][0]]
            case_str = "".join([str(k) for k in case])
            result.append("\\input{" + self.type.label + "/" + case_str + ".tex}\n")
        return result

    def report_latex_all(self):
        for i in range(len(self.dist_expr_p)):
            self.report_latex(i)

    def report_poly(self, ii, j):
        i = self.dist_expr_p[ii][2]
        [uy, uxu] = self.prepare_zUy_UxU(i, j)
        uy = self.bruhat.split_strict_Bruhat(uy, n_coef=-1)
        uxu = self.bruhat.split_strict_Bruhat(uxu)

        truncate = self.group.truncate_u_sr

        result = []
        result += [[self.tori[self.para.w0w.index(self.dist_expr_p[ii][0][0])], uxu[2]]]

        Ux_left = truncate(uxu[0])
        Ux_right = truncate(uxu[3])
        Ux = Ux_left + Ux_right
        Ux = self.group.invert(Ux)
        Ux = self.group.canonic_u(Ux)
        Ux = truncate(Ux)

        U = Ux + uy[0]
        U = self.group.canonic_u(U)
        U = truncate(U)
        U = self.group.simplify_params(U)

        poly = []
        for u in U:
            poly += [u[2]]

        result += [poly]
        return result
